use crate::{
    commands::Command,
    utils::{input_to_reader, output_to_writer, table_path_to_writer},
};
use seq_io::{fasta::Record, parallel::parallel_fasta};

#[derive(clap::ArgEnum, Clone, Debug, PartialEq)]
pub enum Strand {
    Forward,
    Reverse,
    Both,
}

#[derive(serde::Serialize, Debug)]
struct Row {
    orf_id: String,
    seq_id: String,
    start: usize,
    stop: Option<usize>,
    wraps: usize,
    length: usize,
    ratio: f64,
}

pub fn orfs(cmd: &Command) -> anyhow::Result<()> {
    match cmd {
        Command::Orfs {
            input,
            output,
            min_length,
            start_codons,
            stop_codons,
            include_stop,
            min_wraps,
            max_wraps,
            min_ratio,
            strand,
            no_stop_required,
            table,
            threads,
        } => {
            let reader = input_to_reader(input)?;
            let mut writer = output_to_writer(output)?;
            let mut table_writer = table_path_to_writer(table);

            // Step 1: Find all stop and start codons by frame
            let start_codons = start_codons.split(',').collect::<Vec<_>>();
            let stop_codons = stop_codons.split(',').collect::<Vec<_>>();

            parallel_fasta(
                reader,
                *threads,
                64,
                |record, orfs: &mut (Vec<circkit::orfs::Orf>, Vec<circkit::orfs::Orf>, Vec<u8>)| {
                    // runs in worker
                    let normalized = match needletail::sequence::normalize(record.seq(), false) {
                        Some(x) => x,
                        None => record.seq().to_vec(),
                    };

                    let (starts, stops) = circkit::orfs::start_stop_codon_indices_by_frame_naive(
                        std::str::from_utf8(&normalized).unwrap(),
                        &start_codons,
                        &stop_codons,
                    );

                    let mut all_orfs =
                        circkit::orfs::find_orfs_with_indices(normalized.len(), starts, stops);

                    // length filtering, stop codon requirement (with optional bypass), and wrap filtering
                    all_orfs.retain(|orf| {
                        (orf.length - 3 >= *min_length)
                            && (*no_stop_required || orf.stop.is_some())
                            && (*min_wraps <= orf.wraps)
                            && (orf.wraps <= *max_wraps)
                            && (orf.length as f64 / normalized.len() as f64 >= *min_ratio)
                    });

                    orfs.0 = circkit::orfs::longest_orfs(&mut all_orfs);

                    orfs.1 = if *strand == Strand::Both || *strand == Strand::Reverse {
                        orfs.2 = bio::alphabets::dna::revcomp(&normalized);
                        let (starts, stops) =
                            circkit::orfs::start_stop_codon_indices_by_frame_naive(
                                std::str::from_utf8(&orfs.2).unwrap(),
                                &start_codons,
                                &stop_codons,
                            );

                        let mut all_rc_orfs =
                            circkit::orfs::find_orfs_with_indices(normalized.len(), starts, stops);
                        all_rc_orfs.retain(|orf| {
                            (orf.length - 3 >= *min_length)
                                && (*no_stop_required || orf.stop.is_some())
                                && (*min_wraps <= orf.wraps)
                                && (orf.wraps <= *max_wraps)
                                && (orf.length as f64 / normalized.len() as f64 >= *min_ratio)
                        });
                        circkit::orfs::longest_orfs(&mut all_rc_orfs)
                    } else {
                        orfs.2.clear();
                        Vec::new()
                    };
                },
                |record, orfs| {
                    let head = std::str::from_utf8(record.head()).expect(
                        "Could not convert FASTA record header to UTF-8. Are you sure it's ASCII?",
                    );

                    for orf in &orfs.0 {
                        writer.write_all(b">").unwrap();
                        writer.write_all(record.head()).unwrap();
                        writer.write_all(b" ORF").unwrap();
                        writer.write_all(orf.start.to_string().as_bytes()).unwrap();
                        writer.write_all(b"\n").unwrap();
                        writer
                            .write_all(
                                orf.seq_with_opts(&record.full_seq(), *include_stop)
                                    .as_bytes(),
                            )
                            .unwrap();
                        writer.write_all(b"\n").unwrap();

                        // write the table file if it was requested
                        if let Some(ref mut table_writer) = table_writer {
                            table_writer
                                .serialize(Row {
                                    orf_id: format!("{} ORF{}", head, orf.start),
                                    seq_id: head.to_string(),
                                    start: orf.start,
                                    stop: orf.stop,
                                    wraps: orf.wraps,
                                    length: orf.length
                                        - match *include_stop {
                                            true => 0,
                                            false => 3,
                                        },
                                    ratio: orf.length as f64 / record.full_seq().len() as f64,
                                })
                                .expect("failed to write to table");
                        }
                    }
                    for orf in &orfs.1 {
                        writer.write_all(b">").unwrap();
                        writer.write_all(record.head()).unwrap();
                        writer.write_all(b" ORF").unwrap();
                        writer.write_all(orf.start.to_string().as_bytes()).unwrap();
                        writer.write_all(b" RC\n").unwrap();
                        writer
                            .write_all(orf.seq_with_opts(&orfs.2, *include_stop).as_bytes())
                            .unwrap();
                        writer.write_all(b"\n").unwrap();

                        // write the table file if it was requested
                        if let Some(ref mut table_writer) = table_writer {
                            table_writer
                                .serialize(Row {
                                    orf_id: format!("{} ORF{} RC", head, orf.start),
                                    seq_id: head.to_string(),
                                    start: &orfs.2.len() - 1 - orf.start, // reverse complement coordinates back to forward strand
                                    stop: match orf.stop {
                                        Some(x) => Some(&orfs.2.len() - 1 - x),
                                        None => None,
                                    },
                                    wraps: orf.wraps,
                                    length: orf.length
                                        - match *include_stop {
                                            true => 0,
                                            false => 3,
                                        },
                                    ratio: orf.length as f64 / record.full_seq().len() as f64,
                                })
                                .expect("failed to write to table");
                        }
                    }

                    // Some(value) will stop the reader, and the value will be returned.
                    // In the case of never stopping, we need to give the compiler a hint about the
                    // type parameter, thus the special 'turbofish' notation is needed,
                    // hoping on progress here: https://github.com/rust-lang/rust/issues/27336
                    None::<()>
                },
            )?;
            writer.flush()?;
            if let Some(mut table_writer) = table_writer {
                table_writer.flush()?;
            }
        }
        _ => panic!("input command is not for orfs"),
    }
    Ok(())
}
