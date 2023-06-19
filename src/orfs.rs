use crate::{
    commands::Command,
    utils::{input_to_reader, output_to_writer},
};
use seq_io::{fasta::Record, parallel::parallel_fasta};

#[derive(clap::ArgEnum, Clone, Debug)]
pub enum Strand {
    Forward,
    Reverse,
    Both,
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
            strand,
            no_stop_required,
            threads,
        } => {
            let reader = input_to_reader(input)?;
            let mut writer = output_to_writer(output)?;

            // Step 1: Find all stop and start codons by frame
            let start_codons = start_codons.split(',').collect::<Vec<_>>();
            let stop_codons = stop_codons.split(',').collect::<Vec<_>>();

            parallel_fasta(
                reader,
                *threads,
                64,
                |record, orfs: &mut (Vec<circkit::orfs::Orf>, Vec<circkit::orfs::Orf>)| {
                    // runs in worker
                    let normalized = match needletail::sequence::normalize(record.seq(), false) {
                        Some(x) => x,
                        None => record.seq().to_vec(),
                    };

                    let (starts, stops) = circkit::orfs::start_stop_codon_indices_by_frame_naive(
                        &std::str::from_utf8(&normalized).unwrap(),
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
                    });

                    orfs.0 = circkit::orfs::longest_orfs(&mut all_orfs);
                },
                |record, orfs| {
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
                    }

                    // Some(value) will stop the reader, and the value will be returned.
                    // In the case of never stopping, we need to give the compiler a hint about the
                    // type parameter, thus the special 'turbofish' notation is needed,
                    // hoping on progress here: https://github.com/rust-lang/rust/issues/27336
                    None::<()>
                },
            )?;
            writer.flush()?;
        }
        _ => panic!("input command is not for orfs"),
    }
    Ok(())
}
