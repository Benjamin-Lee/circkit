use crate::{
    commands::Command,
    utils::{
        input_to_reader, normalized_sequence, output_to_writer, process_fasta, table_path_to_writer,
    },
};
use anyhow::Context;
use seq_io::fasta::Record;
use std::io::Write;

#[derive(clap::ValueEnum, Clone, Debug, PartialEq)]
pub enum Strand {
    Forward,
    Reverse,
    Both,
}

#[derive(serde::Serialize, Debug)]
struct Row<'a> {
    orf_id: String,
    seq_id: &'a str,
    start: usize,
    stop: Option<usize>,
    length: usize,
    wraps: usize,
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
            table_format,
            threads,
            processing,
        } => {
            let settings =
                processing.resolve(*threads, (*threads as usize).saturating_mul(2).max(2))?;
            let reader = input_to_reader(input)?;
            let mut writer = output_to_writer(output)?;
            let mut table_writer = table_path_to_writer(table, *table_format, output)?;

            // Step 1: Find all stop and start codons by frame
            let start_codons = start_codons.split(',').collect::<Vec<_>>();
            let stop_codons = stop_codons.split(',').collect::<Vec<_>>();
            let matcher = circkit::orfs::CodonMatcher::new(&start_codons, &stop_codons);

            process_fasta(
                reader,
                settings,
                |record, orfs: &mut (Vec<circkit::orfs::Orf>, Vec<circkit::orfs::Orf>, Vec<u8>)| {
                    // runs in worker
                    let normalized = normalized_sequence(record.seq());
                    orfs.0 = if *strand != Strand::Reverse {
                        let (starts, stops) = matcher.indices(&normalized);

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

                        circkit::orfs::longest_orfs(&mut all_orfs)
                    } else {
                        Vec::new()
                    };

                    orfs.1 = if *strand == Strand::Both || *strand == Strand::Reverse {
                        orfs.2 = bio::alphabets::dna::revcomp(normalized.as_ref());
                        let (starts, stops) = matcher.indices(&orfs.2);

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
                    let head = if table_writer.is_some() {
                        std::str::from_utf8(record.head())
                            .context("FASTA headers must be UTF-8 when writing metadata")?
                    } else {
                        ""
                    };
                    let sequence_length = crate::io::sequence_len(&record);

                    // Join wrapped FASTA lines once per record, rather than once per ORF.
                    let forward_sequence = if orfs.0.is_empty() {
                        None
                    } else {
                        Some(record.full_seq())
                    };
                    for orf in &orfs.0 {
                        writer.write_all(b">")?;
                        writer.write_all(record.head())?;
                        writeln!(writer, "_ORF{}", orf.start)?;
                        orf.write_seq_with_opts(
                            forward_sequence.as_ref().unwrap(),
                            *include_stop,
                            &mut writer,
                        )?;
                        writer.write_all(b"\n")?;

                        // write the table file if it was requested
                        if let Some(ref mut table_writer) = table_writer {
                            table_writer.serialize(Row {
                                orf_id: format!("{}_ORF{}", head, orf.start),
                                seq_id: head,
                                start: orf.start,
                                stop: orf.stop,
                                length: orf.length
                                    - match *include_stop {
                                        true => 0,
                                        false => 3,
                                    },
                                wraps: orf.wraps,
                                ratio: orf.length as f64 / sequence_length as f64,
                            })?;
                        }
                    }
                    for orf in &orfs.1 {
                        writer.write_all(b">")?;
                        writer.write_all(record.head())?;
                        writeln!(writer, "_RC_ORF{}", orf.start)?;
                        orf.write_seq_with_opts(&orfs.2, *include_stop, &mut writer)?;
                        writer.write_all(b"\n")?;

                        // write the table file if it was requested
                        if let Some(ref mut table_writer) = table_writer {
                            table_writer.serialize(Row {
                                orf_id: format!("{}_RC_ORF{}", head, orf.start),
                                seq_id: head,
                                start: &orfs.2.len() - 1 - orf.start, // reverse complement coordinates back to forward strand
                                stop: match orf.stop {
                                    Some(x) => Some(&orfs.2.len() - 1 - x),
                                    None => None,
                                },
                                length: orf.length
                                    - match *include_stop {
                                        true => 0,
                                        false => 3,
                                    },
                                wraps: orf.wraps,
                                ratio: orf.length as f64 / sequence_length as f64,
                            })?;
                        }
                    }

                    Ok(())
                },
            )?;
            writer.finish()?;
            if let Some(table_writer) = table_writer {
                table_writer.finish()?;
            }
        }
        _ => panic!("input command is not for orfs"),
    }
    Ok(())
}
