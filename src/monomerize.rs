use anyhow::bail;
use anyhow::Context;
use seq_io::fasta::Record;
use std::io::Write;

use crate::{
    commands::Command,
    utils::{
        input_to_reader, normalized_sequence, output_to_writer, process_fasta, table_path_to_writer,
    },
};

#[derive(serde::Serialize)]
struct Row<'a> {
    id: &'a str,
    original_length: usize,
    monomer_length: usize,
}

pub fn monomerize(cmd: &Command) -> anyhow::Result<()> {
    match cmd {
        Command::Monomerize {
            input,
            output,
            sensitive,
            seed_length,
            max_mismatch,
            min_identity,
            min_overlap,
            min_overlap_percent,
            min_length,
            max_length,
            keep_all,
            table,
            table_format,
            threads,
            mismatch_chunk_size,
            processing,
        } => {
            let settings = processing.resolve(*threads, 64)?;
            // region: some basic sanity checks
            if max_mismatch.is_some() && min_identity.is_some() {
                bail!("cannot specify both max_mismatch and min_identity");
            }

            // make sure the minimum identity is in range
            if let Some(min_identity) = *min_identity {
                if !(0.0..=1.0).contains(&min_identity) {
                    bail!("min_identity must be between 0.0 and 1.0");
                }
            }
            // endregion

            let reader = input_to_reader(input)?;
            let mut writer = output_to_writer(output)?;
            let mut table_writer = table_path_to_writer(table, *table_format)?;
            let mut scratch = Vec::new();

            let mut builder = circkit::monomerize::Monomerizer::builder();

            // set the seed length
            builder.seed_len(usize::try_from(*seed_length).context("seed length is too large")?);
            builder.mismatch_chunk_size(*mismatch_chunk_size);

            // set the maximum mismatch count
            if let Some(max_mismatch) = *max_mismatch {
                builder.overlap_dist(max_mismatch);
            }

            // set the minimum identity
            if let Some(min_identity) = *min_identity {
                builder.overlap_min_identity(min_identity);
            }

            let monomerizer = builder.build().context("configure monomerizer")?;

            process_fasta(
                reader,
                settings,
                |record, idx| {
                    // normalize the sequence
                    let normalized = normalized_sequence(record.seq());

                    // make sure the sequence is at least as long as the seed length and the minimum length
                    if normalized.len() < monomerizer.seed_len || normalized.len() < *min_length {
                        *idx = None;
                        return;
                    }

                    *idx = match sensitive {
                        true => monomerizer.last_monomer_end_index_sensitive(&normalized),
                        false => monomerizer.last_monomer_end_index(&normalized),
                    }
                },
                |record, idx| {
                    let sequence = crate::io::joined_sequence(&record, &mut scratch);
                    let full_length = sequence.len();

                    // region: check the monomer is long enough, either absolute or relative to the original sequence

                    // absolute monomer length
                    if let Some(monomer_length) = *idx {
                        if monomer_length < *min_length
                            || monomer_length > max_length.unwrap_or(usize::MAX)
                        {
                            *idx = None; // reject the monomer
                        }
                    }

                    // absolute overlap length
                    if let Some(min_overlap) = *min_overlap {
                        if let Some(monomer_length) = *idx {
                            if full_length - monomer_length < min_overlap {
                                *idx = None; // reject the monomer
                            }
                        }
                    }

                    // relative overlap length
                    if let Some(min_overlap_percent) = *min_overlap_percent {
                        if let Some(monomer_length) = *idx {
                            // full_length - monomer_length is the length of the overlapping region
                            // a complete monomer would have an overlap ratio of 1.0
                            if (full_length - monomer_length) as f64 / (monomer_length as f64)
                                < min_overlap_percent
                            {
                                *idx = None; // reject the monomer
                            }
                        }
                    }
                    // endregion

                    // when keep_all is true, we write all sequences
                    // otherwise, we only write sequences that have been monomerized (i.e. the monomer index is Some)
                    if (idx.is_some()) || *keep_all {
                        let end_idx = idx.unwrap_or(full_length);
                        writer.write_all(b">")?;
                        writer.write_all(record.head())?;
                        writer.write_all(b"\n")?;
                        writer.write_all(&sequence[..end_idx])?;
                        writer.write_all(b"\n")?;

                        // write the table file if it was requested
                        if let Some(ref mut table_writer) = table_writer {
                            table_writer.serialize(Row {
                                id: std::str::from_utf8(record.head())
                                    .context("FASTA headers must be UTF-8 when writing metadata")?,
                                original_length: full_length,
                                monomer_length: end_idx,
                            })?
                        }
                    }
                    Ok(())
                },
            )?;
            writer.finish()?;
            if let Some(table_writer) = table_writer {
                table_writer.finish()?;
            }
            Ok(())
        }
        _ => panic!("input command is not for monomerize"),
    }
}
