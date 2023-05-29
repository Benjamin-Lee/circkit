use anyhow::bail;
use seq_io::{fasta::Record, parallel::parallel_fasta};

use crate::{
    commands::Command,
    utils::{input_to_reader, output_to_writer},
};

pub fn monomerize(cmd: &Command) -> anyhow::Result<()> {
    match cmd {
        Command::Monomerize {
            input,
            output,
            seed_length,
            max_mismatch,
            min_identity,
            keep_all,
            threads,
            batch_size,
            min_overlap,
            min_overlap_percent,
        } => {
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

            parallel_fasta(
                reader,
                *threads,
                *batch_size,
                |record, idx| {
                    let mut builder = circkit::monomerize::Monomerizer::builder();

                    // set the seed length
                    builder.seed_len((*seed_length).try_into().expect("Seed length is too large"));

                    // set the maximum mismatch count
                    if let Some(max_mismatch) = *max_mismatch {
                        builder.overlap_dist(max_mismatch);
                    }

                    // set the minimum identity
                    if let Some(min_identity) = *min_identity {
                        builder.overlap_min_identity(min_identity);
                    }

                    let m = builder.build().unwrap();

                    // normalize the sequence
                    let normalized = match needletail::sequence::normalize(record.seq(), false) {
                        Some(x) => x,
                        None => record.seq().to_vec(),
                    };

                    // make sure the sequence is at least as long as the seed length
                    if normalized.len() < m.seed_len {
                        *idx = None;
                        return;
                    }

                    *idx = m.last_monomer_end_index(&normalized);
                },
                |record, idx| {
                    // get the full sequence
                    let full_seq = record.full_seq();

                    // region: check the monomer is long enough, either absolute or relative to the original sequence
                    // absolute length
                    if let Some(min_overlap) = *min_overlap {
                        if let Some(last_monomer_end_index) = *idx {
                            if last_monomer_end_index < min_overlap {
                                *idx = None; // reject the monomer
                            }
                        }
                    }

                    // relative length
                    if let Some(min_overlap_percent) = *min_overlap_percent {
                        if let Some(last_monomer_end_index) = *idx {
                            let monomer_length = last_monomer_end_index;
                            // dbg!(record.id().unwrap());
                            // dbg!(
                            //     (full_seq.len() - monomer_length) as f64 / (monomer_length as f64)
                            // );
                            // dbg!(min_overlap_percent);
                            // full_seq.len() - monomer_length is the length of the overlapping region
                            // a complete monomer would have an overlap ratio of 1.0
                            if (full_seq.len() - monomer_length) as f64 / (monomer_length as f64)
                                < min_overlap_percent
                            {
                                *idx = None; // reject the monomer
                            }
                        }
                    }
                    // endregion

                    // when keep_all is true, we write all sequences
                    // otherwise, we only write sequences that have been monomerized (i.e. the monomer index is Some)
                    if idx.is_some() || *keep_all {
                        let end_idx = idx.unwrap_or(full_seq.len());
                        writer.write_all(b">").unwrap();
                        writer.write_all(record.id().unwrap().as_bytes()).unwrap();
                        writer.write_all(b"\n").unwrap();
                        writer.write_all(&full_seq[..end_idx]).unwrap();
                        writer.write_all(b"\n").unwrap();
                    }
                    None::<()>
                },
            )?;
            writer.flush()?;

            Ok(())
        }
        _ => panic!("input command is not for monomerize"),
    }
}
