use anyhow::bail;
use seq_io::{fasta::Record, parallel::parallel_fasta};

use crate::{
    commands::Command,
    utils::{input_to_reader, output_to_writer},
};

pub fn monomerize2(cmd: &Command) -> anyhow::Result<()> {
    match cmd {
        Command::Monomerize2 {
            input,
            output,
            normalize,
            seed_length,
            max_mismatch,
            min_identity,
            keep_all,
            threads,
            ..
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
                u32::try_from(*threads)?,
                64,
                |record, seq| {
                    let original_seq = record.full_seq();
                    let original_len = original_seq.len();
                    let mut builder = circkit::monomerize2::Monomerizer::builder();

                    // set the seed length
                    builder.seed_len((*seed_length).try_into().expect("Seed length is too large"));

                    if let Some(max_mismatch) = *max_mismatch {
                        builder.overlap_dist(max_mismatch);
                    }

                    let m = builder.build().unwrap();
                    *seq = m.monomerize(&original_seq).to_owned();

                    // if requested, keep the monomer even if it was unit length originally
                    if !*keep_all && seq.len() == original_len {
                        *seq = vec![];
                        return;
                    }

                    // optionally normalize the sequence to save on IO
                    if *normalize {
                        *seq = circkit::normalize(seq, circkit::normalize::Alphabet::Dna);
                    }
                },
                |record, seq| {
                    if !seq.is_empty() {
                        writer.write_all(b">").unwrap();
                        writer.write_all(record.id().unwrap().as_bytes()).unwrap();
                        writer.write_all(b"\n").unwrap();
                        writer.write_all(seq).unwrap();
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
