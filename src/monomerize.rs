use anyhow::bail;
use needletail::parser::write_fasta;

use crate::{
    commands::Command,
    utils::{to_reader, to_writer},
};

pub fn monomerize(cmd: &Command) -> anyhow::Result<()> {
    match cmd {
        Command::Monomerize {
            input,
            output,
            normalize,
            seed_length,
            max_mismatch,
            min_identity,
            keep_all,
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

            let mut reader = to_reader(input)?;
            let mut writer = to_writer(output)?;

            while let Some(r) = reader.next() {
                let record = r?;

                let seq = &record.seq();

                let monomer = if *normalize {
                    circkit::normalize(
                        circkit::monomerize(
                            seq,
                            usize::try_from(*seed_length)?,
                            max_mismatch.unwrap_or(0),
                        ),
                        circkit::normalize::Alphabet::Dna,
                    )
                } else {
                    circkit::monomerize(
                        seq,
                        usize::try_from(*seed_length)?,
                        max_mismatch.unwrap_or(0),
                    )
                    .to_owned()
                };

                // if requested, keep the monomer even if it was unit length originally
                if !*keep_all && monomer.len() == seq.len() {
                    continue;
                }

                write_fasta(
                    record.id(),
                    &monomer,
                    &mut writer,
                    needletail::parser::LineEnding::Unix,
                )?;
            }

            // Clean up before exiting
            writer.flush()?;
            Ok(())
        }
        _ => panic!("input command is not for monomerize"),
    }
}
