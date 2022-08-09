use crate::commands::Command;
use crate::utils::{to_reader, to_writer};
use needletail::{parser::write_fasta, Sequence};

/// Rotationally canonicalize sequences.
pub fn normalize(cmd: &Command) -> anyhow::Result<()> {
    match cmd {
        Command::Normalize { input, output } => {
            let mut reader = to_reader(input)?;
            let mut writer = to_writer(output)?;

            while let Some(r) = reader.next() {
                let record = r?;

                write_fasta(
                    record.id(),
                    &circkit::normalize(
                        &record.seq().normalize(false),
                        circkit::normalize::Alphabet::Dna,
                    ),
                    &mut writer,
                    needletail::parser::LineEnding::Unix,
                )?;
            }

            // Clean up before exiting
            writer.flush()?;
        }
        _ => panic!("input command is not for normalize"),
    };
    Ok(())
}
