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
        } => {
            let mut reader = to_reader(&input)?;
            let mut writer = to_writer(&output)?;

            while let Some(r) = reader.next() {
                let record = r?;

                write_fasta(
                    record.id(),
                    circkit::monomerize(&record.seq(), 10, 1),
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
