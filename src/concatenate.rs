use crate::{
    commands::Command,
    utils::{input_to_reader, output_to_writer},
};
use seq_io::fasta::Record;

/// Concatenate sequences to themselves.
///
/// This can be useful when using circular sequences with tools that don't directly support circular sequences.
pub fn concatenate(cmd: &Command) -> anyhow::Result<()> {
    match cmd {
        Command::Cat { input, output } => {
            let mut reader = input_to_reader(input)?;
            let mut writer = output_to_writer(output)?;

            while let Some(Ok(record)) = reader.next() {
                let full_seq = record.full_seq();
                writer.write_all(b">")?;
                writer.write_all(record.head())?;
                writer.write_all(b"\n")?;
                writer.write_all(&full_seq)?;
                writer.write_all(&full_seq)?;
                writer.write_all(b"\n")?;
            }

            writer.flush()?;

            Ok(())
        }
        _ => panic!("Wrong command"),
    }
}

pub fn deconcatenate(cmd: &Command) -> anyhow::Result<()> {
    match cmd {
        Command::Decat { input, output } => {
            let mut reader = input_to_reader(input)?;
            let mut writer = output_to_writer(output)?;

            while let Some(Ok(record)) = reader.next() {
                let full_seq = record.full_seq();
                writer.write_all(b">")?;
                writer.write_all(record.head())?;
                writer.write_all(b"\n")?;
                writer.write_all(&full_seq[..full_seq.len() / 2])?;
                writer.write_all(b"\n")?;
            }

            writer.flush()?;

            Ok(())
        }
        _ => panic!("Wrong command"),
    }
}
