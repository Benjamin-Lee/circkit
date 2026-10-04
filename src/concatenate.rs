use crate::{
    commands::Command,
    utils::{input_to_reader, output_to_writer},
};
use seq_io::fasta::Record;
use std::io::Write;

/// Concatenate sequences to themselves.
///
/// This can be useful when using circular sequences with tools that don't directly support circular sequences.
pub fn concatenate(cmd: &Command) -> anyhow::Result<()> {
    match cmd {
        Command::Cat { input, output } => {
            let mut reader = input_to_reader(input)?;
            let mut writer = output_to_writer(output)?;
            let mut scratch = Vec::new();

            while let Some(record) = reader.next() {
                let record = record?;
                writer.write_all(b">")?;
                writer.write_all(record.head())?;
                writer.write_all(b"\n")?;
                let sequence = crate::io::joined_sequence(&record, &mut scratch);
                writer.write_all(sequence)?;
                writer.write_all(sequence)?;
                writer.write_all(b"\n")?;
            }

            writer.finish()?;

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
            let mut scratch = Vec::new();

            while let Some(record) = reader.next() {
                let record = record?;
                writer.write_all(b">")?;
                writer.write_all(record.head())?;
                writer.write_all(b"\n")?;
                let sequence = crate::io::joined_sequence(&record, &mut scratch);
                writer.write_all(&sequence[..sequence.len() / 2])?;
                writer.write_all(b"\n")?;
            }

            writer.finish()?;

            Ok(())
        }
        _ => panic!("Wrong command"),
    }
}
