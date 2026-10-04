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

            while let Some(record) = reader.next() {
                let record = record?;
                writer.write_all(b">")?;
                writer.write_all(record.head())?;
                writer.write_all(b"\n")?;
                for _ in 0..2 {
                    for line in record.seq_lines() {
                        writer.write_all(line)?;
                    }
                }
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

            while let Some(record) = reader.next() {
                let record = record?;
                writer.write_all(b">")?;
                writer.write_all(record.head())?;
                writer.write_all(b"\n")?;
                crate::io::write_sequence_range(
                    &mut writer,
                    &record,
                    0..crate::io::sequence_len(&record) / 2,
                )?;
                writer.write_all(b"\n")?;
            }

            writer.finish()?;

            Ok(())
        }
        _ => panic!("Wrong command"),
    }
}
