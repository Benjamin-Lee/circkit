use anyhow::bail;
use seq_io::fasta::Record;
use std::io::Write;

use crate::{
    commands::Command,
    utils::{input_to_reader, output_to_writer},
};

pub fn rotate(cmd: &Command) -> anyhow::Result<()> {
    match cmd {
        Command::Rotate {
            input,
            output,
            bases,
            percent,
        } => {
            let mut reader = input_to_reader(input)?;
            let mut writer = output_to_writer(output)?;

            while let Some(record) = reader.next() {
                let record = record?;
                let length = crate::io::sequence_len(&record);

                let new_start_index = match percent {
                    Some(percent) => {
                        let shift = (length as f64 * percent).floor();
                        if !shift.is_finite()
                            || shift < i64::MIN as f64
                            || shift >= -(i64::MIN as f64)
                        {
                            return Err(crate::diagnostics::argument_error(
                                "rotation percentage is too large for this record",
                            ));
                        }
                        shift as i64
                    }
                    None => match bases {
                        Some(bases) => *bases,
                        None => bail!("provide either --bases or --percent"),
                    },
                };

                writer.write_all(b">")?;
                writer.write_all(record.head())?;
                writer.write_all(b"\n")?;

                let rotation_index = match (length, new_start_index >= 0) {
                    (0, _) => 0,
                    (_, true) => length - (new_start_index as u64 % length as u64) as usize,
                    (_, false) => (new_start_index.unsigned_abs() % length as u64) as usize,
                };

                crate::io::write_sequence_range(&mut writer, &record, rotation_index..length)?;
                crate::io::write_sequence_range(&mut writer, &record, 0..rotation_index)?;
                writer.write_all(b"\n")?;
            }
            writer.finish()?;
            Ok(())
        }
        _ => panic!("This should never happen"),
    }
}
