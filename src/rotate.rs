use anyhow::bail;
use seq_io::fasta::Record;

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

            // ensure bases and percent aren't 0
            if bases == &Some(0) || percent == &Some(0.0) {
                bail!("Rotation by 0 is not allowed");
            }

            while let Some(Ok(record)) = reader.next() {
                let full_seq = record.full_seq();

                let new_start_index = match percent {
                    Some(percent) => f64::floor(full_seq.len() as f64 * percent) as i64,
                    None => bases.expect("Must provide either --bases or --percent"),
                };

                writer.write_all(b">").unwrap();
                writer.write_all(record.head()).unwrap();
                writer.write_all(b"\n").unwrap();

                let rotation_index = match new_start_index >= 0 {
                    true => full_seq.len() - (new_start_index as usize % full_seq.len()),
                    false => new_start_index.abs() as usize % full_seq.len(),
                };

                writer.write_all(&full_seq[rotation_index..]).unwrap();
                writer.write_all(&full_seq[..rotation_index]).unwrap();
                writer.write_all(b"\n").unwrap();
            }

            Ok(())
        }
        _ => panic!("This should never happen"),
    }
}
