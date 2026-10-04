use crate::{
    commands::Command,
    utils::{canonicalization_input, normalized_sequence, output_to_writer, process_fasta},
};
use seq_io::fasta::Record;
use std::io::Write;

pub fn canonicalize(cmd: &Command) -> anyhow::Result<()> {
    match cmd {
        Command::Canonicalize {
            input,
            output,
            threads,
            processing,
            rotation,
        } => {
            let (reader, settings) = canonicalization_input(input, output, processing, *threads)?;
            let reuse_buffers = settings.serial;
            let rotation = circkit::canonicalize::RotationOptions {
                duval_max_len: rotation.rotation_cutoff,
            };
            let mut writer = output_to_writer(output)?;

            process_fasta(
                reader,
                settings,
                |record, data: &mut (Vec<u8>, Vec<u8>)| {
                    // runs in worker

                    let normalized = normalized_sequence(record.seq());
                    if reuse_buffers {
                        rotation.canonicalize_into(&normalized, &mut data.0, &mut data.1);
                    } else {
                        data.0 = rotation.canonicalize(&normalized);
                        data.1 = Vec::new();
                    }
                },
                |record, data| {
                    // runs in main thread
                    writer.write_all(b">")?;
                    writer.write_all(record.head())?;
                    writer.write_all(b"\n")?;
                    writer.write_all(&data.0)?;
                    writer.write_all(b"\n")?;

                    Ok(())
                },
            )?;
            writer.finish()?;
        }
        _ => panic!("input command is not for canonicalize"),
    }
    Ok(())
}
