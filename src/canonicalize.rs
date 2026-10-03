use crate::{
    commands::Command,
    utils::{input_to_reader, normalized_sequence, output_to_writer, process_fasta},
};
use seq_io::fasta::Record;

pub fn canonicalize(cmd: &Command) -> anyhow::Result<()> {
    match cmd {
        Command::Canonicalize {
            input,
            output,
            threads,
            processing,
            rotation,
        } => {
            let settings = processing.resolve(*threads, 64)?;
            let reuse_buffers = settings.serial;
            let rotation = circkit::canonicalize::RotationOptions {
                duval_max_len: rotation.rotation_cutoff,
            };
            let reader = input_to_reader(input)?;
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
                    writer.write_all(b">").unwrap();
                    writer.write_all(record.head()).unwrap();
                    writer.write_all(b"\n").unwrap();
                    writer.write_all(&data.0).unwrap();
                    writer.write_all(b"\n").unwrap();

                    // Some(value) will stop the reader, and the value will be returned.
                    // In the case of never stopping, we need to give the compiler a hint about the
                    // type parameter, thus the special 'turbofish' notation is needed,
                    // hoping on progress here: https://github.com/rust-lang/rust/issues/27336
                    None::<()>
                },
            )?;
            writer.flush()?;
        }
        _ => panic!("input command is not for canonicalize"),
    }
    Ok(())
}
