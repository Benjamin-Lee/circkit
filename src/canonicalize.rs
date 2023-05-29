use crate::{
    commands::Command,
    utils::{input_to_reader, output_to_writer},
};
use seq_io::{fasta::Record, parallel::parallel_fasta};

pub fn canonicalize(cmd: &Command) -> anyhow::Result<()> {
    match cmd {
        Command::Canonicalize {
            input,
            output,
            threads,
        } => {
            let reader = input_to_reader(input)?;
            let mut writer = output_to_writer(output)?;

            parallel_fasta(
                reader,
                *threads,
                64,
                |record, seq| {
                    // runs in worker

                    let normalized = match needletail::sequence::normalize(record.seq(), false) {
                        Some(x) => x,
                        None => record.seq().to_vec(),
                    };

                    *seq = circkit::canonicalize(&normalized);
                },
                |record, seq| {
                    // runs in main thread
                    writer.write_all(b">").unwrap();
                    writer.write_all(record.head()).unwrap();
                    writer.write_all(b"\n").unwrap();
                    writer.write_all(seq).unwrap();
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
