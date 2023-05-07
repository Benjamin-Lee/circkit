use crate::{
    commands::Command,
    utils::{input_to_reader, output_to_writer},
};
use seq_io::{fasta::Record, parallel::parallel_fasta};

pub fn normalize(cmd: &Command) -> anyhow::Result<()> {
    match cmd {
        Command::Normalize { input, output } => {
            let reader = input_to_reader(input)?;
            let mut writer = output_to_writer(output)?;

            parallel_fasta(
                reader,
                8,
                64,
                |record, seq| {
                    // runs in worker
                    *seq =
                        circkit::normalize(&record.full_seq(), circkit::normalize::Alphabet::Dna);
                },
                |record, seq| {
                    // runs in main thread
                    writer.write_all(b">").unwrap();
                    writer.write_all(record.id().unwrap().as_bytes()).unwrap();
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
        _ => panic!("input command is not for normalize"),
    }
    Ok(())
}
