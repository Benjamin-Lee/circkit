use crate::{
    commands::Command,
    utils::{input_to_reader, output_to_writer},
};
use seq_io::{fasta::Record, parallel::parallel_fasta};

pub fn orfs(cmd: &Command) -> anyhow::Result<()> {
    match cmd {
        Command::Orfs {
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
                |record, orfs| {
                    // runs in worker
                    let normalized = match needletail::sequence::normalize(record.seq(), false) {
                        Some(x) => x,
                        None => record.seq().to_vec(),
                    };

                    let mut all_orfs =
                        circkit::orfs::find_orfs(std::str::from_utf8(&normalized).unwrap());
                    *orfs = circkit::orfs::longest_orfs(&mut all_orfs);
                },
                |record, orfs| {
                    for orf in orfs {
                        writer.write_all(b">").unwrap();
                        writer.write_all(record.head()).unwrap();
                        writer.write_all(b" ORF").unwrap();
                        writer.write_all(orf.start.to_string().as_bytes()).unwrap();
                        writer.write_all(b"\n").unwrap();
                        writer
                            .write_all(orf.seq(&record.full_seq()).as_bytes())
                            .unwrap();
                        writer.write_all(b"\n").unwrap();
                    }

                    // Some(value) will stop the reader, and the value will be returned.
                    // In the case of never stopping, we need to give the compiler a hint about the
                    // type parameter, thus the special 'turbofish' notation is needed,
                    // hoping on progress here: https://github.com/rust-lang/rust/issues/27336
                    None::<()>
                },
            )?;
            writer.flush()?;
        }
        _ => panic!("input command is not for orfs"),
    }
    Ok(())
}
