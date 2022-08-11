use crate::commands::Command;
use seq_io::{
    fasta::{Reader, Record},
    parallel::parallel_fasta,
};
use std::{
    fs::File,
    io::{prelude::*, stdin, stdout, BufReader, BufWriter},
    path::PathBuf,
};

fn input_to_reader(input: &Option<PathBuf>) -> anyhow::Result<Reader<Box<dyn Read + Send>>> {
    match input {
        Some(input) => {
            let fp_bufreader = BufReader::new(File::open(input)?);
            let niffed = niffler::send::get_reader(Box::new(fp_bufreader))?.0;
            let reader = Reader::new(niffed);
            Ok(reader)
        }
        None => {
            let stdin_bufreader = BufReader::new(stdin());
            let niffed = niffler::send::get_reader(Box::new(stdin_bufreader))?.0;
            let reader = Reader::new(niffed);
            Ok(reader)
        }
    }
}

fn output_to_writer(output: &Option<PathBuf>) -> anyhow::Result<Box<dyn Write>> {
    match output {
        Some(output) => {
            let fp_bufwriter = BufWriter::new(File::create(output)?);
            Ok(Box::new(fp_bufwriter))
        }
        None => {
            let stdout_bufwriter = BufWriter::new(stdout());
            Ok(Box::new(stdout_bufwriter))
        }
    }
}

pub fn normalize2(cmd: &Command) -> anyhow::Result<()> {
    match cmd {
        Command::Normalize2 { input, output } => {
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
            )
            .unwrap();
        }
        _ => {
            anyhow::bail!("foo is not implemented for {:?}", cmd);
        }
    }
    Ok(())
}
