use crate::commands::Command;
use needletail::{parse_fastx_file, parse_fastx_stdin, parser::write_fasta, Sequence};
use std::fs::File;
use std::io::{prelude::*, stdout, BufWriter};

/// Testing
pub fn normalize(cmd: &Command) -> anyhow::Result<()> {
    let mut reader;
    let mut writer: Box<dyn Write>;

    match cmd {
        Command::Normalize { input, output } => {
            reader = match input {
                Some(filename) => parse_fastx_file(&filename).expect("valid path/file"),
                None => parse_fastx_stdin().expect("valid stdin"),
            };
            writer = match output {
                Some(filename) => {
                    let file = File::create(filename).expect("valid path/file");
                    Box::new(BufWriter::new(file))
                }
                None => Box::new(BufWriter::new(stdout())),
            };
        }
        _ => panic!("input command is not for normalize"),
    };

    while let Some(r) = reader.next() {
        let record = r?;

        write_fasta(
            record.id(),
            &circkit::normalize::normalize(
                &record.seq().normalize(false),
                circkit::normalize::Alphabet::Dna,
            ),
            &mut writer,
            needletail::parser::LineEnding::Unix,
        )?;
    }

    // Clean up before exiting
    writer.flush()?;
    Ok(())
}
