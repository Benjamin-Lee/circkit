use anyhow::{bail, Result};
use needletail::{parse_fastx_file, parse_fastx_stdin, FastxReader};
use seq_io::fasta::Reader;
use std::{
    fs::File,
    io::{prelude::*, stdin, stdout, BufReader, BufWriter},
    path::PathBuf,
};

pub fn to_reader(input: &Option<PathBuf>) -> Result<Box<dyn FastxReader>> {
    match input {
        Some(filename) => Ok(parse_fastx_file(&filename).expect("valid path/file")),
        None => {
            if atty::is(atty::Stream::Stdin) {
                bail!("No stdin detected. Did you mean to include a file argument?");
            }
            Ok(parse_fastx_stdin().expect("valid stdin"))
        }
    }
}
pub fn to_writer(output: &Option<PathBuf>) -> Result<Box<dyn Write>> {
    match output {
        Some(filename) => {
            let file = File::create(filename)?;
            if !file.metadata()?.is_file() {
                bail!(
                    "Input {} is a directory. Did you mean to include a file in the directory?",
                    filename.display()
                );
            }
            Ok(Box::new(BufWriter::new(file)))
        }
        None => Ok(Box::new(BufWriter::new(stdout()))),
    }
}

pub fn input_to_reader(input: &Option<PathBuf>) -> anyhow::Result<Reader<Box<dyn Read + Send>>> {
    match input {
        Some(input) => {
            let fp_bufreader = BufReader::new(File::open(input)?);
            let niffed = niffler::send::get_reader(Box::new(fp_bufreader))?.0;
            let reader = Reader::new(niffed);
            Ok(reader)
        }
        None => {
            if atty::is(atty::Stream::Stdin) {
                bail!("No stdin detected. Did you mean to include a file argument?");
            }
            let stdin_bufreader = BufReader::new(stdin());
            let niffed = niffler::send::get_reader(Box::new(stdin_bufreader))?.0;
            let reader = Reader::new(niffed);
            Ok(reader)
        }
    }
}

pub fn output_to_writer(output: &Option<PathBuf>) -> anyhow::Result<Box<dyn Write>> {
    match output {
        Some(output) => {
            // match the suffix of outout to see if it should be compressed
            let suffix = output.extension().unwrap().to_str().unwrap();

            let compression_format = match suffix {
                "gz" => niffler::send::compression::Format::Gzip,
                "bz2" => niffler::send::compression::Format::Bzip,
                "xz" => niffler::send::compression::Format::Lzma,
                "zst" => niffler::send::compression::Format::Zstd,
                _ => niffler::send::compression::Format::No,
            };

            let fp_bufwriter = BufWriter::new(File::create(output)?);
            let niffed = niffler::send::get_writer(
                Box::new(fp_bufwriter),
                compression_format,
                match compression_format {
                    niffler::send::compression::Format::Gzip => niffler::compression::Level::Six,
                    niffler::send::compression::Format::Bzip => niffler::compression::Level::Nine,
                    niffler::send::compression::Format::Lzma => niffler::compression::Level::Six,
                    niffler::send::compression::Format::Zstd => niffler::compression::Level::One,
                    niffler::send::compression::Format::No => niffler::compression::Level::One,
                },
            )?;
            Ok(niffed)
        }
        None => {
            let stdout_bufwriter = BufWriter::new(stdout());
            Ok(Box::new(stdout_bufwriter))
        }
    }
}
