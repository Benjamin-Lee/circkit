use anyhow::bail;
use seq_io::fasta::Reader;
use std::{
    fs::File,
    io::{prelude::*, stdin, stdout, BufReader, BufWriter},
    path::PathBuf,
};

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
            let suffix = output.extension().unwrap_or_default().to_str().unwrap();

            let compression_format = match suffix {
                "gz" => niffler::send::compression::Format::Gzip,
                "bz2" => niffler::send::compression::Format::Bzip,
                "xz" => niffler::send::compression::Format::Lzma,
                "zst" => niffler::send::compression::Format::Zstd,
                _ => niffler::send::compression::Format::No,
            };

            let outfile = match File::create(output) {
                Ok(file) => file,
                Err(_) => {
                    bail!(
                        "Could not create output file {}. Are you sure it's not actually a directory?",
                        output.display()
                    );
                }
            };

            let fp_bufwriter = BufWriter::new(outfile);
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

pub fn table_path_to_writer(table: &Option<PathBuf>) -> Option<csv::Writer<File>> {
    table.as_ref().map(|path| {
        csv::WriterBuilder::new()
            .delimiter(match path.extension().and_then(|x| x.to_str()) {
                Some("tsv") => b'\t',
                _ => b',',
            })
            .from_path(path)
            .expect("Could not create output table.")
    })
}
