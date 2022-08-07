use anyhow::{bail, Result};
use needletail::{parse_fastx_file, parse_fastx_stdin, FastxReader};
use std::{
    fs::File,
    io::{stdout, BufWriter, Write},
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
