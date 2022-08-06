use anyhow::{bail, Result};
use bio::io::fasta;
use std::{
    fs::File,
    io::{self, stdout, BufRead, BufWriter, Read, Write},
    path::PathBuf,
};

pub fn get_bufwriter(output: &Option<PathBuf>) -> BufWriter<Box<dyn Write>> {
    const BUFFER_CAPACITY: usize = 64 * 1024;

    // decide where to write to
    let stdout = stdout();
    let file;
    match output {
        Some(path) => {
            file = File::create(path).expect("Unable to create output file");
            BufWriter::with_capacity(BUFFER_CAPACITY, Box::new(file))
        }
        None => BufWriter::with_capacity(BUFFER_CAPACITY, Box::new(stdout.lock())),
    }
}

pub fn get_reader(input: &Option<PathBuf>) -> Box<dyn BufRead> {
    let stdin = io::stdin();
    match input {
        Some(path) => Box::new(io::BufReader::new(
            File::open(path).expect("Unable to create output file"),
        )),
        None => Box::new(stdin.lock()),
    }
}

// #[derive(Debug)]
// pub enum Input {
//     Stdin,
//     File(PathBuf),
// }

// impl Input {
//     ///hi
//     pub fn from(path: Option<PathBuf>) -> Result<Self> {
//         let input = match path {
//             Some(path) => Input::File(path),
//             None => Input::Stdin,
//         };

//         match input {
//             Input::Stdin => {
//                 if atty::is(atty::Stream::Stdin) {
//                     bail!("No stdin detected. Did you mean to include a file argument?");
//                 }
//             }
//             Input::File(ref path) => {
//                 let file = File::open(path)?;
//                 if !file.metadata()?.is_file() {
//                     bail!(
//                         "Input {} is a directory. Did you mean to include a file in the directory?",
//                         path.display()
//                     );
//                 }
//             }
//         }

//         Ok(input)
//     }
// }
