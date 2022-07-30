use bio::io::fasta;
use std::{
    fs::File,
    io::{self, stdout, BufWriter, Read, Write},
    path::PathBuf,
};

fn get_bufwriter(output: &Option<PathBuf>) -> BufWriter<Box<dyn Write>> {
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

fn get_reader(input: &Option<PathBuf>) -> Box<dyn Read> {
    let stdin = io::stdin();
    match input {
        Some(path) => Box::new(File::open(path).expect("Unable to create output file")),
        None => Box::new(stdin.lock()),
    }
}

/// Macro to bail early from functions where the input is empty and stdin is a tty (not really readable).
macro_rules! check_tty {
    ($input: ident) => {
        if $input.is_none() && atty::is(atty::Stream::Stdin) {
            warn!("No input file specified and stdin is not detected");
            return Ok(());
        }
    };
}

/// Concatenate sequences to themselves.
///
/// This can be useful when using circular sequences with tools that don't directly support circular sequences.
pub fn concatenate(input: &Option<PathBuf>, output: &Option<PathBuf>) -> Result<(), io::Error> {
    check_tty!(input);

    let mut records = fasta::Reader::new(get_reader(input)).records();

    let mut handle = get_bufwriter(output);

    while let Some(Ok(record)) = records.next() {
        handle.write_all(b">")?;
        handle.write_all(record.id().as_bytes())?;

        // if there is a description, write it too
        if let Some(desc) = record.desc() {
            handle.write_all(b" ")?;
            handle.write_all(desc.as_bytes())?;
        }

        handle.write_all(b"\n")?;
        handle.write_all(record.seq())?;
        handle.write_all(record.seq())?;
        handle.write_all(b"\n")?;
    }

    // make sure to flush the buffer
    handle.flush()?;

    Ok(())
}

pub fn deconcatenate(input: &Option<PathBuf>, output: &Option<PathBuf>) -> io::Result<()> {
    check_tty!(input);

    let mut records = fasta::Reader::new(get_reader(input)).records();

    let mut handle = get_bufwriter(output);

    while let Some(Ok(record)) = records.next() {
        // Check if sequence is even length and warn if not
        if record.seq().len() % 2 != 0 {
            warn!(
                "Sequence \"{}\" is not even length ({}). Are you sure the monomer is concatenated to itself?",
                record.id(),
                record.seq().len(),
            );
        }

        handle.write_all(b">")?;
        handle.write_all(record.id().as_bytes())?;

        // if there is a description, write it too
        if let Some(desc) = record.desc() {
            handle.write_all(b" ")?;
            handle.write_all(desc.as_bytes())?;
        }

        handle.write_all(b"\n")?;
        handle.write_all(&record.seq()[0..record.seq().len() / 2])?;
        handle.write_all(b"\n")?;
    }

    // make sure to flush the buffer
    handle.flush()?;

    Ok(())
}
