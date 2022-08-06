use crate::utils::{get_bufwriter, get_reader};
use bio::io::fasta;
use log::warn;
use std::io;
use std::{io::Write, path::PathBuf};

/// Concatenate sequences to themselves.
///
/// This can be useful when using circular sequences with tools that don't directly support circular sequences.
pub fn concatenate(input: &Option<PathBuf>, output: &Option<PathBuf>) -> io::Result<()> {
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
