//! Streaming input, buffered output, and explicit compression finalization.
use crate::commands::TableFormat;
use anyhow::{Context, Result};
use niffler::send::compression::Format;
use seq_io::fasta::{Reader, Record, RefRecord};
use serde::Serialize;
use std::{
    fs::File,
    io::{self, BufReader, BufWriter, IsTerminal, Read, StdoutLock, Write},
    path::{Path, PathBuf},
};

pub type FastaReader = Reader<Box<dyn Read + Send>>;

pub fn is_stdio(path: &Path) -> bool {
    path.as_os_str() == "-"
}

fn destination_path(path: &Path) -> io::Result<PathBuf> {
    let absolute = std::path::absolute(path)?;
    if let Some(parent) = absolute.parent() {
        if let (Ok(parent), Some(name)) = (std::fs::canonicalize(parent), absolute.file_name()) {
            return Ok(parent.join(name));
        }
    }
    Ok(absolute)
}

fn same_path(first: &Path, second: &Path) -> io::Result<bool> {
    if destination_path(first)? == destination_path(second)? {
        return Ok(true);
    }
    match same_file::is_same_file(first, second) {
        Ok(same) => Ok(same),
        Err(error) if error.kind() == io::ErrorKind::NotFound => Ok(false),
        Err(error) => Err(error),
    }
}

pub fn validate_paths(
    input: &Option<PathBuf>,
    output: &Option<PathBuf>,
    table: Option<&PathBuf>,
) -> Result<()> {
    let input = input.as_deref().filter(|path| !is_stdio(path));
    let output = output.as_deref().filter(|path| !is_stdio(path));
    let table_file = table.map(PathBuf::as_path).filter(|path| !is_stdio(path));
    if table.is_some_and(|path| is_stdio(path)) && output.is_none() {
        return Err(crate::diagnostics::argument_error(
            "FASTA and metadata cannot share stdout; use -o FILE with --table -",
        ));
    }
    for (first, second, message) in [
        (
            input,
            output,
            "input and FASTA output must be different files",
        ),
        (
            input,
            table_file,
            "input and metadata output must be different files",
        ),
        (
            output,
            table_file,
            "FASTA and metadata outputs must be different files",
        ),
    ] {
        if let (Some(first), Some(second)) = (first, second) {
            if same_path(first, second).context("check input/output file paths")? {
                return Err(crate::diagnostics::argument_error(message));
            }
        }
    }
    Ok(())
}

pub fn input_to_reader(input: &Option<PathBuf>) -> Result<FastaReader> {
    Ok(input_to_reader_with_format(input)?.0)
}

pub(crate) fn input_to_reader_with_format(
    input: &Option<PathBuf>,
) -> Result<(FastaReader, Format)> {
    let mut source: Box<dyn Read + Send> = match input.as_deref().filter(|path| !is_stdio(path)) {
        Some(path) => Box::new(BufReader::new(
            File::open(path).with_context(|| format!("open input {}", path.display()))?,
        )),
        None => {
            if io::stdin().is_terminal() {
                anyhow::bail!(
                    "No stdin detected. Supply an input file or pipe FASTA into circkit."
                );
            }
            Box::new(BufReader::new(io::stdin()))
        }
    };
    // niffler requires five bytes to sniff. Short and empty FASTA streams are
    // valid input too; retain the prefix and let the FASTA parser inspect them.
    let mut prefix = Vec::with_capacity(5);
    source
        .by_ref()
        .take(5)
        .read_to_end(&mut prefix)
        .context("read input prefix")?;
    let short = prefix.len() < 5;
    let source: Box<dyn Read + Send> = Box::new(io::Cursor::new(prefix).chain(source));
    let (reader, format) = if short {
        (source, Format::No)
    } else {
        niffler::send::get_reader(source).context("detect input compression")?
    };
    Ok((Reader::new(reader), format))
}

pub(crate) fn output_compression_format(output: &Option<PathBuf>) -> Format {
    match output
        .as_ref()
        .and_then(|path| path.extension())
        .and_then(|suffix| suffix.to_str())
    {
        Some("gz") => Format::Gzip,
        Some("bz2") => Format::Bzip,
        Some("xz") => Format::Lzma,
        Some("zst") => Format::Zstd,
        _ => Format::No,
    }
}

enum Sink {
    File(File),
    Stdout(StdoutLock<'static>),
    #[cfg(test)]
    Test(Box<dyn Write>),
}

impl Write for Sink {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        match self {
            Self::File(file) => file.write(bytes),
            Self::Stdout(stdout) => stdout.write(bytes),
            #[cfg(test)]
            Self::Test(writer) => writer.write(bytes),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match self {
            Self::File(file) => file.flush(),
            Self::Stdout(stdout) => stdout.flush(),
            #[cfg(test)]
            Self::Test(writer) => writer.flush(),
        }
    }
}

enum Encoder {
    Plain(Sink),
    Gzip(flate2::write::GzEncoder<Sink>),
    Bzip(bzip2::write::BzEncoder<Sink>),
    Xz(xz2::write::XzEncoder<Sink>),
    Zstd(zstd::stream::write::Encoder<'static, Sink>),
}

impl Encoder {
    fn finish(self) -> io::Result<Sink> {
        match self {
            Self::Plain(sink) => Ok(sink),
            Self::Gzip(writer) => writer.finish(),
            Self::Bzip(writer) => writer.finish(),
            Self::Xz(writer) => writer.finish(),
            Self::Zstd(writer) => writer.finish(),
        }
    }
}

impl Write for Encoder {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        match self {
            Self::Plain(writer) => writer.write(bytes),
            Self::Gzip(writer) => writer.write(bytes),
            Self::Bzip(writer) => writer.write(bytes),
            Self::Xz(writer) => writer.write(bytes),
            Self::Zstd(writer) => writer.write(bytes),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match self {
            Self::Plain(writer) => writer.flush(),
            Self::Gzip(writer) => writer.flush(),
            Self::Bzip(writer) => writer.flush(),
            Self::Xz(writer) => writer.flush(),
            Self::Zstd(writer) => writer.flush(),
        }
    }
}

#[derive(Debug)]
struct OutputError {
    destination: String,
    stdout: bool,
    source: io::Error,
}

impl std::fmt::Display for OutputError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(formatter, "write {}: {}", self.destination, self.source)
    }
}

impl std::error::Error for OutputError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(&self.source)
    }
}

fn output_error(error: io::Error, destination: &str, stdout: bool) -> io::Error {
    io::Error::new(
        error.kind(),
        OutputError {
            destination: destination.to_owned(),
            stdout,
            source: error,
        },
    )
}

/// Only stdout being closed early is a successful pipe termination. A broken
/// file/metadata pipe must still fail, even if the primary output uses stdout.
pub fn is_stdout_broken_pipe(error: &anyhow::Error) -> bool {
    error.chain().any(|source| {
        source.downcast_ref::<io::Error>().is_some_and(|error| {
            error.kind() == io::ErrorKind::BrokenPipe
                && error
                    .get_ref()
                    .and_then(|inner| inner.downcast_ref::<OutputError>())
                    .is_some_and(|output| output.stdout)
        }) || source.downcast_ref::<OutputError>().is_some_and(|output| {
            output.stdout && output.source.kind() == io::ErrorKind::BrokenPipe
        })
    })
}

pub struct OutputWriter {
    buffer: BufWriter<Encoder>,
    destination: String,
    stdout: bool,
}

impl Write for OutputWriter {
    fn write_all(&mut self, bytes: &[u8]) -> io::Result<()> {
        self.buffer
            .write_all(bytes)
            .map_err(|error| output_error(error, &self.destination, self.stdout))
    }

    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.buffer
            .write(bytes)
            .map_err(|error| output_error(error, &self.destination, self.stdout))
    }

    fn flush(&mut self) -> io::Result<()> {
        self.buffer
            .flush()
            .map_err(|error| output_error(error, &self.destination, self.stdout))
    }
}

impl OutputWriter {
    /// Flush input buffers, emit compression trailers, and report final write
    /// failures explicitly rather than relying on Drop, which discards errors.
    pub fn finish(self) -> Result<()> {
        let Self {
            buffer,
            destination,
            stdout,
        } = self;
        let result = (|| {
            let encoder = buffer.into_inner().map_err(|error| error.into_error())?;
            let mut sink = encoder.finish()?;
            sink.flush()
        })();
        result.map_err(|error| output_error(error, &destination, stdout).into())
    }
}

pub fn output_to_writer(output: &Option<PathBuf>) -> Result<OutputWriter> {
    let (sink, destination, stdout) = match output.as_deref().filter(|path| !is_stdio(path)) {
        Some(path) => (
            Sink::File(
                File::create(path).with_context(|| format!("create output {}", path.display()))?,
            ),
            path.display().to_string(),
            false,
        ),
        None => (Sink::Stdout(io::stdout().lock()), "stdout".to_owned(), true),
    };
    let encoder = match output_compression_format(output) {
        Format::No => Encoder::Plain(sink),
        Format::Gzip => Encoder::Gzip(flate2::write::GzEncoder::new(
            sink,
            flate2::Compression::default(),
        )),
        Format::Bzip => Encoder::Bzip(bzip2::write::BzEncoder::new(
            sink,
            bzip2::Compression::best(),
        )),
        Format::Lzma => Encoder::Xz(xz2::write::XzEncoder::new(sink, 6)),
        Format::Zstd => Encoder::Zstd(
            zstd::stream::write::Encoder::new(sink, 1).context("initialize zstd output")?,
        ),
    };
    Ok(OutputWriter {
        buffer: BufWriter::new(encoder),
        destination,
        stdout,
    })
}

// serde_json's Error::source skips the custom io::Error payload. Recover the
// original I/O error so destination and stdout provenance survive serialization.
pub(crate) fn json_write_error(error: serde_json::Error) -> anyhow::Error {
    if error.is_io() {
        io::Error::from(error).into()
    } else {
        error.into()
    }
}

pub enum MetadataWriter {
    Delimited(Box<csv::Writer<OutputWriter>>),
    Jsonl(Box<OutputWriter>),
}

impl MetadataWriter {
    pub fn serialize<T: Serialize>(&mut self, row: T) -> Result<()> {
        match self {
            Self::Delimited(writer) => writer
                .serialize(row)
                .map_err(|error| {
                    if error.is_io_error() {
                        match error.into_kind() {
                            csv::ErrorKind::Io(error) => anyhow::Error::from(error),
                            _ => unreachable!("is_io_error guarantees an I/O error"),
                        }
                    } else {
                        error.into()
                    }
                })
                .context("write metadata record"),
            Self::Jsonl(writer) => {
                serde_json::to_writer(&mut **writer, &row)
                    .map_err(json_write_error)
                    .context("write JSONL metadata record")?;
                writer.write_all(b"\n")?;
                Ok(())
            }
        }
    }

    pub fn finish(self) -> Result<()> {
        match self {
            Self::Delimited(writer) => (*writer)
                .into_inner()
                .map_err(|error| error.into_error())?
                .finish(),
            Self::Jsonl(writer) => (*writer).finish(),
        }
    }
}

pub fn table_path_to_writer(
    table: &Option<PathBuf>,
    format: Option<TableFormat>,
) -> Result<Option<MetadataWriter>> {
    let Some(path) = table else {
        return Ok(None);
    };
    let suffix_path = if output_compression_format(table) == Format::No {
        path.clone()
    } else {
        path.with_extension("")
    };
    let format = format.unwrap_or_else(|| {
        match suffix_path.extension().and_then(|suffix| suffix.to_str()) {
            Some("tsv") => TableFormat::Tsv,
            Some("jsonl") => TableFormat::Jsonl,
            _ => TableFormat::Csv,
        }
    });
    let output = output_to_writer(table).context("open metadata output")?;
    Ok(Some(match format {
        TableFormat::Jsonl => MetadataWriter::Jsonl(Box::new(output)),
        TableFormat::Csv | TableFormat::Tsv => MetadataWriter::Delimited(Box::new(
            csv::WriterBuilder::new()
                .delimiter(if format == TableFormat::Tsv {
                    b'\t'
                } else {
                    b','
                })
                .from_writer(output),
        )),
    }))
}

pub fn sequence_len(record: &RefRecord<'_>) -> usize {
    record.seq_lines().map(|line| line.len()).sum()
}

/// Borrow single-line records and reuse one buffer for wrapped records. Larger
/// contiguous writes are faster than sending each short FASTA line separately.
pub(crate) fn joined_sequence<'a>(record: &'a RefRecord<'_>, scratch: &'a mut Vec<u8>) -> &'a [u8] {
    if record.num_seq_lines() <= 1 {
        record.seq()
    } else {
        scratch.clear();
        for line in record.seq_lines() {
            scratch.extend_from_slice(line);
        }
        scratch
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::{cell::Cell, rc::Rc};

    struct FailingSink(Rc<Cell<bool>>);

    impl Write for FailingSink {
        fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
            if self.0.get() {
                Err(io::Error::new(io::ErrorKind::BrokenPipe, "closed sink"))
            } else {
                Ok(bytes.len())
            }
        }
        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    #[test]
    fn compression_trailer_errors_are_reported_and_pipe_origin_is_preserved() {
        for stdout in [false, true] {
            let fail = Rc::new(Cell::new(false));
            let encoder = Encoder::Gzip(flate2::write::GzEncoder::new(
                Sink::Test(Box::new(FailingSink(fail.clone()))),
                flate2::Compression::default(),
            ));
            let mut writer = OutputWriter {
                buffer: BufWriter::new(encoder),
                destination: "test".to_owned(),
                stdout,
            };
            writer.write_all(b">record\nACGT\n").unwrap();
            writer.flush().unwrap();
            fail.set(true);
            let error = writer.finish().unwrap_err();
            assert_eq!(is_stdout_broken_pipe(&error), stdout);
            assert!(error.to_string().contains("closed sink"));
        }
    }
}
