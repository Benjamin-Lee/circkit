use crate::commands::{Execution, ProcessingOptions};
use anyhow::bail;
use seq_io::fasta::Reader;
use std::borrow::Cow;
use std::{
    fs::File,
    io::{prelude::*, stdin, stdout, BufReader, BufWriter},
    path::PathBuf,
};

/// Borrow already normalized DNA instead of allocating and then discarding a copy.
pub fn normalized_sequence(sequence: &[u8]) -> Cow<'_, [u8]> {
    if already_normalized(sequence) {
        Cow::Borrowed(sequence)
    } else {
        Cow::Owned(
            needletail::sequence::normalize(sequence, false).unwrap_or_else(|| sequence.to_vec()),
        )
    }
}

#[inline]
fn already_normalized(sequence: &[u8]) -> bool {
    #[cfg(all(
        any(target_arch = "x86", target_arch = "x86_64"),
        not(feature = "scalar-normalization")
    ))]
    if std::is_x86_feature_detected!("avx2") && sequence.len() >= 32 {
        // SAFETY: feature detection establishes AVX2 support; the kernel loads
        // only complete 32-byte chunks, including an overlapping load for the tail.
        return unsafe { normalized_avx2(sequence) };
    }
    sequence
        .iter()
        .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'N' | b'-'))
}

#[cfg(all(
    any(target_arch = "x86", target_arch = "x86_64"),
    not(feature = "scalar-normalization")
))]
#[target_feature(enable = "avx2")]
unsafe fn normalized_avx2(sequence: &[u8]) -> bool {
    #[cfg(target_arch = "x86")]
    use std::arch::x86::*;
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::*;

    // Accepted bytes have distinct low nibbles. A shuffle supplies the expected
    // byte for each nibble; an exact comparison rejects every other byte.
    // Unused entries are 0xff, and the shuffle zeroes lanes with bit 7 set.
    let lookup = _mm256_broadcastsi128_si256(_mm_setr_epi8(
        -1, b'A' as i8, -1, b'C' as i8, b'T' as i8, -1, -1, b'G' as i8, -1, -1, -1, -1, -1,
        b'-' as i8, b'N' as i8, -1,
    ));
    let check = |bytes| {
        _mm256_movemask_epi8(_mm256_cmpeq_epi8(_mm256_shuffle_epi8(lookup, bytes), bytes)) == -1
    };
    let mut offset = 0;
    while offset + 32 <= sequence.len() {
        let bytes = _mm256_loadu_si256(sequence[offset..].as_ptr().cast());
        if !check(bytes) {
            return false;
        }
        offset += 32;
    }
    if offset == sequence.len() {
        return true;
    }
    // The overlapping final load is entirely within this >=32-byte slice.
    check(_mm256_loadu_si256(
        sequence[sequence.len() - 32..].as_ptr().cast(),
    ))
}

pub struct FastaSettings {
    pub threads: u32,
    pub queue_depth: usize,
    pub serial: bool,
}

impl ProcessingOptions {
    pub fn resolve(
        &self,
        threads: u32,
        default_queue_depth: usize,
    ) -> anyhow::Result<FastaSettings> {
        if threads == 0 {
            bail!("The number of threads must be at least one");
        }
        if self.execution == Execution::Serial && threads != 1 {
            bail!("Serial execution requires --threads 1");
        }
        let single_cpu = std::thread::available_parallelism().map_or(true, |cpus| cpus.get() == 1);
        Ok(FastaSettings {
            threads,
            queue_depth: self.queue_depth.map_or(default_queue_depth, |n| n.get()),
            serial: self.execution == Execution::Serial
                || (self.execution == Execution::Auto && threads == 1 && single_cpu),
        })
    }
}

/// Serial execution reuses one record's work buffers. Pipeline execution retains
/// overlap between reading, processing, and writing, with bounded queues.
pub fn process_fasta<D, W, F, Out>(
    mut reader: Reader<Box<dyn Read + Send>>,
    settings: FastaSettings,
    work: W,
    mut output: F,
) -> anyhow::Result<Option<Out>>
where
    D: Default + Send,
    W: Send + Sync + Fn(seq_io::fasta::RefRecord<'_>, &mut D),
    F: FnMut(seq_io::fasta::RefRecord<'_>, &mut D) -> Option<Out>,
{
    if settings.serial {
        let mut data = D::default();
        while let Some(record) = reader.next() {
            let record = record?;
            work(record.clone(), &mut data);
            if let Some(value) = output(record, &mut data) {
                return Ok(Some(value));
            }
        }
        Ok(None)
    } else {
        Ok(seq_io::parallel::parallel_fasta(
            reader,
            settings.threads,
            settings.queue_depth,
            work,
            output,
        )?)
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
