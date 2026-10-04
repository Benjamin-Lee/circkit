use crate::commands::{Execution, ProcessingOptions};
use crate::diagnostics::argument_error;
pub use crate::io::{input_to_reader, output_to_writer, table_path_to_writer, FastaReader};
use crate::io::{input_to_reader_with_format, output_compression_format};
use niffler::send::compression::Format;
use seq_io::fasta::Reader;
use std::borrow::Cow;
use std::{io::Read, path::PathBuf};

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
            return Err(argument_error("The number of threads must be at least one"));
        }
        if self.execution == Execution::Serial && threads != 1 {
            return Err(argument_error("Serial execution requires --threads 1"));
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
pub fn process_fasta<D, W, F>(
    mut reader: Reader<Box<dyn Read + Send>>,
    settings: FastaSettings,
    work: W,
    mut output: F,
) -> anyhow::Result<()>
where
    D: Default + Send,
    W: Send + Sync + Fn(seq_io::fasta::RefRecord<'_>, &mut D),
    F: FnMut(seq_io::fasta::RefRecord<'_>, &mut D) -> anyhow::Result<()>,
{
    if settings.serial {
        let mut data = D::default();
        while let Some(record) = reader.next() {
            let record = record?;
            work(record.clone(), &mut data);
            output(record, &mut data)?;
        }
        Ok(())
    } else {
        let error = seq_io::parallel::parallel_fasta(
            reader,
            settings.threads,
            settings.queue_depth,
            work,
            |record, data| output(record, data).err(),
        )?;
        match error {
            Some(error) => Err(error),
            None => Ok(()),
        }
    }
}

/// Plain single-worker canonicalization benefits from reusing buffers without
/// pipeline coordination. Compressed I/O keeps the existing overlap policy.
pub fn canonicalization_input(
    input: &Option<PathBuf>,
    output: &Option<PathBuf>,
    processing: &ProcessingOptions,
    threads: u32,
) -> anyhow::Result<(FastaReader, FastaSettings)> {
    // Validate execution options before opening a file or waiting for stdin.
    let mut settings = processing.resolve(threads, 64)?;
    let (reader, format) = input_to_reader_with_format(input)?;
    if processing.execution == Execution::Auto
        && threads == 1
        && format == Format::No
        && output_compression_format(output) == Format::No
    {
        settings.serial = true;
    }
    Ok((reader, settings))
}
