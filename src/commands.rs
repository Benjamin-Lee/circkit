use clap::{Args, Parser, Subcommand};
use std::{num::NonZeroUsize, path::PathBuf};

use crate::orfs::Strand;

pub fn default_threads() -> u32 {
    std::thread::available_parallelism().map_or(1, |cpus| cpus.get().min(u32::MAX as usize) as u32)
}

#[derive(clap::ValueEnum, Clone, Copy, Debug, PartialEq, Eq)]
pub enum ErrorFormat {
    Text,
    Json,
}

#[derive(clap::ValueEnum, Clone, Copy, Debug, PartialEq, Eq)]
pub enum TableFormat {
    Csv,
    Tsv,
    Jsonl,
}

fn parse_finite(value: &str) -> Result<f64, String> {
    value
        .parse::<f64>()
        .map_err(|_| "expected a number".to_owned())
        .and_then(|number| {
            if number.is_finite() {
                Ok(number)
            } else {
                Err("must be finite".to_owned())
            }
        })
}

fn parse_fraction(value: &str) -> Result<f64, String> {
    let number = parse_finite(value)?;
    if (0.0..=1.0).contains(&number) {
        Ok(number)
    } else {
        Err("must be between 0 and 1".to_owned())
    }
}

fn parse_nonnegative(value: &str) -> Result<f64, String> {
    let number = parse_finite(value)?;
    if number >= 0.0 {
        Ok(number)
    } else {
        Err("must be nonnegative".to_owned())
    }
}

fn parse_bases(value: &str) -> Result<i64, String> {
    let number = value
        .parse::<i64>()
        .map_err(|_| "expected a signed integer".to_owned())?;
    if number != 0 {
        Ok(number)
    } else {
        Err("rotation by zero is not allowed".to_owned())
    }
}

fn parse_percent(value: &str) -> Result<f64, String> {
    let number = parse_finite(value)?;
    if number != 0.0 {
        Ok(number)
    } else {
        Err("rotation by zero is not allowed".to_owned())
    }
}

fn parse_wraps(value: &str) -> Result<usize, String> {
    let number = value
        .parse::<usize>()
        .map_err(|_| "expected an integer".to_owned())?;
    if number <= 3 {
        Ok(number)
    } else {
        Err("must be between 0 and 3".to_owned())
    }
}

fn parse_codons(value: &str) -> Result<String, String> {
    let mut codons = Vec::new();
    for codon in value.split(',') {
        let codon = codon.trim().to_ascii_uppercase();
        if codon.len() != 3
            || !codon
                .bytes()
                .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T' | b'N'))
        {
            return Err("expected comma-separated DNA triplets using A, C, G, T, or N".to_owned());
        }
        codons.push(codon);
    }
    Ok(codons.join(","))
}

#[derive(clap::ValueEnum, Clone, Copy, Debug, PartialEq)]
pub enum Execution {
    Auto,
    Serial,
    Pipeline,
}

impl std::fmt::Display for Execution {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::Auto => "auto",
            Self::Serial => "serial",
            Self::Pipeline => "pipeline",
        })
    }
}

#[derive(Args, Debug)]
pub struct ProcessingOptions {
    /// FASTA buffers queued in the pipeline. Defaults to 64, or max(2, 2*threads) for orfs.
    #[arg(long, alias = "batch-size")]
    pub queue_depth: Option<NonZeroUsize>,

    /// Execution strategy. Auto uses serial processing for one worker on one available CPU,
    /// or for single-worker canonicalize/uniq with uncompressed input and output.
    /// Serial requires --threads 1; pipeline overlaps reading, processing, and writing.
    #[arg(long, value_enum, default_value_t = Execution::Auto)]
    pub execution: Execution,
}

#[derive(Args, Debug)]
pub struct RotationOptions {
    /// Maximum sequence length for allocating Duval search. Zero always uses constant-space search.
    #[arg(long, default_value_t = circkit::canonicalize::RotationOptions::default().duval_max_len)]
    pub rotation_cutoff: usize,
}

#[derive(Parser)]
#[command(name = "circkit", author, version, about, long_about = None,
    after_help = "Examples:\n  circkit canonicalize reads.fasta -o canonical.fasta\n  circkit monomerize - --keep-all --threads 1\n  circkit schema orfs\n\nUse '-' for stdin/stdout. Run 'circkit COMMAND --help' for options.")]
pub struct Cli {
    #[command(subcommand)]
    pub command: Command,
    // Level of verbosity.
    #[command(flatten)]
    pub verbose: clap_verbosity_flag::Verbosity,

    /// Diagnostic format on stderr; FASTA output is unaffected.
    #[arg(long, global = true, value_enum, default_value_t = ErrorFormat::Text)]
    pub error_format: ErrorFormat,
}

#[derive(Subcommand, Debug)]
pub enum Command {
    /// Find monomers of (potentially) circular or multimeric sequences
    Monomerize {
        /// Input FASTA file, or '-' for stdin. Compression detected from bytes [default: stdin]
        #[arg(value_hint = clap::ValueHint::FilePath)]
        input: Option<PathBuf>,

        #[arg(short, long, value_hint = clap::ValueHint::FilePath)]
        /// Output FASTA file, or '-' for stdout. .gz/.bz2/.xz/.zst selects compression [default: stdout]
        output: Option<PathBuf>,

        #[arg(long)]
        /// Whether to check the sequence in reverse when the forward pass monomerization is complete.
        /// This mode can handle mutations in the seed.
        /// Using this flag will roughly double the runtime, since each sequence must now be processed twice.
        sensitive: bool,

        /// The length of the seed to search for (5 to 64 bases).
        /// Must be less than or equal to the length of the sequence but should be much smaller to be meaningful
        #[arg(long, default_value = "10", value_parser = clap::value_parser!(u64).range(5..=64))]
        seed_length: u64,

        // Overlap similarity cutoffs
        #[arg(long, group = "overlap_cutoffs")]
        /// The maximum number of mismatches to allow in the overlap.
        /// Conflicts with --min-identity
        max_mismatch: Option<u64>,

        /// The minimum identity of the overlapping region (a fraction between 0 and 1).
        /// Conflicts with --max-mismatch
        #[arg(long, conflicts_with = "overlap_cutoffs", value_parser = parse_fraction)]
        min_identity: Option<f64>,

        /// Minimum length of the overlap (in nt) required to keep the monomer.
        /// If the overlap is shorter than this, the monomer is discarded unless --keep-all is used, in which case the original sequence (without trimming) is output.
        /// Can be combined with --min-overlap-percent for more stringent filtering.
        #[arg(long)]
        min_overlap: Option<usize>,

        /// Minimum overlap relative to the input sequence, as a finite nonnegative ratio.
        /// A value of 1.0 means that the sequence must be a complete dimer.
        /// Can be used with --min-overlap for more stringent filtering.
        /// If --keep-all is used, sequences with too short of an overlap are still output but as the original sequence.
        #[arg(long, visible_alias = "min-overlap-ratio", value_parser = parse_nonnegative)]
        min_overlap_percent: Option<f64>,

        /// The minimum length of the monomer to keep (in nt).
        #[arg(long, default_value_t = 0)]
        min_length: usize,

        /// The maximum length of the monomer to keep (in nt).
        #[arg(long)]
        max_length: Option<usize>,

        /// Whether to output sequences that did not have any overlap.
        /// These sequences could possibly be circular or multimeric since they failed to monomerize.
        /// Useful for cleaning up datasets in which the sequences are not all monomers (e.g. viroids in GenBank).
        #[arg(short, long)]
        keep_all: bool,

        /// A path for the monomerization metadata for each sequence.
        /// The following columns are output: id, original_length, monomer_length.
        /// Format follows the extension: .tsv, .jsonl, otherwise CSV; compression suffixes are supported.
        /// Use '-' for metadata on stdout with -o FILE.
        /// Note that if no sequences are output, the output table will be an empty file.
        #[arg(long, value_hint = clap::ValueHint::FilePath)]
        table: Option<PathBuf>,

        /// Metadata format; requires --table. Default: .tsv selects TSV, .jsonl selects JSONL, otherwise CSV.
        #[arg(long, value_enum, requires = "table")]
        table_format: Option<TableFormat>,

        /// The number of threads to use.
        /// Defaults to available logical CPUs, respecting OS affinity and CPU quotas.
        #[arg(short, long, default_value_t = default_threads(), value_parser = clap::value_parser!(u32).range(1..))]
        threads: u32,

        /// Bytes checked between early mismatch rejections. Does not change matching criteria.
        #[arg(long, default_value_t = circkit::Monomerizer::default().mismatch_chunk_size)]
        mismatch_chunk_size: NonZeroUsize,

        #[command(flatten)]
        processing: ProcessingOptions,
    },
    /// Concatenate sequences to themselves
    Cat {
        /// Input FASTA file, or '-' for stdin. Compression detected from bytes [default: stdin]
        #[arg(value_hint = clap::ValueHint::FilePath)]
        input: Option<PathBuf>,

        /// Output FASTA file, or '-' for stdout. .gz/.bz2/.xz/.zst selects compression [default: stdout]
        #[arg(short, long, value_hint = clap::ValueHint::FilePath)]
        output: Option<PathBuf>,
    },

    /// Deconcatenate sequences to themselves
    Decat {
        /// Input FASTA file, or '-' for stdin. Compression detected from bytes [default: stdin]
        #[arg(value_hint = clap::ValueHint::FilePath)]
        input: Option<PathBuf>,

        /// Output FASTA file, or '-' for stdout. .gz/.bz2/.xz/.zst selects compression [default: stdout]
        #[arg(short, long, value_hint = clap::ValueHint::FilePath)]
        output: Option<PathBuf>,
    },

    /// Select a consistent rotation and strand for circular sequences.
    #[command(visible_alias = "canon")]
    Canonicalize {
        /// Input FASTA file, or '-' for stdin. Compression detected from bytes [default: stdin]
        #[arg(value_hint = clap::ValueHint::FilePath)]
        input: Option<PathBuf>,

        /// Output FASTA file, or '-' for stdout. .gz/.bz2/.xz/.zst selects compression [default: stdout]
        #[arg(short, long, value_hint = clap::ValueHint::FilePath)]
        output: Option<PathBuf>,

        /// The number of threads to use.
        /// Defaults to available logical CPUs, respecting OS affinity and CPU quotas.
        #[arg(short, long, default_value_t = default_threads(), value_parser = clap::value_parser!(u32).range(1..))]
        threads: u32,

        #[command(flatten)]
        processing: ProcessingOptions,

        #[command(flatten)]
        rotation: RotationOptions,
    },
    /// Deduplicate circular sequences
    Uniq {
        /// Input FASTA file, or '-' for stdin. Compression detected from bytes [default: stdin]
        #[arg(value_hint = clap::ValueHint::FilePath)]
        input: Option<PathBuf>,

        /// Output FASTA file, or '-' for stdout. .gz/.bz2/.xz/.zst selects compression [default: stdout]
        #[arg(short, long, value_hint = clap::ValueHint::FilePath)]
        output: Option<PathBuf>,

        /// Whether output canonicalized circular sequences.
        /// This is faster than canonicalizing separately (perhaps via piping) since the sequences are canonicalized anyway when deduplicating.
        #[arg(short, long, alias = "norm", alias = "canonicalize", alias = "canon")]
        canonicalize: bool,

        /// A path for deduplication metadata.
        /// The following columns are output: id, duplicate_id.
        /// Format follows the extension: .tsv, .jsonl, otherwise CSV; compression suffixes are supported.
        /// Use '-' for metadata on stdout with -o FILE.
        /// Note that if no sequences are output, the output table will be an empty file.
        #[arg(long, value_hint = clap::ValueHint::FilePath)]
        table: Option<PathBuf>,

        /// Metadata format; requires --table. Default: .tsv selects TSV, .jsonl selects JSONL, otherwise CSV.
        #[arg(long, value_enum, requires = "table")]
        table_format: Option<TableFormat>,

        /// The number of threads to use. Defaults to available logical CPUs, respecting OS affinity and CPU quotas.
        #[arg(short, long, default_value_t = default_threads(), value_parser = clap::value_parser!(u32).range(1..))]
        threads: u32,

        #[command(flatten)]
        processing: ProcessingOptions,

        #[command(flatten)]
        rotation: RotationOptions,
    },

    /// Rotate circular sequences to the left or right
    #[command(group(clap::ArgGroup::new("rotation").required(true).args(["bases", "percent"])))]
    Rotate {
        /// Input FASTA file, or '-' for stdin. Compression detected from bytes [default: stdin]
        #[arg(value_hint = clap::ValueHint::FilePath)]
        input: Option<PathBuf>,

        /// Output FASTA file, or '-' for stdout. .gz/.bz2/.xz/.zst selects compression [default: stdout]
        #[arg(short, long, value_hint = clap::ValueHint::FilePath)]
        output: Option<PathBuf>,

        /// The nonzero number of bases to rotate the sequence. Positive numbers rotate to the right, negative numbers rotate to the left.
        /// Rotation by amounts greater than the sequence length are equivalent to rotation by the remainder of the division of the rotation amount by the sequence length.
        /// For example, rotating a sequence of length 100 by 101 bases is equivalent to rotating by 1 base.
        /// This flag is mutually exclusive with --percent.
        #[arg(short, long, allow_hyphen_values = true, value_parser = parse_bases,
            conflicts_with = "percent")]
        bases: Option<i64>,

        /// The finite, nonzero fraction of the sequence to rotate.
        /// This must be expressed as a decimal, e.g. 0.5 for 50%.
        /// This flag is mutually exclusive with --bases.
        #[arg(short, long, allow_hyphen_values = true, value_parser = parse_percent,
            conflicts_with = "bases")]
        percent: Option<f64>,
    },

    /// Find ORFs in circular sequences
    Orfs {
        /// Input FASTA file, or '-' for stdin. Compression detected from bytes [default: stdin]
        #[arg(value_hint = clap::ValueHint::FilePath)]
        input: Option<PathBuf>,

        /// Output FASTA file, or '-' for stdout. .gz/.bz2/.xz/.zst selects compression [default: stdout]
        #[arg(short, long, value_hint = clap::ValueHint::FilePath)]
        output: Option<PathBuf>,

        /// Minimum coding length in nucleotides, excluding the stop codon.
        #[arg(short, long, default_value = "75")]
        min_length: usize,

        /// The start codons to use.
        /// Comma-separated DNA triplets, e.g. "ATG,GTG"; case-insensitive, N matches literal N.
        #[arg(long, default_value = "ATG", value_parser = parse_codons)]
        start_codons: String,

        /// The stop codons to use.
        /// Comma-separated DNA triplets, e.g. "TAA,TAG,TGA"; case-insensitive, N matches literal N.
        #[arg(long, default_value = "TAA,TAG,TGA", value_parser = parse_codons)]
        stop_codons: String,

        /// Whether to include the stop codon in the output sequence
        #[arg(long, action)]
        include_stop: bool,

        /// Whether to require a stop codon in the ORF. Required by default.
        /// If enabled, partial ORFs are allowed (e.g. ATG AAA GTC)
        #[arg(long, action)]
        no_stop_required: bool,

        /// When present, the minimum number of wraps around the origin an ORF must have in order to be output.
        /// Values greater than 0 mean that ORFs must take advantage or sequence circularity.
        #[arg(long, default_value = "0", value_parser = parse_wraps)]
        min_wraps: usize,

        /// When present, the maximum number of wraps around the origin an ORF can have in order to be output.
        /// Setting this to 0 means that this function acts as a traditional ORF finder.
        /// The most possible wraps is 3.
        #[arg(long, default_value = "3", value_parser = parse_wraps)]
        max_wraps: usize,

        /// The strands in which to search for ORFs
        #[arg(long, value_enum, default_value_t = Strand::Both)]
        strand: Strand,

        /// The minimum ORF length to sequence length ratio to keep (finite and nonnegative).
        /// A ratio of 1 means that the ORF is as long as the sequence.
        /// A ratio of 2 means that the ORF would wrap around the origin twice.
        /// The stop codon is included in the length calculation regardless of the --include-stop flag.
        #[arg(long, default_value = "0", value_parser = parse_nonnegative)]
        min_ratio: f64,

        /// A path for the ORF-finding metadata for each sequence.
        /// The following columns are output: orf_id, seq_id, start, stop, wraps, length, and ratio.
        /// Both start and stop are 0-indexed.
        /// When --no-stop-required is used, the stop column may be empty.
        /// Note that the length is the length of the ORF and not the length of the sequence.
        /// The --include-stop flag is taken into account when calculating the length.
        /// The wraps field corresponds to the number of wraps around the origin.
        /// The ratio field is the ratio of the ORF length (including the stop codon regardless of --include-stop) to the sequence length.
        /// Format follows the extension: .tsv, .jsonl, otherwise CSV; compression suffixes are supported.
        /// Use '-' for metadata on stdout with -o FILE.
        /// Note that if no sequences are output, the output table will be an empty file.
        #[arg(long, value_hint = clap::ValueHint::FilePath)]
        table: Option<PathBuf>,

        /// Metadata format; requires --table. Default: .tsv selects TSV, .jsonl selects JSONL, otherwise CSV.
        #[arg(long, value_enum, requires = "table")]
        table_format: Option<TableFormat>,

        /// The number of threads to use.
        /// Defaults to available logical CPUs, respecting OS affinity and CPU quotas.
        #[arg(short, long, default_value_t = default_threads(), value_parser = clap::value_parser!(u32).range(1..))]
        threads: u32,

        #[command(flatten)]
        processing: ProcessingOptions,
    },

    /// Print a versioned JSON command catalog for automation and agents.
    #[command(visible_alias = "describe")]
    Schema {
        /// Restrict the catalog to one command or alias; omit to describe every command.
        command: Option<String>,
    },

    /// Generate shell completion scripts on stdout.
    Completions {
        #[arg(value_enum)]
        shell: clap_complete::Shell,
    },
}

impl Command {
    pub fn validate(&self) -> anyhow::Result<()> {
        use crate::diagnostics::argument_error;
        let paths = match self {
            Self::Monomerize {
                input,
                output,
                table,
                min_length,
                max_length,
                ..
            } => {
                if max_length.is_some_and(|maximum| *min_length > maximum) {
                    return Err(argument_error("--min-length cannot exceed --max-length"));
                }
                Some((input, output, Some(table)))
            }
            Self::Orfs {
                input,
                output,
                table,
                min_wraps,
                max_wraps,
                ..
            } => {
                if min_wraps > max_wraps {
                    return Err(argument_error("--min-wraps cannot exceed --max-wraps"));
                }
                Some((input, output, Some(table)))
            }
            Self::Uniq {
                input,
                output,
                table,
                ..
            } => Some((input, output, Some(table))),
            Self::Cat { input, output }
            | Self::Decat { input, output }
            | Self::Canonicalize { input, output, .. }
            | Self::Rotate { input, output, .. } => Some((input, output, None)),
            Self::Schema { .. } | Self::Completions { .. } => None,
        };
        if let Some((input, output, table)) = paths {
            crate::io::validate_paths(input, output, table.and_then(|table| table.as_ref()))?;
        }
        Ok(())
    }
}
