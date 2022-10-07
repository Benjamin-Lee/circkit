use clap::{Parser, Subcommand};
use std::path::PathBuf;

#[derive(Parser)]
#[clap(name = "circkit", author, version, about, long_about = None)]
pub struct Cli {
    #[clap(subcommand)]
    pub command: Command,
    // Level of verbosity.
    // #[clap(short, long, global = true)]
    // verbose: bool,
}

#[derive(Subcommand, Debug)]
pub enum Command {
    /// Find monomers of (potentially) circular or multimeric sequences
    Monomerize {
        /// Input FASTA file. May be gzip, bzip, or xz compressed [default: stdin]
        input: Option<PathBuf>,
        #[clap(short, long)]
        /// Output FASTA file path [default: stdout]
        output: Option<PathBuf>,
        /// Whether to canonicalize the monomers. This is faster than normalizing separately since it skips reading the sequences back into memory
        #[clap(short, long)]
        normalize: bool,
        /// The length of the seed to search for. Must be less than or equal to the length of the sequence but should be much smaller to be meaningful
        #[clap(long, default_value = "10", value_parser = clap::value_parser!(u64).range(5..=64))]
        seed_length: u64,
        #[clap(long, group = "max_mismatch")]
        /// The maximum number of mismatches to allow in the overlap. Conflicts with --min-identity
        max_mismatch: Option<u64>,
        /// The maximum percentage of the overlap that can be mismatched. Conflicts with --max-mismatch
        #[clap(long, conflicts_with = "max_mismatch")]
        min_identity: Option<f64>,

        /// Minimum length of the overlap (in nt) to keep the monomer. If the overlap is less than this, the monomer is discarded
        #[clap(long)]
        min_overlap: Option<u64>,

        /// Minimum length of the overlap (relative to the input sequence) to require
        #[clap(long)]
        min_overlap_percent: Option<f64>,

        /// Whether to output monomers that did not have any overlap. These may not be circular or multimeric
        #[clap(short, long)]
        keep_all: bool,
    },
    /// Find monomers of (potentially) circular or multimeric sequences
    Monomerize2 {
        /// Input FASTA file. May be gzip, bzip, or xz compressed [default: stdin]
        input: Option<PathBuf>,
        #[clap(short, long)]
        /// Output FASTA file path [default: stdout]
        output: Option<PathBuf>,
        /// Whether to canonicalize the monomers.
        /// This is faster than normalizing separately (perhaps via piping) since it skips reading the sequences back into memory.
        /// Because this is a common use case, it is available as an option.
        /// See also --uniq
        #[clap(short, long, alias = "norm", alias = "canonicalize", alias = "canon")]
        normalize: bool,
        /// The length of the seed to search for. Must be less than or equal to the length of the sequence but should be much smaller to be meaningful
        #[clap(long, default_value = "10", value_parser = clap::value_parser!(u64).range(5..=64))]
        seed_length: u64,

        // Overlap similarity cutoffs
        #[clap(long, group = "overlap_cutoffs")]
        /// The maximum number of mismatches to allow in the overlap. Conflicts with --min-identity
        max_mismatch: Option<u64>,
        /// The maximum percentage of the overlap that can be mismatched. Conflicts with --max-mismatch
        #[clap(long, conflicts_with = "overlap_cutoffs")]
        min_identity: Option<f64>,

        /// Minimum length of the overlap (in nt) to keep the monomer. If the overlap is less than this, the monomer is discarded
        #[clap(long)]
        min_overlap: Option<u64>,

        /// Minimum length of the overlap (relative to the input sequence) to require
        #[clap(long)]
        min_overlap_percent: Option<f64>,

        /// Whether to output sequences that did not have any overlap.
        /// These sequences could possibly be circular or multimeric since they failed to monomerize.
        #[clap(short, long)]
        keep_all: bool,

        /// The number of threads to use. If not specified, the number of logical cores is used.
        #[clap(short, long, default_value_t = num_cpus::get())]
        threads: usize,
    },
    /// concatenate sequences to themselves
    #[clap(visible_alias = "concat", visible_alias = "concatenate")]
    Cat {
        /// Input FASTA file. May be gzip, bzip, or xz compressed [default: stdin]
        input: Option<PathBuf>,
        /// Output FASTA file path [default: stdout]
        #[clap(short, long)]
        output: Option<PathBuf>,
    },

    /// deconcatenate sequences to themselves
    /// (this is the reverse of `cat`)
    #[clap(
        visible_alias = "deconcat",
        visible_alias = "deconcatenate",
        visible_alias = "uncat",
        visible_alias = "unconcatenate"
    )]
    Decat {
        /// Input FASTA file. May be gzip, bzip, or xz compressed [default: stdin]
        input: Option<PathBuf>,
        /// Output FASTA file path [default: stdout]
        #[clap(short, long)]
        output: Option<PathBuf>,
    },

    /// Normalize circular sequences.
    #[clap(
        alias = "rotcanon",
        visible_alias = "canonicalize",
        visible_alias = "canon"
    )]
    Normalize {
        /// Input FASTA file. May be gzip, bzip, or xz compressed [default: stdin]
        input: Option<PathBuf>,
        /// Output FASTA file path [default: stdout]
        #[clap(short, long)]
        output: Option<PathBuf>,
    },
    Normalize2 {
        /// Input FASTA file. May be gzip, bzip, or xz compressed [default: stdin]
        input: Option<PathBuf>,
        /// Output FASTA file path [default: stdout]
        #[clap(short, long)]
        output: Option<PathBuf>,
    },
    /// Deduplicate circular sequences
    Uniq {
        /// Input FASTA file. May be gzip, bzip, or xz compressed [default: stdin]
        input: Option<PathBuf>,
        /// Output FASTA file path [default: stdout]
        #[clap(short, long)]
        output: Option<PathBuf>,
        /// Whether output normalized circular sequences.
        /// This is faster than normalizing separately (perhaps via piping) since the sequences are normalized anyway when deduplicating.
        #[clap(short, long, alias = "norm", alias = "canonicalize", alias = "canon")]
        normalize: bool,
    },
}
