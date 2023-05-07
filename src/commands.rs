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
        /// The length of the seed to search for. Must be less than or equal to the length of the sequence but should be much smaller to be meaningful
        #[clap(long, default_value = "10", value_parser = clap::value_parser!(u64).range(5..=64))]
        seed_length: u64,

        // Overlap similarity cutoffs
        #[clap(long, group = "overlap_cutoffs")]
        /// The maximum number of mismatches to allow in the overlap. Conflicts with --min-identity
        max_mismatch: Option<u64>,
        /// The minimum identity the overlapping region before being considered mismatched. Conflicts with --max-mismatch
        #[clap(long, conflicts_with = "overlap_cutoffs")]
        min_identity: Option<f64>,

        /// Minimum length of the overlap (in nt) required to keep the monomer. If the overlap is shorter than this, the monomer is discarded unless --keep-all is used, in which case the original sequence (without trimming) is output. Can be combined with --min-overlap-percent for more stringent filtering.
        #[clap(long)]
        min_overlap: Option<usize>,

        /// Minimum length of the overlap (relative to the input sequence) to require. Can be used with --min-overlap for more stringent filtering. If --keep-all is used, sequences with too short of an overlap are still output but as the original sequence.
        #[clap(long)]
        min_overlap_percent: Option<f64>,

        /// Whether to output sequences that did not have any overlap.
        /// These sequences could possibly be circular or multimeric since they failed to monomerize.
        /// Useful for cleaning up  datasets in which the sequences are not all monomers (e.g. viroids in GenBank).
        #[clap(short, long)]
        keep_all: bool,

        /// The number of threads to use. If not specified, the number of logical cores is used.
        #[clap(short, long, default_value_t = num_cpus::get().try_into().unwrap())]
        threads: u32,

        #[clap(long, hidden = true, default_value_t = 64)]
        batch_size: usize,
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
        /// The number of threads to use. If not specified, the number of logical cores is used.
        #[clap(short, long, default_value_t = num_cpus::get().try_into().unwrap())]
        threads: u32,
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
        /// The number of threads to use. If not specified, the number of logical cores is used.
        #[clap(short, long, default_value_t = num_cpus::get().try_into().unwrap())]
        threads: u32,
    },
}
