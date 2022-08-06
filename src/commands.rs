use clap::{Parser, Subcommand};
use std::path::PathBuf;

#[derive(Parser)]
#[clap(name = "circkit", author, version, about, long_about = None)]
pub struct Cli {
    #[clap(subcommand)]
    pub command: Option<Command>,
    // Level of verbosity.
    // #[clap(short, long, global = true)]
    // verbose: bool,
}
#[derive(Subcommand, Debug)]
pub enum Command {
    /// find monomers of (potentially) circular or multimeric sequences
    Monomerize {
        #[clap(short, long, action)]
        overlap_required: bool,
        /// normalize monomers
        #[clap(short, long, action)]
        normalize: bool,
    },
    /// concatenate sequences to themselves
    #[clap(visible_alias = "concat", visible_alias = "concatenate")]
    Cat {
        /// FASTA file
        input: Option<PathBuf>,
        /// Output file. If not specified, output will be written to stdout.
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
        /// FASTA file
        input: Option<PathBuf>,
        /// Output file. If not specified, output will be written to stdout.
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
        /// FASTA file
        input: Option<PathBuf>,
        /// Output file. If not specified, output will be written to stdout.
        #[clap(short, long)]
        output: Option<PathBuf>,
    },
}
