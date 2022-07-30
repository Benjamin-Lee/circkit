use circkit::concatenate::{concatenate, deconcatenate};
use clap::{Parser, Subcommand};
use env_logger::Builder;
use human_panic::setup_panic;
use log::LevelFilter;
use std::path::PathBuf;

#[derive(Parser)]
#[clap(author, version, about, long_about = None)]
struct Cli {
    #[clap(subcommand)]
    command: Option<Commands>,
    // Level of verbosity.
    // #[clap(short, long, global = true)]
    // verbose: bool,
}
#[derive(Subcommand)]
enum Commands {
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
}

fn main() {
    setup_panic!();

    let cli = Cli::parse();

    let mut builder = Builder::new();

    builder.filter(None, LevelFilter::Info).init();

    match &cli.command {
        Some(Commands::Monomerize { .. }) => (),
        Some(Commands::Cat { input, output }) => {
            concatenate(input, output).unwrap();
        }
        Some(Commands::Decat { input, output }) => {
            deconcatenate(input, output).expect("deconcatenation failed");
        }
        None => {}
    }
}
