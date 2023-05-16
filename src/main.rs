use circkit_cli::{
    canonicalize::canonicalize,
    commands::{Cli, Command},
    // concatenate::{concatenate, deconcatenate},
    monomerize::monomerize,
    uniq::uniq,
};
use clap::Parser;
use env_logger::Builder;
use human_panic::setup_panic;
use log::LevelFilter;

fn main() -> anyhow::Result<()> {
    setup_panic!();

    let cli = Cli::parse();

    let mut builder = Builder::new();
    builder.filter(None, LevelFilter::Info).init();
    match &cli.command {
        Command::Monomerize { .. } => monomerize(&cli.command)?,
        Command::Cat { .. } => {
            // concatenate(input, output).expect("concatenation failed");
        }
        Command::Decat { .. } => {
            // deconcatenate(input, output).expect("deconcatenation failed");
        }
        Command::Canonicalize { .. } => {
            canonicalize(&cli.command)?;
        }
        Command::Uniq { .. } => {
            uniq(&cli.command)?;
        }
    }
    Ok(())
}
