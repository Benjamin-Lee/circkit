use circkit_cli::{
    canonicalize::canonicalize,
    commands::{Cli, Command},
    concatenate::{concatenate, deconcatenate},
    monomerize::monomerize,
    orfs::orfs,
    rotate::rotate,
    uniq::uniq,
};
use clap::Parser;
use human_panic::setup_panic;

fn main() -> anyhow::Result<()> {
    setup_panic!();

    let cli = Cli::parse();

    env_logger::Builder::new()
        .filter_level(cli.verbose.log_level_filter())
        .init();

    match &cli.command {
        Command::Monomerize { .. } => monomerize(&cli.command)?,
        Command::Cat { .. } => {
            concatenate(&cli.command)?;
        }
        Command::Decat { .. } => {
            deconcatenate(&cli.command)?;
        }
        Command::Canonicalize { .. } => {
            canonicalize(&cli.command)?;
        }
        Command::Uniq { .. } => {
            uniq(&cli.command)?;
        }
        Command::Rotate { .. } => rotate(&cli.command)?,
        Command::Orfs { .. } => orfs(&cli.command)?,
    }
    Ok(())
}
