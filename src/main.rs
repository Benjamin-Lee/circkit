use circkit_cli::{
    commands::{Cli, Command},
    concatenate::{concatenate, deconcatenate},
    normalize::normalize,
    utils::get_reader,
};

use clap::{Parser, Subcommand};
use env_logger::Builder;
use human_panic::setup_panic;
use log::LevelFilter;
use std::path::PathBuf;

fn main() -> anyhow::Result<()> {
    setup_panic!();

    let cli = Cli::parse();

    let mut builder = Builder::new();
    builder.filter(None, LevelFilter::Info).init();

    match &cli.command {
        Some(Command::Monomerize { .. }) => (),
        Some(Command::Cat { input, output }) => {
            concatenate(input, output).expect("concatenation failed");
        }
        Some(Command::Decat { input, output }) => {
            deconcatenate(input, output).expect("deconcatenation failed");
        }
        Some(Command::Normalize { .. }) => {
            normalize(&cli.command.expect("reqs"))?;
            // normalize(get_reader(input), output).expect("normalization failed");
        }
        None => {}
    }
    Ok(())
}
