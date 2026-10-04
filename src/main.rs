use circkit_cli::{
    canonicalize::canonicalize,
    commands::{Cli, Command},
    concatenate::{concatenate, deconcatenate},
    monomerize::monomerize,
    orfs::orfs,
    rotate::rotate,
    uniq::uniq,
};
use circkit_cli::{commands::ErrorFormat, diagnostics, discovery};
use clap::Parser;
use std::process::ExitCode;

fn main() -> ExitCode {
    let arguments: Vec<_> = std::env::args_os().collect();
    let format = diagnostics::requested_format(&arguments);
    let cli = match Cli::try_parse_from(arguments) {
        Ok(cli) => cli,
        Err(error) => {
            let exit_code = error.exit_code() as u8;
            if exit_code == 0 || format == ErrorFormat::Text {
                let _ = error.print();
            } else {
                diagnostics::emit(format, "invalid_arguments", &error.to_string(), exit_code);
            }
            return ExitCode::from(exit_code);
        }
    };

    env_logger::Builder::new()
        .filter_level(cli.verbose.log_level_filter())
        .init();

    let result = cli.command.validate().and_then(|()| run(&cli.command));
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) if circkit_cli::io::is_stdout_broken_pipe(&error) => ExitCode::SUCCESS,
        Err(error) => ExitCode::from(diagnostics::runtime(cli.error_format, &error)),
    }
}

fn run(command: &Command) -> anyhow::Result<()> {
    match command {
        Command::Monomerize { .. } => monomerize(command)?,
        Command::Cat { .. } => {
            concatenate(command)?;
        }
        Command::Decat { .. } => {
            deconcatenate(command)?;
        }
        Command::Canonicalize { .. } => {
            canonicalize(command)?;
        }
        Command::Uniq { .. } => {
            uniq(command)?;
        }
        Command::Rotate { .. } => rotate(command)?,
        Command::Orfs { .. } => orfs(command)?,
        Command::Schema { command } => discovery::schema(command.as_deref())?,
        Command::Completions { shell } => discovery::completions(*shell)?,
    }
    Ok(())
}
