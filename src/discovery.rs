//! Introspection from the same command tree used by argument parsing.
use crate::{
    commands::Cli,
    diagnostics::argument_error,
    io::{json_write_error, output_to_writer},
};
use clap::{ArgAction, CommandFactory};
use serde_json::{json, Value};
use std::{any::TypeId, io::Write, num::NonZeroUsize, path::PathBuf};

fn argument_type(argument: &clap::Arg) -> &'static str {
    match argument.get_action() {
        ArgAction::SetTrue | ArgAction::SetFalse => return "boolean",
        ArgAction::Count => return "integer",
        _ => {}
    }
    let kind = argument.get_value_parser().type_id();
    if kind == TypeId::of::<PathBuf>() {
        "path"
    } else if kind == TypeId::of::<u64>()
        || kind == TypeId::of::<u32>()
        || kind == TypeId::of::<usize>()
        || kind == TypeId::of::<NonZeroUsize>()
        || kind == TypeId::of::<i64>()
    {
        "integer"
    } else if kind == TypeId::of::<f64>() {
        "number"
    } else {
        "string"
    }
}

fn command_schema(command: &clap::Command) -> Value {
    let arguments: Vec<_> = command.get_arguments().filter(|argument| !argument.is_hide_set()).map(|argument| {
        let mut conflicts: Vec<_> = command.get_arg_conflicts_with(argument).into_iter()
            .filter(|other| other.get_id() != argument.get_id()).map(|other| other.get_id().to_string()).collect();
        conflicts.sort();
        conflicts.dedup();
        let range = argument.get_num_args();
        json!({
            "name": argument.get_id().as_str(),
            "long": argument.get_long(),
            "short": argument.get_short().map(|short| short.to_string()),
            "aliases": argument.get_all_aliases().unwrap_or_default(),
            "position": argument.get_index(),
            "type": argument_type(argument),
            "action": format!("{:?}", argument.get_action()),
            "required": argument.is_required_set(),
            "global": argument.is_global_set(),
            "defaults": argument.get_default_values().iter().map(|value| value.to_string_lossy()).collect::<Vec<_>>(),
            "choices": argument.get_possible_values().iter().filter(|value| !value.is_hide_set()).map(|value| value.get_name()).collect::<Vec<_>>(),
            "min_values": range.map(|range| range.min_values()),
            "max_values": range.and_then(|range| (range.max_values() != usize::MAX).then_some(range.max_values())),
            "conflicts_with": conflicts,
            "help": argument.get_long_help().or_else(|| argument.get_help()).map(|help| help.to_string()),
        })
    }).collect();
    let groups: Vec<_> = command
        .get_groups()
        .map(|group| {
            json!({
                "name": group.get_id().as_str(),
                "arguments": group.get_args().map(|argument| argument.as_str()).collect::<Vec<_>>(),
                "required": group.is_required_set(),
                "multiple": group.clone().is_multiple(),
            })
        })
        .collect();
    let mut built = command.clone();
    json!({
        "name": command.get_name(),
        "aliases": command.get_all_aliases().collect::<Vec<_>>(),
        "about": command.get_long_about().or_else(|| command.get_about()).map(|about| about.to_string()),
        "usage": built.render_usage().to_string(),
        "arguments": arguments,
        "groups": groups,
        "metadata_fields": match command.get_name() {
            "monomerize" => vec!["id", "original_length", "monomer_length"],
            "uniq" => vec!["id", "duplicate_id"],
            "orfs" => vec!["orf_id", "seq_id", "start", "stop", "length", "wraps", "ratio"],
            _ => vec![],
        },
    })
}

pub fn schema(requested: Option<&str>) -> anyhow::Result<()> {
    let mut root = Cli::command();
    root.build();
    let commands: Vec<_> = if let Some(name) = requested {
        let command = root.find_subcommand(name).ok_or_else(|| {
            argument_error(format!(
                "unknown command {name:?}; run 'circkit schema' to list commands"
            ))
        })?;
        vec![command_schema(command)]
    } else {
        root.get_subcommands()
            .filter(|command| !command.is_hide_set() && command.get_name() != "help")
            .map(command_schema)
            .collect()
    };
    let catalog = json!({
        "schema_version": 1,
        "program": "circkit",
        "version": env!("CARGO_PKG_VERSION"),
        "commands": commands,
        "streams": {
            "input": "FASTA; '-' or omitted input reads stdin; compression is detected from bytes",
            "output": "FASTA; '-' or omitted output writes stdout; file suffix selects gzip/bzip2/xz/zstd",
            "metadata": "CSV, TSV, or JSONL; '--table -' writes stdout when FASTA uses a named output file",
            "diagnostics": "stderr; '--error-format json' selects versioned JSON errors",
        },
        "exit_codes": { "0": "success, help/version, or stdout closed by a downstream consumer", "1": "input/output or processing failure", "2": "invalid arguments" },
        "metadata_formats": ["csv", "tsv", "jsonl"],
        "ordering": "Multiple workers may reorder FASTA batches and choose different duplicate representatives; use --threads 1 for stable order",
    });
    let mut writer = output_to_writer(&None)?;
    serde_json::to_writer_pretty(&mut writer, &catalog).map_err(json_write_error)?;
    writer.write_all(b"\n")?;
    writer.finish()
}

pub fn completions(shell: clap_complete::Shell) -> anyhow::Result<()> {
    let mut bytes = Vec::new();
    clap_complete::generate(shell, &mut Cli::command(), "circkit", &mut bytes);
    let mut writer = output_to_writer(&None)?;
    writer.write_all(&bytes)?;
    writer.finish()
}
