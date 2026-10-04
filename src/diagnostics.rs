use crate::commands::ErrorFormat;
use std::{ffi::OsString, io::Write};

#[derive(Debug)]
pub struct InvalidArguments(pub String);

impl std::fmt::Display for InvalidArguments {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter.write_str(&self.0)
    }
}

impl std::error::Error for InvalidArguments {}

pub fn argument_error(message: impl Into<String>) -> anyhow::Error {
    InvalidArguments(message.into()).into()
}

/// Recover the explicitly requested diagnostic format even when argument parsing
/// fails. A '--' delimiter ends flag parsing, including this scan.
pub fn requested_format(arguments: &[OsString]) -> ErrorFormat {
    let mut format = ErrorFormat::Text;
    let mut arguments = arguments.iter().skip(1);
    while let Some(argument) = arguments.next() {
        if argument == "--" {
            break;
        }
        let value = if argument == "--error-format" {
            arguments.next().and_then(|value| value.to_str())
        } else {
            argument
                .to_str()
                .and_then(|value| value.strip_prefix("--error-format="))
        };
        match value {
            Some("json") => format = ErrorFormat::Json,
            Some("text") => format = ErrorFormat::Text,
            _ => {}
        }
    }
    format
}

pub fn emit(format: ErrorFormat, code: &str, message: &str, exit_code: u8) {
    let mut stderr = std::io::stderr().lock();
    match format {
        ErrorFormat::Text => {
            let _ = writeln!(stderr, "error: {message}");
        }
        ErrorFormat::Json => {
            let error = serde_json::json!({
                "schema_version": 1,
                "error": { "code": code, "message": message, "exit_code": exit_code }
            });
            let _ = serde_json::to_writer(&mut stderr, &error);
            let _ = stderr.write_all(b"\n");
        }
    }
}

pub fn runtime(format: ErrorFormat, error: &anyhow::Error) -> u8 {
    let usage = error.chain().any(|source| source.is::<InvalidArguments>());
    let exit_code = if usage { 2 } else { 1 };
    let code = if usage {
        "invalid_arguments"
    } else if error.chain().any(|source| {
        source
            .downcast_ref::<seq_io::fasta::Error>()
            .is_some_and(|error| !matches!(error, seq_io::fasta::Error::Io(_)))
    }) {
        "invalid_fasta"
    } else if error
        .chain()
        .any(|source| source.is::<std::str::Utf8Error>())
    {
        "invalid_header"
    } else if error.chain().any(|source| source.is::<std::io::Error>()) {
        "io_error"
    } else {
        "processing_error"
    };
    emit(format, code, &format!("{error:#}"), exit_code);
    exit_code
}
