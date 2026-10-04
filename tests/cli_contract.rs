use assert_cmd::Command;
use serde_json::Value;
use std::{collections::BTreeMap, io::Read};

fn cli(args: &[&str], input: &[u8]) -> std::process::Output {
    Command::cargo_bin("circkit")
        .unwrap()
        .args(args)
        .write_stdin(input.to_vec())
        .output()
        .unwrap()
}

fn json_error(output: &std::process::Output, exit: i32, code: &str) -> Value {
    assert_eq!(output.status.code(), Some(exit), "{output:?}");
    let error: Value = serde_json::from_slice(&output.stderr).unwrap();
    assert_eq!(error["schema_version"], 1);
    assert_eq!(error["error"]["code"], code);
    assert_eq!(error["error"]["exit_code"], exit);
    assert!(!error["error"]["message"].as_str().unwrap().is_empty());
    error
}

#[test]
fn malformed_fasta_fails_in_every_command_and_execution_mode() {
    for operation in [
        "cat",
        "decat",
        "rotate",
        "canonicalize",
        "uniq",
        "monomerize",
        "orfs",
    ] {
        for mode in ["serial", "pipeline"] {
            let mut args = vec![operation, "--error-format", "json"];
            if operation == "rotate" {
                args.extend(["--bases", "1"]);
            } else if !["cat", "decat"].contains(&operation) {
                args.extend(["--threads", "1", "--execution", mode]);
            }
            json_error(&cli(&args, b"this is not FASTA\n"), 1, "invalid_fasta");
        }
    }
}

#[test]
fn invalid_arguments_have_machine_readable_errors_before_io() {
    let directory = assert_fs::TempDir::new().unwrap();
    let output = directory.path().join("preserve.fasta");
    for arguments in [
        vec!["rotate"],
        vec!["rotate", "--bases", "0"],
        vec!["rotate", "--percent", "nan"],
        vec!["rotate", "--bases", "1", "--percent", "0.5"],
        vec!["monomerize", "--min-identity", "1.1"],
        vec!["monomerize", "--min-overlap-percent", "-0.5"],
        vec!["monomerize", "--min-length", "10", "--max-length", "5"],
        vec!["orfs", "--min-ratio", "inf"],
        vec!["orfs", "--min-wraps", "3", "--max-wraps", "1"],
        vec!["orfs", "--max-wraps", "4"],
        vec!["orfs", "--start-codons", "AT"],
        vec!["orfs", "--stop-codons", "TAA,"],
        vec!["orfs", "--table-format", "jsonl"],
        vec!["canonicalize", "--threads", "2", "--execution", "serial"],
    ] {
        for first in [false, true] {
            std::fs::write(&output, b"keep me").unwrap();
            let mut command = Command::cargo_bin("circkit").unwrap();
            if first {
                command.args(["--error-format=json"]);
            }
            command.args(&arguments).arg("-o").arg(&output);
            if !first {
                command.args(["--error-format", "json"]);
            }
            json_error(&command.output().unwrap(), 2, "invalid_arguments");
            assert_eq!(std::fs::read(&output).unwrap(), b"keep me");
        }
    }
}

#[test]
fn short_empty_wrapped_streams_and_signed_rotations_work() {
    for operation in [
        "cat",
        "decat",
        "rotate",
        "canon",
        "uniq",
        "monomerize",
        "orfs",
    ] {
        let mut args = vec![operation, "-", "-o", "-"];
        if operation == "rotate" {
            args.extend(["--bases", "1"]);
        } else if !["cat", "decat"].contains(&operation) {
            args.extend(["--threads", "1"]);
        }
        if operation == "monomerize" {
            args.push("--keep-all");
        }
        for input in [b"".as_slice(), b">x\n", b">x\n\n"] {
            let result = cli(&args, input);
            assert!(result.status.success(), "{operation}: {result:?}");
            assert!(result.stderr.is_empty());
        }
    }
    let input = b">one description\r\nAC\r\nGT\r\nA\r\n>empty\r\n\r\n";
    for (bases, expected) in [
        ("1", "AACGT"),
        ("-1", "CGTAA"),
        ("6", "AACGT"),
        ("-9223372036854775808", "TAACG"),
    ] {
        let result = cli(&["rotate", "-", "-o", "-", "--bases", bases], input);
        assert!(result.status.success(), "{result:?}");
        assert_eq!(
            result.stdout,
            format!(">one description\n{expected}\n>empty\n\n").as_bytes()
        );
    }
    for (percent, expected) in [("0.5", "TAACG"), ("-0.5", "TAACG")] {
        let result = cli(&["rotate", "--percent", percent], input);
        assert!(result.status.success(), "{result:?}");
        assert_eq!(
            result.stdout,
            format!(">one description\n{expected}\n>empty\n\n").as_bytes()
        );
    }
    json_error(
        &cli(
            &["rotate", "--percent", "1e300", "--error-format", "json"],
            input,
        ),
        2,
        "invalid_arguments",
    );
    assert_eq!(
        cli(&["cat"], input).stdout,
        b">one description\nACGTAACGTA\n>empty\n\n"
    );
    assert_eq!(
        cli(&["decat"], input).stdout,
        b">one description\nAC\n>empty\n\n"
    );
}

#[test]
fn aliases_global_flags_and_codons_preserve_compatibility() {
    let input = b">a\nATGAAATAAATGAAATAA\n>b\nATGAAATAAATGAAATAA\n";
    for alias in ["--norm", "--canon", "--canonicalize"] {
        let result = cli(
            &["uniq", alias, "--threads", "1", "--batch-size", "1", "-q"],
            input,
        );
        assert!(result.status.success(), "{alias}: {result:?}");
        assert!(result.stderr.is_empty());
    }
    let result = cli(
        &[
            "monomerize",
            "--min-overlap-ratio",
            "1.0",
            "--keep-all",
            "--threads",
            "1",
        ],
        input,
    );
    assert!(result.status.success());
    let expected = cli(
        &[
            "orfs",
            "--threads",
            "1",
            "--min-length",
            "0",
            "--start-codons",
            "ATG,GTG",
        ],
        input,
    );
    let actual = cli(
        &[
            "orfs",
            "--threads",
            "1",
            "--min-length",
            "0",
            "--start-codons",
            "atg, gtg",
            "--stop-codons",
            "taa, tag, tga",
        ],
        input,
    );
    assert!(actual.status.success(), "{actual:?}");
    assert_eq!(expected.stdout, actual.stdout);
}

#[test]
fn schema_describes_actual_choices_aliases_and_rotation_requirements() {
    let all = cli(&["schema"], b"");
    assert!(all.status.success(), "{all:?}");
    assert!(all.stderr.is_empty());
    let catalog: Value = serde_json::from_slice(&all.stdout).unwrap();
    assert_eq!(catalog["schema_version"], 1);
    assert_eq!(catalog["commands"].as_array().unwrap().len(), 9);
    for (name, canonical) in [
        ("canon", "canonicalize"),
        ("rotate", "rotate"),
        ("describe", "schema"),
    ] {
        let result = cli(&["describe", name], b"");
        assert!(result.status.success(), "{result:?}");
        let schema: Value = serde_json::from_slice(&result.stdout).unwrap();
        assert_eq!(schema["commands"][0]["name"], canonical);
        if name == "rotate" {
            let groups = schema["commands"][0]["groups"].as_array().unwrap();
            let rotation = groups
                .iter()
                .find(|group| group["name"] == "rotation")
                .unwrap();
            assert_eq!(rotation["required"], true);
            assert_eq!(rotation["multiple"], false);
        }
    }
    let orfs = catalog["commands"]
        .as_array()
        .unwrap()
        .iter()
        .find(|cmd| cmd["name"] == "orfs")
        .unwrap();
    let strand = orfs["arguments"]
        .as_array()
        .unwrap()
        .iter()
        .find(|arg| arg["name"] == "strand")
        .unwrap();
    assert_eq!(
        strand["choices"],
        serde_json::json!(["forward", "reverse", "both"])
    );
    json_error(
        &cli(
            &["schema", "no-such-command", "--error-format", "json"],
            b"",
        ),
        2,
        "invalid_arguments",
    );
}

#[test]
fn completions_and_help_are_clean_and_available_without_input() {
    for shell in ["bash", "fish", "zsh", "powershell", "elvish"] {
        let output = cli(&["completions", shell], b"");
        assert!(output.status.success(), "{output:?}");
        assert!(output.stderr.is_empty());
        let script = String::from_utf8(output.stdout).unwrap();
        assert!(script.contains("circkit"));
        assert!(script.contains("table-format"));
    }
    for args in [vec!["--help"], vec!["--version"], vec!["orfs", "--help"]] {
        let output = cli(&args, b"");
        assert!(output.status.success());
        assert!(output.stderr.is_empty());
        assert!(!output.stdout.is_empty());
    }
}

#[test]
fn metadata_formats_preserve_fields_and_values_and_finish_compression() {
    let directory = assert_fs::TempDir::new().unwrap();
    let input =
        b">one description\nATGAAATAAATGAAATAA\n>duplicate description\nATGAAATAAATGAAATAA\n";
    for operation in ["monomerize", "uniq", "orfs"] {
        let schema_output = cli(&["schema", operation], b"");
        let schema: Value = serde_json::from_slice(&schema_output.stdout).unwrap();
        let expected_fields = schema["commands"][0]["metadata_fields"].as_array().unwrap();
        let mut reference = Vec::new();
        for extension in ["csv", "tsv", "jsonl", "jsonl.gz", "tsv.gz"] {
            let table = directory.path().join(format!("table.{extension}"));
            let mut command = Command::cargo_bin("circkit").unwrap();
            command
                .args([operation, "--threads", "1", "--table"])
                .arg(&table)
                .write_stdin(input.to_vec());
            if operation == "monomerize" {
                command.arg("--keep-all");
            }
            if operation == "orfs" {
                command.args(["--min-length", "0"]);
            }
            command.assert().success();
            let bytes = if extension.ends_with(".gz") {
                let mut bytes = Vec::new();
                flate2::read::MultiGzDecoder::new(std::fs::File::open(&table).unwrap())
                    .read_to_end(&mut bytes)
                    .unwrap();
                bytes
            } else {
                std::fs::read(&table).unwrap()
            };
            let rows: Vec<BTreeMap<String, String>> = if extension.starts_with("jsonl") {
                String::from_utf8(bytes)
                    .unwrap()
                    .lines()
                    .map(|line| {
                        let row: BTreeMap<String, Value> = serde_json::from_str(line).unwrap();
                        row.into_iter()
                            .map(|(key, value)| {
                                (
                                    key,
                                    match value {
                                        Value::String(value) => value,
                                        Value::Null => String::new(),
                                        value => value.to_string(),
                                    },
                                )
                            })
                            .collect()
                    })
                    .collect()
            } else {
                csv::ReaderBuilder::new()
                    .delimiter(if extension.starts_with("tsv") {
                        b'\t'
                    } else {
                        b','
                    })
                    .from_reader(bytes.as_slice())
                    .deserialize()
                    .map(Result::unwrap)
                    .collect()
            };
            assert!(!rows.is_empty(), "{operation}");
            let mut fields: Vec<_> = expected_fields
                .iter()
                .map(|field| field.as_str().unwrap().to_owned())
                .collect();
            fields.sort();
            assert_eq!(rows[0].keys().cloned().collect::<Vec<_>>(), fields);
            if extension == "csv" {
                reference = rows;
            } else {
                assert_eq!(reference, rows, "{operation} {extension}");
            }
        }
    }
    let output = directory.path().join("output.fasta");
    let result = Command::cargo_bin("circkit")
        .unwrap()
        .args([
            "orfs",
            "--threads",
            "1",
            "--min-length",
            "0",
            "--strand",
            "forward",
            "--no-stop-required",
            "--table",
            "-",
            "--table-format",
            "jsonl",
            "-o",
        ])
        .arg(output)
        .write_stdin(">partial\nATGAAA\n")
        .output()
        .unwrap();
    assert!(result.status.success(), "{result:?}");
    let row: Value = serde_json::from_slice(&result.stdout).unwrap();
    assert_eq!(row["stop"], Value::Null);
}

#[test]
fn conflicting_paths_are_rejected_without_truncating_source_files() {
    let directory = assert_fs::TempDir::new().unwrap();
    let input = directory.path().join("input.fasta");
    let data = b">original\nACGTACGT\n";
    std::fs::write(&input, data).unwrap();
    let hardlink = directory.path().join("hardlink.fasta");
    std::fs::hard_link(&input, &hardlink).unwrap();
    let mut aliases = vec![input.clone(), hardlink];
    #[cfg(unix)]
    {
        let symlink = directory.path().join("symlink.fasta");
        std::os::unix::fs::symlink(&input, &symlink).unwrap();
        aliases.push(symlink);
    }
    for alias in aliases {
        for target in ["-o", "--table"] {
            let result = Command::cargo_bin("circkit")
                .unwrap()
                .args(["uniq", "--error-format", "json"])
                .arg(&input)
                .arg(target)
                .arg(&alias)
                .output()
                .unwrap();
            json_error(&result, 2, "invalid_arguments");
            assert_eq!(std::fs::read(&input).unwrap(), data);
        }
    }
    let result = cli(&["uniq", "--table", "-", "--error-format", "json"], data);
    json_error(&result, 2, "invalid_arguments");
    let output = directory.path().join("new.fasta");
    let result = Command::cargo_bin("circkit")
        .unwrap()
        .args(["uniq", "--error-format", "json"])
        .arg(&input)
        .arg("-o")
        .arg(&output)
        .arg("--table")
        .arg(directory.path().join("./new.fasta"))
        .output()
        .unwrap();
    json_error(&result, 2, "invalid_arguments");
    assert!(!output.exists());
}

#[cfg(target_os = "linux")]
#[test]
fn file_and_metadata_write_errors_propagate_through_the_pipeline() {
    for (operation, flags) in [
        ("cat", vec![]),
        ("decat", vec![]),
        ("rotate", vec!["--bases", "1"]),
        (
            "canonicalize",
            vec!["--threads", "1", "--execution", "pipeline"],
        ),
        ("uniq", vec!["--threads", "1", "--execution", "serial"]),
        ("monomerize", vec!["--keep-all", "--threads", "1"]),
        ("orfs", vec!["--min-length", "0", "--threads", "1"]),
    ] {
        let mut args = vec![operation, "--error-format", "json", "-o", "/dev/full"];
        args.extend(flags);
        let output = cli(&args, b">one\nATGAAATAAATGAAATAA\n");
        let error = json_error(&output, 1, "io_error");
        assert!(error["error"]["message"]
            .as_str()
            .unwrap()
            .contains("/dev/full"));
    }
    let output = cli(
        &[
            "uniq",
            "--threads",
            "1",
            "--table",
            "/dev/full",
            "--error-format",
            "json",
        ],
        b">one\nACGTACGT\n>two\nACGTACGT\n",
    );
    json_error(&output, 1, "io_error");
}

#[cfg(unix)]
#[test]
fn closed_stdout_is_quiet_success_in_serial_and_pipeline_commands() {
    use std::process::{Command as Process, Stdio};
    let directory = assert_fs::TempDir::new().unwrap();
    let input = directory.path().join("large.fasta");
    std::fs::write(&input, b">one\nATGAAATAAATGAAATAA\n".repeat(10000)).unwrap();
    for operation in [
        "cat",
        "rotate",
        "canonicalize",
        "orfs",
        "schema",
        "completions",
    ] {
        let mut process = Process::new(assert_cmd::cargo::cargo_bin("circkit"));
        process
            .arg(operation)
            .stdout(Stdio::piped())
            .stderr(Stdio::piped());
        if operation == "completions" {
            process.arg("bash");
        } else if operation != "schema" {
            process.arg(&input);
            if operation == "rotate" {
                process.args(["--bases", "1"]);
            } else if operation != "cat" {
                process.args(["--threads", "1", "--execution", "pipeline"]);
            }
            if operation == "orfs" {
                process.args(["--min-length", "0"]);
            }
        }
        let mut child = process.spawn().unwrap();
        drop(child.stdout.take());
        let result = child.wait_with_output().unwrap();
        assert!(result.status.success(), "{operation}: {result:?}");
        assert!(result.stderr.is_empty(), "{operation}: {result:?}");
    }
}

#[cfg(unix)]
#[test]
fn byte_headers_and_native_paths_work_without_metadata() {
    let directory = assert_fs::TempDir::new().unwrap();
    // APFS enforces UTF-8 filenames. Linux also exercises arbitrary path bytes;
    // non-UTF-8 FASTA headers are valid on both filesystems.
    #[cfg(target_os = "macos")]
    let input = directory.path().join("input séquences.fasta");
    #[cfg(not(target_os = "macos"))]
    let input = {
        use std::os::unix::ffi::OsStringExt;
        directory
            .path()
            .join(std::ffi::OsString::from_vec(b"input \xff.fasta".to_vec()))
    };
    let output = directory.path().join("output with spaces.fasta");
    std::fs::write(&input, b">id\xff\nACGT\n").unwrap();
    Command::cargo_bin("circkit")
        .unwrap()
        .arg("uniq")
        .arg(&input)
        .args(["--threads", "1", "-o"])
        .arg(&output)
        .assert()
        .success();
    assert_eq!(std::fs::read(output).unwrap(), b">id\xff\nACGT\n");
    let result = Command::cargo_bin("circkit")
        .unwrap()
        .arg("uniq")
        .arg(input)
        .args(["--threads", "1", "--error-format", "json", "--table"])
        .arg(directory.path().join("metadata.jsonl"))
        .output()
        .unwrap();
    json_error(&result, 1, "invalid_header");
}

#[test]
fn truncated_compressed_input_is_an_io_error_and_multimember_gzip_is_read() {
    use std::io::Write;
    let mut compressed = Vec::new();
    for data in [b">one\nACGT\n".as_slice(), b">two\nTGCA\n"] {
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        encoder.write_all(data).unwrap();
        compressed.extend(encoder.finish().unwrap());
    }
    let result = cli(&["cat", "-"], &compressed);
    assert!(result.status.success(), "{result:?}");
    assert_eq!(result.stdout, b">one\nACGTACGT\n>two\nTGCATGCA\n");
    compressed.truncate(compressed.len() - 4);
    for args in [
        vec!["cat"],
        vec!["canonicalize", "--threads", "1", "--execution", "pipeline"],
    ] {
        let mut args = args;
        args.extend(["--error-format", "json"]);
        json_error(&cli(&args, &compressed), 1, "io_error");
    }
}

#[cfg(unix)]
#[test]
fn closed_metadata_stdout_preserves_pipe_handling_for_csv_and_jsonl() {
    use std::process::{Command as Process, Stdio};
    let directory = assert_fs::TempDir::new().unwrap();
    let input = directory.path().join("large.fasta");
    std::fs::write(&input, b">one\nATGAAATAAATGAAATAA\n".repeat(10000)).unwrap();
    for format in ["csv", "jsonl"] {
        let mut child = Process::new(assert_cmd::cargo::cargo_bin("circkit"))
            .arg("monomerize")
            .arg(&input)
            .args([
                "--threads",
                "1",
                "--table",
                "-",
                "--table-format",
                format,
                "-o",
            ])
            .arg(directory.path().join("out.fasta"))
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .unwrap();
        drop(child.stdout.take());
        let result = child.wait_with_output().unwrap();
        assert!(result.status.success(), "{format}: {result:?}");
        assert!(result.stderr.is_empty(), "{format}: {result:?}");
    }
}

#[test]
fn orf_strand_selection_restricts_fasta_and_metadata_in_both_execution_modes() {
    let directory = assert_fs::TempDir::new().unwrap();
    let input = b">forward\nATGAAATAA\n>reverse\nTTATTTCAT\n";
    for mode in ["serial", "pipeline"] {
        for (strand, headers) in [
            ("forward", vec!["forward_ORF0"]),
            ("reverse", vec!["reverse_RC_ORF0"]),
            ("both", vec!["forward_ORF0", "reverse_RC_ORF0"]),
        ] {
            let table = directory.path().join("orfs.jsonl");
            let output = Command::cargo_bin("circkit")
                .unwrap()
                .args([
                    "orfs",
                    "--threads",
                    "1",
                    "--execution",
                    mode,
                    "--min-length",
                    "0",
                    "--max-wraps",
                    "0",
                    "--include-stop",
                    "--strand",
                    strand,
                    "--table",
                ])
                .arg(&table)
                .write_stdin(input.to_vec())
                .output()
                .unwrap();
            assert!(output.status.success(), "{strand}: {output:?}");
            let fasta = String::from_utf8(output.stdout).unwrap();
            let actual: Vec<_> = fasta
                .lines()
                .filter_map(|line| line.strip_prefix('>'))
                .collect();
            assert_eq!(actual, headers, "{mode} {strand}");
            for sequence in fasta.lines().filter(|line| !line.starts_with('>')) {
                assert_eq!(sequence, "ATGAAATAA");
            }
            let metadata = std::fs::read_to_string(table).unwrap();
            let actual: Vec<String> = metadata
                .lines()
                .map(|line| {
                    let row: Value = serde_json::from_str(line).unwrap();
                    row["orf_id"].as_str().unwrap().to_owned()
                })
                .collect();
            assert_eq!(actual, headers, "{mode} {strand} metadata");
        }
    }
}
