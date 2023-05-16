use assert_cmd::prelude::*; // Add methods on commands
use assert_fs::prelude::*; // Add methods on paths
use predicates::prelude::*; // Used for writing assertions
use rstest::rstest; // Parameterized tests
use std::process::Command; // Run programs
mod common;

#[rstest]
/// This test is a simple check that the different compression formats are supported for output.
/// The test data itself is rather simple
fn compressed_output(
    #[values("canonicalize", "uniq", "monomerize")] command: &str,
    #[values("gz", "bz2", "xz", "zst")] extension: &str,
    #[values("1", "2", "4")] threads: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("circkit")?;

    let file = std::path::Path::new("tests/examples")
        .join("compressed_output")
        .join("in.fasta");

    let output_dir = assert_fs::TempDir::new()?;
    let output = output_dir.child(format!("out.fasta.{}", extension));
    let output_decompressed = output_dir.child("out.fasta");

    cmd.arg(command)
        .arg(file.as_path())
        .arg("--threads")
        .arg(threads.to_string())
        .arg("-o")
        .arg(output.path());

    // when we're testing canonicalize and uniq, using this flag should usually produce the same output
    if command == "uniq" {
        cmd.arg("--canonicalize");
    }

    // when we're testing monomerize, using this flag should result a no-op for the test data
    if command == "monomerize" {
        cmd.arg("--keep-all");
    }

    // it should succeed and not print anything to stdout or stderr
    cmd.assert()
        .success()
        .stdout(predicate::str::is_empty())
        .stderr(predicate::str::is_empty());

    let mut decompress_cmd = Command::new(match extension {
        "gz" => "gzip",
        "bz2" => "bzip2",
        "xz" => "xz",
        "zst" => "zstd",
        _ => panic!("invalid extension"),
    });
    decompress_cmd.arg("-d").arg(output.path());
    decompress_cmd.assert().success();

    let known_good = std::path::Path::new("tests/examples")
        .join("compressed_output")
        .join("out.fasta");

    assert!(common::sequences_are_identical(
        known_good.to_str().unwrap(),
        output_decompressed.path().to_str().unwrap()
    ));

    Ok(())
}

#[rstest]
/// This test is a simple check that the different compression formats are supported for input.
fn compressed_input(
    #[values("canonicalize", "uniq", "monomerize")] command: &str,
    #[values("gz", "bz2", "xz", "zst")] extension: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("circkit")?;

    let file = std::path::Path::new("tests/examples")
        .join("compressed_input")
        .join(format!("in.fasta.{}", extension));

    let output = assert_fs::NamedTempFile::new("out.fasta")?;

    cmd.arg(command)
        .arg(file.as_path())
        .arg("-o")
        .arg(output.path());

    // when we're testing canonicalize and uniq, using this flag should usually produce the same output
    if command == "uniq" {
        cmd.arg("--canonicalize");
    }

    // when we're testing monomerize, using this flag should result a no-op for the test data
    if command == "monomerize" {
        cmd.arg("--keep-all");
    }

    // it should succeed and not print anything to stdout or stderr
    cmd.assert()
        .success()
        .stdout(predicate::str::is_empty())
        .stderr(predicate::str::is_empty());

    let known_good = std::path::Path::new("tests/examples")
        .join("compressed_input")
        .join("out.fasta");

    assert!(common::sequences_are_identical(
        known_good.to_str().unwrap(),
        output.path().to_str().unwrap()
    ));

    Ok(())
}
