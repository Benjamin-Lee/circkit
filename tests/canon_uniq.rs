use assert_cmd::prelude::*; // Add methods on commands
use assert_fs::prelude::*; // Add methods on paths
use predicates::prelude::*; // Used for writing assertions
use rstest::rstest; // Parameterized tests
use std::process::Command; // Run programs
mod common;

#[test]
fn file_doesnt_exist() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("circkit")?;

    cmd.arg("canonicalize").arg("test/file/doesnt/exist");
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("No such file or directory"));
    Ok(())
}

#[test]
fn simple_fasta_file_to_stdout() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("circkit")?;

    let file = assert_fs::NamedTempFile::new("simple.fasta")?;
    file.write_str(concat!(">seq1", "\n", "ATGCA"))?;

    cmd.arg("canonicalize").arg(file.path());
    cmd.assert()
        .success()
        .stdout(predicate::str::contains(">seq1\nAATGC"));
    Ok(())
}

#[rstest]
// base case
#[case("canonicalize", "simple")]
#[case("uniq", "simple")]
// a base case but where there are multiple sequences in the file
#[case("canonicalize", "multiple_sequences")]
#[case("uniq", "multiple_sequences")]
// the same as above, but with the sequences split across multiple lines
#[case("canonicalize", "multiple_sequences_split_lines")]
#[case("uniq", "multiple_sequences_split_lines")]
// an example where the input file has an RNA sequence
#[case("canonicalize", "rna_input")]
#[case("uniq", "rna_input")]
// an example where the same sequence is repeated multiple times but in different rotations and polarities
#[case("uniq", "repeated")]
fn fasta_files(
    #[case] command: &str,
    #[case] directory: &str,
    #[values("1", "2", "4")] threads: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("circkit")?;

    let file = std::path::Path::new("tests/examples")
        .join(directory)
        .join("in.fasta");

    let output = assert_fs::NamedTempFile::new("out.fasta")?;

    cmd.arg(command)
        .arg(file.as_path())
        .arg("--threads")
        .arg(threads)
        .arg("-o")
        .arg(output.path());

    // when we're testing canonicalize and uniq, using this flag should usually produce the same output
    if command == "uniq" {
        cmd.arg("--canonicalize");
    }

    // it should succeed and not print anything to stdout or stderr
    cmd.assert()
        .success()
        .stdout(predicate::str::is_empty())
        .stderr(predicate::str::is_empty());

    let known_good = std::path::Path::new("tests/examples")
        .join(directory)
        .join("out.fasta");

    assert!(common::sequences_are_identical(
        known_good.to_str().unwrap(),
        output.path().to_str().unwrap()
    ));

    Ok(())
}
