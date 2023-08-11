use assert_cmd::prelude::*; // Add methods on commands
use predicates::prelude::*; // Used for writing assertions
use rstest::rstest; // Parameterized tests
use std::process::Command; // Run programs
mod common;

#[rstest]
#[case("nim_cated_Soil_microbial_communities_from_permafrost_in_Bonanza_Creek__Alaska")]
fn fasta_files(
    #[case] directory: &str,
    #[values("1", "2", "4")] threads: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("circkit")?;

    let file = std::path::Path::new("tests/examples")
        .join(directory)
        .join("in.fasta.xz");

    let output = assert_fs::NamedTempFile::new("out.fasta")?;

    cmd.arg("monomerize")
        .arg(file.as_path())
        .arg("--threads")
        .arg(threads)
        .arg("--min-identity")
        .arg("0.95")
        .arg("-o")
        .arg(output.path());

    // it should succeed and not print anything to stdout or stderr
    cmd.assert()
        .success()
        .stdout(predicate::str::is_empty())
        .stderr(predicate::str::is_empty());

    let known_good = std::path::Path::new("tests/examples")
        .join(directory)
        .join("out.fasta");
    assert!(common::sequences_are_subset(
        output.path().to_str().unwrap(),
        known_good.to_str().unwrap()
    ));

    Ok(())
}

#[rstest]
fn min_overlap(#[values("1", "2", "4")] threads: &str) -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("circkit")?;

    cmd.arg("monomerize")
        .arg("--threads")
        .arg(threads)
        .arg("--min-overlap")
        .arg("81");

    common::check_fasta("min_overlap", &mut cmd)?;

    Ok(())
}

#[rstest]
#[case("0.51")]
#[case("1.0")]
#[case("1.5")]
fn min_overlap_percent(
    #[case] cutoff: &str,
    #[values("1", "2", "4")] threads: &str,
) -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("circkit")?;

    cmd.arg("monomerize")
        .arg("--threads")
        .arg(threads)
        .arg("--min-overlap-percent")
        .arg(cutoff);

    common::check_fasta(format!("min_overlap_percent_{}", cutoff).as_str(), &mut cmd)?;

    Ok(())
}
