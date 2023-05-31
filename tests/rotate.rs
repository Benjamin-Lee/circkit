use assert_cmd::prelude::*; // Add methods on commands
use std::process::Command; // Run programs
mod common;
use rstest::rstest;

#[rstest]
#[case("rotate_5", 5)]
#[case("rotate_minus_5", -5)]
fn rotate_bases(#[case] directory: &str, #[case] bases: i64) -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("circkit")?;
    cmd.arg("rotate");
    cmd.arg("--bases").arg(bases.to_string());
    common::check_fasta(directory, &mut cmd)?;
    Ok(())
}

#[rstest]
#[case("rotate_0.25", 0.25)]
#[case("rotate_0.5", 0.5)]
fn rotate_percent(#[case] directory: &str, #[case] percent: f64) -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("circkit")?;
    cmd.arg("rotate");
    cmd.arg("--percent").arg(percent.to_string());
    common::check_fasta(directory, &mut cmd)?;
    Ok(())
}
