use assert_cmd::prelude::*; // Add methods on commands
use std::process::Command; // Run programs
mod common;

#[test]
fn cat() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("circkit")?;
    cmd.arg("cat");
    common::check_fasta("cat", &mut cmd)?;
    Ok(())
}

#[test]
fn decat() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("circkit")?;
    cmd.arg("decat");
    common::check_fasta("decat", &mut cmd)?;
    Ok(())
}
