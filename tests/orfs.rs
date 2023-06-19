use assert_cmd::prelude::*;
use assert_cmd::prelude::*; // Add methods on commands
use bio::io::fasta;
use predicates::prelude::*;
use predicates::prelude::*; // Used for writing assertions
use rstest::rstest;
use std::collections::HashSet;
use std::process::Command; // Run programs // Used for writing assertions // Add methods on commands // Parameterized tests

pub fn orfs_are_identical(file1: &str, file2: &str) -> bool {
    let mut seqs1 = HashSet::new();

    let mut reader1 = fasta::Reader::from_file(file1).unwrap().records();
    let mut reader2 = fasta::Reader::from_file(file2).unwrap().records();

    while let Some(Ok(record)) = reader1.next() {
        seqs1.insert(record.seq().to_vec());
    }

    while let Some(Ok(record)) = reader2.next() {
        if !seqs1.contains(record.seq()) {
            println!(
                ">{}\n{} not found in {}",
                record.desc().unwrap(),
                std::str::from_utf8(record.seq()).unwrap(),
                file1
            );
            return false;
        }
    }
    true
}

#[rstest]
fn orfipy(#[values("1", "2", "4")] threads: &str) -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("circkit")?;

    // we're reusing the input data from monomerization to compare the outputs from Orfipy
    let in_directory =
        "nim_cated_Soil_microbial_communities_from_permafrost_in_Bonanza_Creek__Alaska";
    let out_directory =
        "orfs_cated_Soil_microbial_communities_from_permafrost_in_Bonanza_Creek__Alaska";

    let file = std::path::Path::new("tests/examples")
        .join(in_directory)
        .join("in.fasta.xz");
    let known_good = std::path::Path::new("tests/examples")
        .join(out_directory)
        .join("out.fasta");

    let output = assert_fs::NamedTempFile::new("out.fasta")?;

    cmd.arg("orfs")
        .arg(file.as_path())
        .arg("--threads")
        .arg(threads.to_string())
        .arg("--start-codons")
        .arg("ATG,CTG,TTG")
        .arg("--max-wraps")
        .arg("0")
        .arg("--include-stop")
        .arg("-o")
        .arg(output.path());
    cmd.assert().success();

    let mut seqs_good = HashSet::new();
    let mut seqs_circkit = HashSet::new();

    let mut reader_good = fasta::Reader::from_file(known_good).unwrap().records();
    let mut reader_circkit = fasta::Reader::from_file(output).unwrap().records();

    while let Some(Ok(record)) = reader_good.next() {
        seqs_good.insert(record.seq().to_vec());
    }
    while let Some(Ok(record)) = reader_circkit.next() {
        seqs_circkit.insert(record.seq().to_vec());
    }

    assert_eq!(seqs_good, seqs_circkit);

    Ok(())
}
