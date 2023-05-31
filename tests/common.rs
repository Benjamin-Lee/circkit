use assert_cmd::prelude::*;
use bio::io::fasta;
use predicates::prelude::*;
use std::collections::HashMap;
use std::process::Command; // Run programs // Used for writing assertions // Add methods on commands

pub fn check_fasta(directory: &str, cmd: &mut Command) -> anyhow::Result<()> {
    let file = std::path::Path::new("tests/examples")
        .join(directory)
        .join("in.fasta");

    let output = assert_fs::NamedTempFile::new("out.fasta")?;

    cmd.arg(file.as_path()).arg("-o").arg(output.path());

    // it should succeed and not print anything to stdout or stderr
    cmd.assert()
        .success()
        .stdout(predicate::str::is_empty())
        .stderr(predicate::str::is_empty());

    let known_good = std::path::Path::new("tests/examples")
        .join(directory)
        .join("out.fasta");
    assert!(sequences_are_identical(
        output.path().to_str().unwrap(),
        known_good.to_str().unwrap()
    ));

    Ok(())
}

pub fn sequences_are_identical(file1: &str, file2: &str) -> bool {
    let mut seqs1 = HashMap::new();
    let mut seqs2 = HashMap::new();

    let mut reader1 = fasta::Reader::from_file(file1).unwrap().records();
    let mut reader2 = fasta::Reader::from_file(file2).unwrap().records();

    while let Some(Ok(record)) = reader1.next() {
        seqs1.insert(record.id().to_string(), record.seq().to_vec());
    }

    while let Some(Ok(record)) = reader2.next() {
        seqs2.insert(record.id().to_string(), record.seq().to_vec());
    }

    // println!("{:?}", seqs1);
    // println!("{:?}", seqs2);

    // iterate over the keys in seqs1 and check if they are the same in seqs2, printing the sequences if not
    for key in seqs1.keys() {
        if !seqs2.contains_key(key) {
            println!("{} not found in {}", key, file2);
            return false;
        }
        if seqs1.get(key).unwrap() != seqs2.get(key).unwrap() {
            println!(
                "{} not the same in {}.\nFile1:\n{}\nFile2:\n{}\n",
                key,
                file2,
                String::from_utf8_lossy(seqs1.get(key).unwrap()),
                String::from_utf8_lossy(seqs2.get(key).unwrap())
            );
            return false;
        }
    }

    // iterate over the keys in seqs2 and check if they are the same in seqs1, printing the sequences if not
    for key in seqs2.keys() {
        if !seqs1.contains_key(key) {
            println!("{} not found in {}", key, file1);
            return false;
        }
        if seqs1.get(key).unwrap() != seqs2.get(key).unwrap() {
            println!(
                "{} not the same in {}.\nFile1:\n{}\nFile2:\n{}\n",
                key,
                file1,
                String::from_utf8_lossy(seqs1.get(key).unwrap()),
                String::from_utf8_lossy(seqs2.get(key).unwrap())
            );
            return false;
        }
    }

    seqs1 == seqs2
}

/// Check that the sequences in file2 are a subset of the sequences in file1
pub fn sequences_are_subset(superset: &str, subset: &str) -> bool {
    let mut superset_seqs = HashMap::new();
    let mut subset_seqs = HashMap::new();

    let mut superset_reader = fasta::Reader::from_file(superset).unwrap().records();
    let mut subset_reader = fasta::Reader::from_file(subset).unwrap().records();

    while let Some(Ok(record)) = superset_reader.next() {
        superset_seqs.insert(record.id().to_string(), record.seq().to_vec());
    }

    while let Some(Ok(record)) = subset_reader.next() {
        subset_seqs.insert(record.id().to_string(), record.seq().to_vec());
    }

    let mut return_value = true;

    // iterate over the keys in seqs2 and check if they are the same in seqs1, printing the sequences if not
    for key in subset_seqs.keys() {
        if !superset_seqs.contains_key(key) {
            println!("{} not found in {}", key, superset);
            return_value = false;
        }
    }

    // also check that the sequences in the superset are the same or substrings of the sequences in subset
    for (key, seq) in subset_seqs.iter() {
        let superseq_seq = String::from_utf8(superset_seqs.get(key).unwrap().to_owned()).unwrap();
        let subset_seq = String::from_utf8(seq.to_owned()).unwrap();
        if superseq_seq != subset_seq && !subset_seq.starts_with(std::str::from_utf8(seq).unwrap())
        {
            println!(
                "{} not the same in {}.\nSuperset:\n{}\nSubset:\n{}\n",
                key,
                superset,
                String::from_utf8_lossy(superset_seqs.get(key).unwrap()),
                String::from_utf8_lossy(seq)
            );
            return_value = false;
        }
    }

    return_value
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hashmaps_are_identical() {
        assert!(sequences_are_identical(
            "tests/examples/simple/in.fasta",
            "tests/examples/simple/out.fasta"
        ));
        assert!(sequences_are_subset(
            "tests/examples/simple/in.fasta",
            "tests/examples/simple/out.fasta"
        ));
    }
}
