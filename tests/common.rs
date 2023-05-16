use bio::io::fasta;
use std::collections::HashMap;

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

    println!("{:?}", seqs1);
    println!("{:?}", seqs2);

    seqs1 == seqs2
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
    }
}
