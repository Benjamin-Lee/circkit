use bio::alignment::distance::simd::*;
use bio::pattern_matching::shift_and;
use log::warn;

// pub fn check_overlap(seq: &[u8]) -> Option<&[u8]> {
//     // slice last ten bases of the record
//     let seed = &seq[seq.len() - 10..];
//     let remaining = &seq[..seq.len() - 10];

//     let shiftand = shift_and::ShiftAnd::new(seed);
//     let occ = shiftand.find_all(remaining).next();

//     // if there isn't a match to the seed, return None
//     if occ.is_none() {
//         return None;
//     }

//     let occ = occ.unwrap();

//     let potential_overlap = &seq[0..occ + seed.len()];

//     let dist = hamming(
//         potential_overlap,
//         &seq[seq.len() - potential_overlap.len()..],
//     );

//     // decide what to return
//     if dist < 10 {
//         check_overlap(&seq[occ..])
//     } else {
//         None
//     }
// }

pub fn monomerize(seq: &[u8], seed_length: usize, overlap_max_mismatch: u64) -> &[u8] {
    // if the sequence is less than 10 bases long, return None
    if seq.len() < seed_length {
        warn!("Sequence is less than seed length");
        return seq;
    }

    // slice last n bases of the record
    let seed = &seq[seq.len() - seed_length..];

    // we start out assuming the seed isn't even present anywhere else
    let mut last_match: Option<usize> = None;

    // create the matcher
    let shift_and = shift_and::ShiftAnd::new(seed);
    let searcher = shift_and.find_all(&seq[..seq.len() - seed_length]);

    // loop over the searcher (noting that the final n bases are not included)
    for occ in searcher {
        // slice the potential overlap from the sequence
        let potential_overlap = &seq[last_match.unwrap_or(0)..occ + seed.len()];

        // extend backward from the end of the sequence to the start of the potential overlap
        let seed_extension = &seq[seq.len() - potential_overlap.len()..];

        // compare the potential overlap to the seed
        let dist = hamming(potential_overlap, seed_extension);

        // decide whether the overlap is good enough to be a monomer
        if dist <= overlap_max_mismatch {
            last_match = Some(occ);
        }
    }

    match last_match {
        Some(last_match) => &seq[last_match + seed.len()..],
        None => seq,
    }
}

#[cfg(test)]
mod tests {
    use crate::monomerize::monomerize;
    use pretty_assertions::assert_eq;

    #[test]
    fn no_overlap_at_all() {
        assert_eq!(
            monomerize(b"TTTTTTTTTTTTAAAAAAAAAA", 10, 1),
            b"TTTTTTTTTTTTAAAAAAAAAA"
        );
    }
    #[test]
    fn identical_repeat() {
        assert_eq!(monomerize(b"ATGCATGC", 2, 1), b"ATGC");
    }

    #[test]
    fn big_with_mismatch_overlap() {
        let x = b"TTTTTGGTTTTTAAAAAAAAAATTTTTTTTTTTTAAAAAAAAAA";
        assert_eq!(monomerize(x, 10, 3), b"TTTTTTTTTTTTAAAAAAAAAA");
        assert_eq!(monomerize(x, 10, 2), b"TTTTTTTTTTTTAAAAAAAAAA");
        assert_eq!(monomerize(x, 10, 1), x);
        assert_eq!(monomerize(x, 10, 0), x);
    }
    #[test]
    fn repeat_at_beginning() {
        assert_eq!(monomerize(b"AAAAATTTTTAAAAATTTTT", 10, 0), b"AAAAATTTTT");
    }

    #[test]
    fn multimer() {
        assert_eq!(
            monomerize(b"AAAAATTTTTAAAAATTTTTAAAAATTTTT", 10, 0),
            b"AAAAATTTTT"
        );
    }

    #[test]
    fn dimer() {
        let input = b"TAAAAAAAAAAAAATAAAAAAAAAAAAA";
        let output = b"TAAAAAAAAAAAAA";

        for seed_len in 6..=10 {
            assert_eq!(monomerize(input, seed_len, 0), output);
            println!("\n");
        }
    }
}
