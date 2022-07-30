// use bio::alignment::distance::simd::*;
// use bio::io::fasta;
// use bio::pattern_matching::shift_and;
// use std::io;
// use std::str::from_utf8;

// fn check_overlap(seq: &[u8]) -> Option<&[u8]> {
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

// pub fn monomerize() -> bool {
//     true
// }

pub mod concatenate;
#[macro_use]
extern crate log;
