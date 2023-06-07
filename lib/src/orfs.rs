use std::{collections::HashSet, vec};

use aho_corasick::AhoCorasick;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Orf {
    /// The index of the start codon's first nucleotide.
    pub start: usize,
    /// The index of the stop codon's first nucleotide. Might be `None` if there is no stop codon.
    pub stop: Option<usize>,
    pub frame_shift: usize,
    /// The length of the ORF in nucleotides, including the start and stop codons.
    pub length: usize,
}

// this function converts an Orf into a string with the ORF sequence
impl Orf {
    pub fn seq(&self, seq: &[u8]) -> String {
        // use a cyclical iterator to get the nucleotides, starting at the start codon
        let nucleotides = seq
            .iter()
            .cycle()
            .skip(self.start)
            .take(self.length)
            .copied()
            .collect::<Vec<_>>();
        String::from_utf8(nucleotides).unwrap()
    }
}

pub fn find_orfs(seq: &str) -> Vec<Orf> {
    // Step 1: Find all stop and start codons by frame
    let start_codons = ["ATG"];
    let stop_codons = ["TAA", "TAG", "TGA"];

    let mut start_codon_indices_by_frame = vec![Vec::new(), Vec::new(), Vec::new()];
    let mut stop_codon_indices_by_frame = vec![Vec::new(), Vec::new(), Vec::new()];

    // read from env whether to use aho-corasick or not
    if std::env::var("AHO_CORASICK").is_err() {
        for i in 0..seq.len() - 2 {
            let codon = &seq[i..i + 3];
            if start_codons.contains(&codon) {
                start_codon_indices_by_frame[i % 3].push(i);
            } else if stop_codons.contains(&codon) {
                stop_codon_indices_by_frame[i % 3].push(i);
            }
        }
    } else {
        let patterns = start_codons
            .iter()
            .chain(stop_codons.iter())
            .map(|s| s.as_bytes())
            .collect::<Vec<_>>();
        let ac = AhoCorasick::new(patterns).unwrap();

        for mat in ac.find_overlapping_iter(seq) {
            let i = mat.start();
            if mat.pattern().as_usize() < start_codons.len() {
                start_codon_indices_by_frame[i % 3].push(i);
            } else {
                stop_codon_indices_by_frame[i % 3].push(i);
            }
        }
    }

    // Handle the last two codons wrapping around
    let penultimate_codon = format!("{}{}", &seq[seq.len() - 2..], &seq[..1]);
    debug_assert!(penultimate_codon.len() == 3);
    if start_codons.contains(&penultimate_codon.as_str()) {
        start_codon_indices_by_frame[(seq.len() - 2) % 3].push(seq.len() - 2);
    } else if stop_codons.contains(&penultimate_codon.as_str()) {
        stop_codon_indices_by_frame[(seq.len() - 2) % 3].push(seq.len() - 2);
    }

    let ultimate_codon = format!("{}{}", &seq[seq.len() - 1..], &seq[..2]);
    debug_assert!(ultimate_codon.len() == 3);
    if start_codons.contains(&ultimate_codon.as_str()) {
        start_codon_indices_by_frame[(seq.len() - 1) % 3].push(seq.len() - 1);
    } else if stop_codons.contains(&ultimate_codon.as_str()) {
        stop_codon_indices_by_frame[(seq.len() - 1) % 3].push(seq.len() - 1);
    }

    // Step 3: Find the longest ORF for each start codon
    let mut orfs = Vec::new();

    for &start_codon_index in start_codon_indices_by_frame.iter().flatten() {
        let current_frame = start_codon_index % 3;
        // let mut orf_seq = String::new();
        let mut orf_length = 0; // in nucleotides. We do this instead of storing the sequence to avoid allocating a new string for each ORF

        // Find the next stop codon in the current frame
        // The normal case has the sequence being a multiple of 3, so any wrap around is in the same frame
        if seq.len() % 3 == 0 {
            // Find greater than or equal to the start codon index
            let stop_codon_index = stop_codon_indices_by_frame[current_frame]
                .iter()
                .copied()
                .find(|&i| i > start_codon_index)
                .or_else(|| {
                    stop_codon_indices_by_frame[current_frame]
                        .iter()
                        .copied()
                        .find(|&i| i < start_codon_index)
                });

            orfs.push(Orf {
                start: start_codon_index,
                stop: stop_codon_index,
                frame_shift: 0,
                length: match stop_codon_index {
                    Some(stop) => match stop >= start_codon_index {
                        // Case 1: The stop codon is after the start codon
                        true => stop - start_codon_index + 3,
                        // Case 2: The stop codon is before the start codon, so the ORF wraps around
                        false => stop + seq.len() - start_codon_index + 3,
                    },
                    None => seq.len(), // No stop codon, so the ORF is the entire sequence
                },
            });
            continue;
        }

        // The sequence is not a multiple of 3, so we need to check the other frames
        // Find the next stop codon in the current frame
        let stop_codon_index = stop_codon_indices_by_frame[current_frame]
            .iter()
            .copied()
            .find(|&i| i >= start_codon_index);

        // Easy case: There is a stop codon in the current frame later in the sequence
        if let Some(stop) = stop_codon_index {
            orfs.push(Orf {
                start: start_codon_index,
                stop: Some(stop),
                frame_shift: 0,
                length: stop - start_codon_index + 3,
            });
            continue;
        } else {
            orf_length += seq.len() - start_codon_index;
        }

        // Hard case: There is no stop codon in the current frame later in the sequence
        // We need to check the next frame (first round)
        // TODO: check performance of short form: (current_frame + (3 - seq.len() % 3)) % 3;
        let current_frame = if seq.len() % 3 == 2 {
            (current_frame + 1) % 3
        } else {
            (current_frame + 2) % 3
        };

        // Find the next stop codon in the new current frame
        // Check if there is a stop codon in the next frame
        if let Some(stop) = stop_codon_indices_by_frame[current_frame].get(0).copied() {
            orfs.push(Orf {
                start: start_codon_index,
                stop: Some(stop),
                frame_shift: 1,
                length: orf_length + stop + 3,
            });
            continue;
        } else {
            orf_length += seq.len();
        }

        // There is no stop codon in the next frame, so we need to check the next frame (second round)
        let current_frame = if seq.len() % 3 == 2 {
            (current_frame + 1) % 3
        } else {
            (current_frame + 2) % 3
        };

        // Check if there is a stop codon in the next frame
        if let Some(stop) = stop_codon_indices_by_frame[current_frame].get(0).copied() {
            orfs.push(Orf {
                start: start_codon_index,
                stop: Some(stop),
                frame_shift: 2,
                length: orf_length + stop + 3,
            });
            continue;
        } else {
            orf_length += seq.len();
        }

        // Finally, we're back to the original frame, so we need to check up to the start codon
        let current_frame = if seq.len() % 3 == 2 {
            (current_frame + 1) % 3
        } else {
            (current_frame + 2) % 3
        };

        // Check if there is a stop codon in the next frame
        if let Some(stop) = stop_codon_indices_by_frame[current_frame].get(0).copied() {
            orfs.push(Orf {
                start: start_codon_index,
                stop: Some(stop),
                frame_shift: 3,
                length: orf_length + stop + 3,
            });
            continue;
        } else {
            orf_length += start_codon_index;
            // There is no stop codon in the next frame and we have already checked all frames
            // This means that the ORF wraps around the sequence infinitely
            orfs.push(Orf {
                start: start_codon_index,
                stop: None,
                frame_shift: 3,
                length: orf_length,
            });
        }
    }

    orfs
}

/// For each stop codon, keep only the longest ORF
pub fn longest_orfs(orfs: &mut Vec<Orf>) -> Vec<Orf> {
    // For each stop codon, keep only the longest ORF
    orfs.sort_by_key(|orf| orf.length); // TODO: make this an unstable sort for performance (if it makes a difference)
    orfs.reverse();
    let mut longest_orfs = Vec::new();
    let mut seen_stop_codons = HashSet::new(); // TODO: check performance of HashSet vs. Vec vs alternative hasher
    for orf in orfs {
        if !seen_stop_codons.contains(&orf.stop) {
            seen_stop_codons.insert(orf.stop);
            longest_orfs.push(*orf);
        }
    }
    longest_orfs
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn wrap_around_once_mod_0() {
        let seq = "GCATAAGCAATG";
        //                        ^^^
        //               123^^^
        let orfs = find_orfs(seq);
        assert_eq!(
            orfs,
            vec![Orf {
                start: 9,
                stop: Some(3),
                frame_shift: 0,
                length: 9
            }]
        );
    }

    #[test]
    fn wrap_around_once_mod_1() {
        let seq = "GCATAAGATG";
        //                      ^^^
        //               123^^^
        let orfs = find_orfs(seq);
        assert_eq!(
            orfs,
            vec![Orf {
                start: 7,
                stop: Some(3),
                frame_shift: 1,
                length: 9
            }]
        );
    }

    #[test]
    fn wrap_around_once_mod_2() {
        let seq = "GCATAAGCATG";
        //                       ^^^
        //               123^^^
        let orfs = find_orfs(seq);
        assert_eq!(
            orfs,
            vec![Orf {
                start: 8,
                stop: Some(3),
                frame_shift: 1,
                length: 9
            }]
        );
    }

    #[test]
    fn wrap_around_twice_mod_1() {
        let seq = "ATGAAAAAAAAAA";
        //               ^^^1231231231
        //               2312312312312
        //               3^^^

        let orfs = find_orfs(seq);
        assert_eq!(
            orfs,
            vec![Orf {
                start: 0,
                stop: Some(1),
                frame_shift: 2,
                length: 30
            },]
        );
    }

    #[test]
    fn wrap_around_twice_mod_2() {
        let seq = "AATGCATAAAA";
        //                ^^^1231231
        //               23123123123
        //               123123^^^

        let orfs = find_orfs(seq);
        assert_eq!(
            orfs,
            vec![Orf {
                start: 1,
                stop: Some(6),
                frame_shift: 2,
                length: 30
            },]
        );
    }

    #[test]
    fn wrap_around_three_times() {
        // This sequence only a stop codon in the same frame as the start codon but requires wrapping around the sequence three times
        let seq: &str = "GGTCGGAGAATTGGGTCAGTTTCGGGCTTAAAAACTCTGACTTGTCATGCTCGTGGCGTCCCTACCG";
        //                                                             ^^^123123123123123123
        //               1231231231231231231231231231231231231231231231231231231231231231231
        //               2312312312312312312312312312312312312312312312312312312312312312312
        //               3123123123123123123123123123^^^
        let orfs = find_orfs(seq);
        assert_eq!(orfs.len(), 1);
        let orf = orfs[0];
        assert_eq!(orf.start, 46);
        assert_eq!(orf.stop.unwrap(), 28);
        assert_eq!(orf.seq(seq.as_bytes()), "ATGCTCGTGGCGTCCCTACCGGGTCGGAGAATTGGGTCAGTTTCGGGCTTAAAAACTCTGACTTGTCATGCTCGTGGCGTCCCTACCGGGTCGGAGAATTGGGTCAGTTTCGGGCTTAAAAACTCTGACTTGTCATGCTCGTGGCGTCCCTACCGGGTCGGAGAATTGGGTCAGTTTCGGGCTTAA");
        assert_eq!(orf.frame_shift, 3);
    }

    #[test]
    fn longest_orf_small() {
        let seq = "ATGATGTAG";
        //               ^^^123^^^

        let mut orfs = find_orfs(seq);
        let longest_orfs = longest_orfs(&mut orfs);

        assert_eq!(
            longest_orfs,
            vec![Orf {
                start: 0,
                stop: Some(6),
                frame_shift: 0,
                length: 9
            }]
        );
    }

    /// Use proptest to ensure that the orfs are the same as the ones generated by Rust-Bio
    mod fuzz {
        use super::*;
        use bio::seq_analysis::orf::Finder;
        use proptest::prelude::*;
        proptest! {
            #[test]
            fn test_bio_orfs(seq in "[ATGC]{3,300}") {
                let finder = Finder::new(vec![b"ATG"], vec![b"TAA", b"TAG", b"TGA"], 0);
                let bio_orfs: Vec<bio::seq_analysis::orf::Orf> = finder.find_all(seq.as_bytes()).into_iter().collect::<Vec<_>>();
                let circkit_orfs: Vec<Orf> = find_orfs(&seq);

                // for each bio orf, make sure there is a circkit orf that is the same
                for bio_orf in &bio_orfs {
                    let circkit_orf = circkit_orfs.iter().find(|circkit_orf| {
                        circkit_orf.start == bio_orf.start && circkit_orf.stop.unwrap() + 3 == bio_orf.end
                    });

                    prop_assert!(circkit_orf.is_some(), "Circkit ORF: {:?} not found in bio orf: {:?}", circkit_orf, bio_orf);
                    if circkit_orf.is_some() {
                        prop_assert_eq!((circkit_orf.unwrap().start % 3) as i8, bio_orf.offset, "Circkit {:?} not in same offset as Bio {:?}", circkit_orf, bio_orf);
                    }
                }

                // ensure that each ORF that was found by circkit but not bio is a wrapped ORF
                for circkit_orf in circkit_orfs {
                    let bio_orf = &bio_orfs.iter().find(|bio_orf| {
                        circkit_orf.start == bio_orf.start && circkit_orf.stop.unwrap_or_default() + 3 == bio_orf.end
                    });
                    if bio_orf.is_none() {
                        // one reason bio might miss an ORF is if it wrapped around a divisible sequence
                        if seq.len() % 3 == 0 {
                            if let Some(stop) = circkit_orf.stop {
                                // the seq.len() - stop <= 2 is there because we might have a wrap the is less than an entire codon's length
                                prop_assert!(circkit_orf.start > stop || seq.len() - stop <= 2, "{:?} should have a start position greater than the stop position", circkit_orf);
                            }
                        }

                        // otherwise, we had a frame shift wrap around
                        if seq.len() % 3 != 0 {
                            prop_assert!(circkit_orf.frame_shift > 0 || seq.len() - circkit_orf.stop.unwrap_or_default() <= 2, "{:?} should have a frame shift", circkit_orf);
                        }
                    }
                }
            }

            #[test]
            fn test_bio_orfs_repeated(seq in "[ATGC]{3,300}") {
                // Quadruple the sequence and make sure that the ORFs are the same
                let dup_seq = format!("{}{}{}{}", seq, seq, seq, seq);
                let finder = Finder::new(vec![b"ATG"], vec![b"TAA", b"TAG", b"TGA"], 0);
                let bio_orfs: Vec<bio::seq_analysis::orf::Orf> = finder.find_all(dup_seq.as_bytes()).into_iter().collect::<Vec<_>>();
                let circkit_orfs: Vec<Orf> = find_orfs(&seq);

                // for each bio orf, make sure there is a circkit orf that is the same
                for bio_orf in bio_orfs {
                    let bio_orf_seq = dup_seq[bio_orf.start..bio_orf.end].to_owned();
                    let circkit_orf = circkit_orfs.iter().find(|circkit_orf| {
                        let circkit_orf_seq = circkit_orf.seq(seq.as_bytes());
                        circkit_orf_seq == bio_orf_seq
                    });

                    prop_assert!(circkit_orf.is_some(), "Circkit ORF not found for bio {:?}", bio_orf);
                }
            }

            #[test]
            fn orf_never_has_start_codon_as_stop_codon(seq in "[ATGC]{3,300}") {
                let orfs = find_orfs(&seq);
                for orf in orfs {
                    prop_assert_ne!(Some(orf.start), orf.stop, "ORF: {:?} has start codon as stop codon", orf);
                }
            }

            #[test]
            fn orf_is_divisible_by_three(seq in "[ATGC]{3,300}") {
                let orfs = find_orfs(&seq);
                for orf in orfs {
                    prop_assert_eq!((orf.length) % 3, 0, "ORF: {:?} is not divisible by 3", orf);
                }
            }

            #[test]
            fn longest_orfs_are_subset_of_all_orfs(seq in "[ATGC]{3,300}") {
                let mut orfs = find_orfs(&seq);
                let longest_orfs = longest_orfs(&mut orfs);
                for orf in &longest_orfs {
                    prop_assert!(orfs.contains(&orf), "Longest ORF: {:?} is not in all ORFs: {:?}", orf, orfs);
                }
                prop_assert!(longest_orfs.len() <= orfs.len(), "Longest ORFs: {:?} is not a subset of all ORFs: {:?}", longest_orfs, orfs);
            }

        }
    }
}
