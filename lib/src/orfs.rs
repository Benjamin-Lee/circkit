use std::{collections::HashSet, vec};

use aho_corasick::AhoCorasick;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Orf {
    /// The index of the start codon's first nucleotide. Zero-indexed.
    pub start: usize,
    /// The index of the stop codon's first nucleotide. Might be `None` if there is no stop codon. Zero-indexed.
    pub stop: Option<usize>,
    /// How many times did this ORF wrap around the origin
    pub wraps: usize,
    /// The length of the ORF in nucleotides, including the start and stop codons.
    pub length: usize,
}

// this function converts an Orf into a string with the ORF sequence
impl Orf {
    pub fn seq_with_opts(&self, seq: &[u8], include_stop: bool) -> String {
        let length = self.length - if include_stop { 0 } else { 3 };
        let mut nucleotides = Vec::with_capacity(length);
        self.write_sequence(seq, length, &mut nucleotides).unwrap();
        String::from_utf8(nucleotides).unwrap()
    }
    pub fn seq(&self, seq: &[u8]) -> String {
        self.seq_with_opts(seq, true)
    }

    /// Write circular sequence chunks directly, without materializing an ORF string.
    pub fn write_seq_with_opts<W: std::io::Write + ?Sized>(
        &self,
        seq: &[u8],
        include_stop: bool,
        writer: &mut W,
    ) -> std::io::Result<()> {
        self.write_sequence(seq, self.length - if include_stop { 0 } else { 3 }, writer)
    }

    fn write_sequence<W: std::io::Write + ?Sized>(
        &self,
        seq: &[u8],
        mut remaining: usize,
        writer: &mut W,
    ) -> std::io::Result<()> {
        if seq.is_empty() {
            return Ok(());
        }
        let mut offset = self.start % seq.len();
        while remaining > 0 {
            let count = remaining.min(seq.len() - offset);
            writer.write_all(&seq[offset..offset + count])?;
            remaining -= count;
            offset = 0;
        }
        Ok(())
    }
}

pub fn find_orfs(seq: &str) -> Vec<Orf> {
    let (starts, stops) = CodonMatcher::default().indices(seq.as_bytes());
    find_orfs_with_indices(seq.len(), starts, stops)
}

/// Reusable exact codon lookup. The common A/C/G/T alphabet uses a 64-entry table.
pub struct CodonMatcher {
    table: [u8; 64],
    fallback: Option<CodonSets>,
}

struct CodonSets {
    starts: Vec<Vec<u8>>,
    stops: Vec<Vec<u8>>,
}

impl Default for CodonMatcher {
    fn default() -> Self {
        let mut table = [0; 64];
        table[14] = 1; // ATG
        table[48] = 2; // TAA
        table[50] = 2; // TAG
        table[56] = 2; // TGA
        Self {
            table,
            fallback: None,
        }
    }
}

impl CodonMatcher {
    pub fn new(starts: &[&str], stops: &[&str]) -> Self {
        let mut matcher = Self {
            table: [0; 64],
            fallback: None,
        };
        for (codons, kind) in [(starts, 1), (stops, 2)] {
            for codon in codons {
                if let Some(index) = Self::code(codon.as_bytes()) {
                    matcher.table[index] |= kind;
                } else {
                    matcher.fallback = Some(CodonSets {
                        starts: starts.iter().map(|s| s.as_bytes().to_vec()).collect(),
                        stops: stops.iter().map(|s| s.as_bytes().to_vec()).collect(),
                    });
                    return matcher;
                }
            }
        }
        matcher
    }

    #[inline]
    fn code(codon: &[u8]) -> Option<usize> {
        fn base(b: u8) -> Option<usize> {
            match b {
                b'A' => Some(0),
                b'C' => Some(1),
                b'G' => Some(2),
                b'T' => Some(3),
                _ => None,
            }
        }
        if codon.len() != 3 {
            return None;
        }
        Some((base(codon[0])? << 4) | (base(codon[1])? << 2) | base(codon[2])?)
    }

    #[inline]
    fn kind(&self, codon: &[u8]) -> u8 {
        if let Some(codons) = &self.fallback {
            u8::from(codons.starts.iter().any(|s| s == codon))
                | (u8::from(codons.stops.iter().any(|s| s == codon)) << 1)
        } else {
            Self::code(codon)
                .map(|index| self.table[index])
                .unwrap_or(0)
        }
    }

    pub fn indices(&self, seq: &[u8]) -> (Vec<Vec<usize>>, Vec<Vec<usize>>) {
        let mut starts = vec![Vec::new(), Vec::new(), Vec::new()];
        let mut stops = vec![Vec::new(), Vec::new(), Vec::new()];
        if seq.len() < 3 {
            return (starts, stops);
        }
        for (index, codon) in seq.windows(3).enumerate() {
            let kind = self.kind(codon);
            if kind & 1 != 0 {
                starts[index % 3].push(index);
            } else if kind & 2 != 0 {
                stops[index % 3].push(index);
            }
        }
        let n = seq.len();
        for (index, codon) in [
            (n - 2, [seq[n - 2], seq[n - 1], seq[0]]),
            (n - 1, [seq[n - 1], seq[0], seq[1]]),
        ] {
            let kind = self.kind(&codon);
            if kind & 1 != 0 {
                starts[index % 3].push(index);
            }
            if kind & 2 != 0 {
                stops[index % 3].push(index);
            }
        }
        (starts, stops)
    }
}

/// A helper function to add the last two codons to a computed codon index
pub fn add_last_codons(seq: &str, codons: &[&str], codon_indices_by_frame: &mut [Vec<usize>]) {
    // Handle the last two codons wrapping around
    let penultimate_codon = format!("{}{}", &seq[seq.len() - 2..], &seq[..1]);
    debug_assert!(penultimate_codon.len() == 3);
    if codons.contains(&penultimate_codon.as_str()) {
        codon_indices_by_frame[(seq.len() - 2) % 3].push(seq.len() - 2);
    }

    let ultimate_codon = format!("{}{}", &seq[seq.len() - 1..], &seq[..2]);
    debug_assert!(ultimate_codon.len() == 3);
    if codons.contains(&ultimate_codon.as_str()) {
        codon_indices_by_frame[(seq.len() - 1) % 3].push(seq.len() - 1);
    }
}

/// Computes the indices of start and stop codons in a sequence using a simple sliding window
pub fn start_stop_codon_indices_by_frame_naive(
    seq: &str,
    start_codons: &[&str],
    stop_codons: &[&str],
) -> (Vec<Vec<usize>>, Vec<Vec<usize>>) {
    let mut start_codon_indices_by_frame = vec![Vec::new(), Vec::new(), Vec::new()];
    let mut stop_codon_indices_by_frame = vec![Vec::new(), Vec::new(), Vec::new()];

    for i in 0..seq.len() - 2 {
        let codon = &seq[i..i + 3];
        if start_codons.contains(&codon) {
            start_codon_indices_by_frame[i % 3].push(i);
        } else if stop_codons.contains(&codon) {
            stop_codon_indices_by_frame[i % 3].push(i);
        }
    }

    add_last_codons(seq, start_codons, &mut start_codon_indices_by_frame);
    add_last_codons(seq, stop_codons, &mut stop_codon_indices_by_frame);

    (start_codon_indices_by_frame, stop_codon_indices_by_frame)
}

/// Compute the indices of all start and stop codons in the sequence by frame using a Rust iterator
pub fn start_stop_codon_indices_by_frame_iter(
    seq: &str,
    start_codons: &[&str],
    stop_codons: &[&str],
) -> (Vec<Vec<usize>>, Vec<Vec<usize>>) {
    let mut start_codon_indices_by_frame = vec![Vec::new(), Vec::new(), Vec::new()];
    let mut stop_codon_indices_by_frame = vec![Vec::new(), Vec::new(), Vec::new()];

    for (i, codon) in seq
        .as_bytes()
        .iter()
        .cycle()
        .take(seq.len() + 2)
        .copied()
        .collect::<Vec<_>>()
        .windows(3)
        .enumerate()
    {
        if start_codons.contains(&std::str::from_utf8(codon).unwrap()) {
            start_codon_indices_by_frame[i % 3].push(i);
        } else if stop_codons.contains(&std::str::from_utf8(codon).unwrap()) {
            stop_codon_indices_by_frame[i % 3].push(i);
        }
    }

    (start_codon_indices_by_frame, stop_codon_indices_by_frame)
}

/// Compute the indices of all start and stop codons in the sequence by frame using the Aho-Corasick algorithm
pub fn start_stop_codon_indices_by_frame_aho_corasick(
    seq: &str,
    start_codons: &[&str],
    stop_codons: &[&str],
    aho_corasick: &AhoCorasick,
) -> (Vec<Vec<usize>>, Vec<Vec<usize>>) {
    let mut start_codon_indices_by_frame = vec![Vec::new(), Vec::new(), Vec::new()];
    let mut stop_codon_indices_by_frame = vec![Vec::new(), Vec::new(), Vec::new()];

    for mat in aho_corasick.find_overlapping_iter(seq) {
        let i = mat.start();
        if mat.pattern().as_usize() < start_codons.len() {
            start_codon_indices_by_frame[i % 3].push(i);
        } else {
            stop_codon_indices_by_frame[i % 3].push(i);
        }
    }

    // Handle the last two codons wrapping around
    add_last_codons(seq, start_codons, &mut start_codon_indices_by_frame);
    add_last_codons(seq, stop_codons, &mut stop_codon_indices_by_frame);

    (start_codon_indices_by_frame, stop_codon_indices_by_frame)
}

/// The actual ORF finding logic.
///
/// This code only works in one polarity, so it will need to be called if you want to find all ORFs in a given sequence.
/// However, it does find all ORFs in the three reading frames.
pub fn find_orfs_with_indices(
    seq_len: usize,
    start_codon_indices_by_frame: Vec<Vec<usize>>,
    stop_codon_indices_by_frame: Vec<Vec<usize>>,
) -> Vec<Orf> {
    // Find the longest ORF for each start codon
    let mut orfs = Vec::with_capacity(start_codon_indices_by_frame.iter().map(Vec::len).sum());
    let sorted: Vec<bool> = stop_codon_indices_by_frame
        .iter()
        .map(|stops| stops.windows(2).all(|pair| pair[0] <= pair[1]))
        .collect();
    let mut cursors = [0; 3];
    let mut previous = [None; 3];
    let mut next_stop = |frame: usize, start: usize, strict: bool| {
        let stops = &stop_codon_indices_by_frame[frame];
        if !sorted[frame] {
            return stops
                .iter()
                .copied()
                .find(|&stop| if strict { stop > start } else { stop >= start });
        }
        if previous[frame].is_some_and(|last| start < last) {
            cursors[frame] =
                stops.partition_point(|&stop| if strict { stop <= start } else { stop < start });
        } else {
            while cursors[frame] < stops.len()
                && (if strict {
                    stops[cursors[frame]] <= start
                } else {
                    stops[cursors[frame]] < start
                })
            {
                cursors[frame] += 1;
            }
        }
        previous[frame] = Some(start);
        stops.get(cursors[frame]).copied()
    };

    for &start_codon_index in start_codon_indices_by_frame.iter().flatten() {
        let current_frame = start_codon_index % 3;
        // let mut orf_seq = String::new();
        let mut orf_length = 0; // in nucleotides. We do this instead of storing the sequence to avoid allocating a new string for each ORF

        // Find the next stop codon in the current frame
        // The normal case has the sequence being a multiple of 3, so any wrap around is in the same frame
        // Keep compatibility with Rust versions before usize::is_multiple_of.
        #[allow(clippy::manual_is_multiple_of)]
        if seq_len % 3 == 0 {
            // Find greater than or equal to the start codon index
            let stop_codon_index =
                next_stop(current_frame, start_codon_index, true).or_else(|| {
                    let stops = &stop_codon_indices_by_frame[current_frame];
                    if sorted[current_frame] {
                        stops.first().copied().filter(|&i| i < start_codon_index)
                    } else {
                        stops.iter().copied().find(|&i| i < start_codon_index)
                    }
                });

            orfs.push(Orf {
                start: start_codon_index,
                stop: stop_codon_index,
                wraps: if stop_codon_index < Some(start_codon_index)
                    || (stop_codon_index.is_some() && seq_len - stop_codon_index.unwrap() < 3)
                {
                    1
                } else {
                    0
                },
                length: match stop_codon_index {
                    Some(stop) => match stop >= start_codon_index {
                        // Case 1: The stop codon is after the start codon
                        true => stop - start_codon_index + 3,
                        // Case 2: The stop codon is before the start codon, so the ORF wraps around
                        false => stop + seq_len - start_codon_index + 3,
                    },
                    None => seq_len, // No stop codon, so the ORF is the entire sequence
                },
            });
            continue;
        }

        // The sequence is not a multiple of 3, so we need to check the other frames
        // Find the next stop codon in the current frame
        let stop_codon_index = next_stop(current_frame, start_codon_index, false);

        // Easy case: There is a stop codon in the current frame later in the sequence
        if let Some(stop) = stop_codon_index {
            orfs.push(Orf {
                start: start_codon_index,
                stop: stop_codon_index,
                wraps: if seq_len - stop >= 3 { 0 } else { 1 },
                length: stop - start_codon_index + 3,
            });
            continue;
        } else {
            orf_length += seq_len - start_codon_index;
        }

        // Hard case: There is no stop codon in the current frame later in the sequence
        // We need to check the next frame (first round)
        // TODO: check performance of short form: (current_frame + (3 - seq.len() % 3)) % 3;
        let current_frame = if seq_len % 3 == 2 {
            (current_frame + 1) % 3
        } else {
            (current_frame + 2) % 3
        };

        // Find the next stop codon in the new current frame
        // Check if there is a stop codon in the next frame
        if let Some(stop) = stop_codon_indices_by_frame[current_frame].first().copied() {
            orfs.push(Orf {
                start: start_codon_index,
                stop: Some(stop),
                wraps: if seq_len - stop >= 3 { 1 } else { 2 },
                length: orf_length + stop + 3,
            });
            continue;
        } else {
            orf_length += seq_len;
        }

        // There is no stop codon in the next frame, so we need to check the next frame (second round)
        let current_frame = if seq_len % 3 == 2 {
            (current_frame + 1) % 3
        } else {
            (current_frame + 2) % 3
        };

        // Check if there is a stop codon in the next frame
        if let Some(stop) = stop_codon_indices_by_frame[current_frame].first().copied() {
            orfs.push(Orf {
                start: start_codon_index,
                stop: Some(stop),
                wraps: if seq_len - stop >= 3 { 2 } else { 3 },
                length: orf_length + stop + 3,
            });
            continue;
        } else {
            orf_length += seq_len;
        }

        // Finally, we're back to the original frame, so we need to check up to the start codon
        let current_frame = if seq_len % 3 == 2 {
            (current_frame + 1) % 3
        } else {
            (current_frame + 2) % 3
        };

        // Check if there is a stop codon in the next frame
        if let Some(stop) = stop_codon_indices_by_frame[current_frame].first().copied() {
            orfs.push(Orf {
                start: start_codon_index,
                stop: Some(stop),
                wraps: 3,
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
                wraps: 3,
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
        if seen_stop_codons.insert(orf.stop) {
            longest_orfs.push(*orf);
        }
    }
    longest_orfs
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn regular_linear_orf() {
        let seq: &str = "AAAATGCCCCCCCCCTAA";
        //                  ^^^123123123^^^
        let orfs = find_orfs(seq);
        assert_eq!(
            orfs,
            vec![Orf {
                start: 3,
                stop: Some(15),
                wraps: 0,
                length: 15
            }]
        );
    }

    #[test]
    fn wrap_around_once_mod_0() {
        let seq: &str = "GCATAAGCAATG";
        //                        ^^^
        //               123^^^
        let orfs = find_orfs(seq);
        assert_eq!(
            orfs,
            vec![Orf {
                start: 9,
                stop: Some(3),
                wraps: 1,
                length: 9
            }]
        );
    }

    #[test]
    fn wrap_around_once_mod_0_partial_stop() {
        let seq = "AATGAAAAAAAAATA";
        //                ^^^123123123^^
        //               ^
        let orfs = find_orfs(seq);
        assert_eq!(
            orfs,
            vec![Orf {
                start: 1,
                stop: Some(13),
                wraps: 1,
                length: 15
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
                wraps: 1,
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
                wraps: 1,
                length: 9
            }]
        );
    }

    #[test]
    fn wrap_around_once_mod_2_partial_stop() {
        let seq = "ATGAAAAAAAAATA";
        //               ^^^123123123^^
        //               ^
        let orfs = find_orfs(seq);
        assert_eq!(
            orfs,
            vec![Orf {
                start: 0,
                stop: Some(12),
                wraps: 1,
                length: 15
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
                wraps: 2,
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
                wraps: 2,
                length: 30
            },]
        );
    }

    #[test]
    fn wrap_around_twice_mod_2_partial_stop() {
        let seq: &str = "ACATACGCATG";
        //                       ^^^
        //               123123123^^
        //               ^
        let orfs = find_orfs(seq);
        assert_eq!(
            orfs,
            vec![Orf {
                start: 8,
                stop: Some(9),
                wraps: 2,
                length: 15
            }]
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
        assert_eq!(orf.wraps, 3);
    }

    #[test]
    fn wrap_around_three_times_mod_2_partial_stop() {
        let seq: &str = "AATGCAAAAAAATA";
        //                ^^^1231231231
        //               23123123123123
        //               123123123123^^
        //               ^

        let orfs = find_orfs(seq);
        assert_eq!(
            orfs,
            vec![Orf {
                start: 1,
                stop: Some(12),
                wraps: 3,
                length: 42
            },]
        );
    }
    #[test]
    fn partial_wrap_end_mod_0() {
        let seq: &str = "AAATGGCATACT";
        //                 ^^^123123^
        //               ^^

        let orfs = find_orfs(seq);

        assert_eq!(
            orfs,
            vec![Orf {
                start: 2,
                stop: Some(11),
                wraps: 1,
                length: 12
            }]
        )
    }

    #[test]
    fn partial_wrap_end() {
        let seq: &str = "ATGGCATA";
        //               ^^^123^^
        //               ^

        let orfs = find_orfs(seq);

        assert_eq!(
            orfs,
            vec![Orf {
                start: 0,
                stop: Some(6),
                wraps: 1,
                length: 9
            }]
        )
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
                wraps: 0,
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
                let bio_orfs: Vec<bio::seq_analysis::orf::Orf> = finder.find_all(seq.as_bytes()).collect::<Vec<_>>();
                let circkit_orfs: Vec<Orf> = find_orfs(&seq);

                // for each bio orf, make sure there is a circkit orf that is the same
                for bio_orf in &bio_orfs {
                    let circkit_orf = circkit_orfs.iter().find(|circkit_orf| {
                        circkit_orf.start == bio_orf.start && circkit_orf.stop.unwrap() + 3 == bio_orf.end
                    });

                    // We must find all ORFs that Rust-Bio does
                    prop_assert!(circkit_orf.is_some(), "Circkit ORF: {:?} not found in bio orf: {:?}", circkit_orf, bio_orf);

                    // We know that Rust-Bio doesn't work on wrapped ORFs, so any one that matches should not be wrapped
                    prop_assert_eq!(circkit_orf.unwrap().wraps, 0, "Circkit ORF: {:?} ({:?}) should not be a wrapped ORF", circkit_orf, circkit_orf.unwrap().seq(seq.as_bytes()));

                    if circkit_orf.is_some() {
                        prop_assert_eq!((circkit_orf.unwrap().start % 3) as i8, bio_orf.offset, "Circkit {:?} not in same offset as Bio {:?}", circkit_orf, bio_orf);
                    }
                }

                // ensure that each ORF that was found by circkit but not Rust-Bio is a wrapped ORF
                for circkit_orf in circkit_orfs {
                    let bio_orf = &bio_orfs.iter().find(|bio_orf| {
                        circkit_orf.start == bio_orf.start && circkit_orf.stop.unwrap_or_default() + 3 == bio_orf.end
                    });

                    if circkit_orf.wraps == 0 {
                        prop_assert!(bio_orf.is_some())
                    }

                    if bio_orf.is_none() {
                        prop_assert_ne!(circkit_orf.wraps, 0, "{:?} should be a wrapped ORF", circkit_orf);
                    }
                }
            }

            #[test]
            fn test_bio_orfs_repeated(seq in "[ATGC]{3,300}") {
                // Quadruple the sequence and make sure that the ORFs are the same
                let dup_seq = format!("{}{}{}{}", seq, seq, seq, seq);
                let finder = Finder::new(vec![b"ATG"], vec![b"TAA", b"TAG", b"TGA"], 0);
                let bio_orfs: Vec<bio::seq_analysis::orf::Orf> = finder.find_all(dup_seq.as_bytes()).collect::<Vec<_>>();
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
                    prop_assert!(orfs.contains(orf), "Longest ORF: {:?} is not in all ORFs: {:?}", orf, orfs);
                }
                prop_assert!(longest_orfs.len() <= orfs.len(), "Longest ORFs: {:?} is not a subset of all ORFs: {:?}", longest_orfs, orfs);
            }

            #[test]
            fn indexing_is_identical(seq in "[ATGC]{3,300}"){
                let start_codons = ["ATG"];
                let stop_codons = ["TAA", "TAG", "TGA"];
                let ac = aho_corasick::AhoCorasick::new(["ATG", "TAA", "TAG", "TGA"]).unwrap();
                prop_assert_eq!(start_stop_codon_indices_by_frame_naive(&seq, &start_codons, &stop_codons), start_stop_codon_indices_by_frame_iter(&seq, &start_codons, &stop_codons));
                prop_assert_eq!(start_stop_codon_indices_by_frame_naive(&seq, &start_codons, &stop_codons), start_stop_codon_indices_by_frame_aho_corasick(&seq, &start_codons, &stop_codons, &ac));
            }
        }
    }
}
