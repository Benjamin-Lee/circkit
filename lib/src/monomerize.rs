use bio::alignment::distance::simd::*;
use bio::alphabets::dna;
use bio::pattern_matching::shift_and;
use log::{debug, warn};

#[derive(Builder, Default, Clone, Copy)]
#[builder(setter(strip_option), build_fn(validate = "Self::validate"))]
pub struct Monomerizer {
    /// The maximum number of mismatches allowed in an overlap. Conflicts with `overlap_min_identity`.
    #[builder(default)]
    pub overlap_dist: Option<u64>,
    /// The minimum percent identity within an overlap that may be considered a match. Conflicts with `overlap_dist`.
    #[builder(default)]
    pub overlap_min_identity: Option<f64>,
    /// The size of the seed to search for in the overlap.
    pub seed_len: usize,
}

impl MonomerizerBuilder {
    fn validate(&self) -> Result<(), String> {
        if self.overlap_dist.is_some() && self.overlap_min_identity.is_some() {
            // there's no support for overlap_dist and overlap_min_identity at the same time yet
            // TODO: allow users to specify both and choose the stricter/looser one
            return Err("Both overlap_dist and overlap_min_identity are set. They are mutually exclusive since they may produce conflicting filtering results.".to_string());
        }

        if let Some(seed_len) = self.seed_len {
            match seed_len {
                1..=63 => {}
                _ => {
                    return Err(format!(
                        "Seed length must be at least 1 and at most 63 but was set to {}.",
                        seed_len
                    ))
                }
            }
        }

        Ok(())
    }
}

impl Monomerizer {
    pub fn builder() -> MonomerizerBuilder {
        MonomerizerBuilder::default()
    }
    /// Compute the index of the last base of the first monomer in the sequence, if found.
    pub fn first_monomer_end_index(self, seq: &[u8]) -> Option<usize> {
        // if the sequence is shorter than the seed, give up
        let seed_len = self.seed_len;
        if seq.len() <= seed_len {
            warn!("Sequence is not longer than seed length");
            return None;
        }

        // slice last n bases of the record
        let seed = &seq[seq.len() - seed_len..];

        // create a seed matcher
        let matcher = shift_and::ShiftAnd::new(seed);

        for occ in matcher.find_all(&seq[..seq.len() - seed_len]) {
            let successor_seed = &seq[..occ + seed_len];
            let starter_seed = &seq[seq.len() - successor_seed.len()..];

            // compare the potential overlap to the seed
            let dist = hamming(starter_seed, successor_seed);

            // compute the maximum distance allowed for the overlap
            let max_dist = match self.overlap_min_identity {
                Some(identity) => {
                    successor_seed.len() as u64
                        - (successor_seed.len() as f64 * identity).floor() as u64
                }
                None => self.overlap_dist.unwrap_or(0),
            };

            debug!(
                "occ: {}, dist: {}, max_dist: {}\nstarter:\t1\t{}\t{}\nsuccessor:\t{}\t{}\t{}\n\n",
                occ,
                dist,
                max_dist,
                std::str::from_utf8(starter_seed).unwrap(),
                occ + starter_seed.len(),
                seq.len() - starter_seed.len(),
                std::str::from_utf8(successor_seed).unwrap(),
                seq.len(),
            );

            // decide whether the overlap is good enough to be a monomer
            if dist <= max_dist {
                return Some(seq.len() - starter_seed.len());
            }
        }
        None
    }
    pub fn last_monomer_end_index(self, seq: &[u8]) -> Option<usize> {
        let mut monomerized = self.first_monomer_end_index(seq);
        debug!("monomerized index (first pass): {:?}\n", monomerized);
        while let Some(monomer_index) = monomerized {
            debug!(
                "new monomer: {}\n",
                std::str::from_utf8(&seq[..monomer_index]).unwrap()
            );
            let new_monomer = self.first_monomer_end_index(&seq[..monomer_index]);
            debug!("new monomer index: {:?}\n", new_monomer);
            if new_monomer.is_none() {
                debug!("no new monomer found");
                break;
            }
            monomerized = new_monomer;
        }
        debug!(
            "Final monomer: {:?}\n{}",
            monomerized,
            std::str::from_utf8(&seq[..monomerized.unwrap_or(seq.len())]).unwrap()
        );
        debug!("------------------\n");
        monomerized
    }

    pub fn last_monomer_end_index_sensitive(&self, seq: &[u8]) -> Option<usize> {
        // First, we monomerize as normal
        let monomer_index = self.last_monomer_end_index(seq);
        let monomer = &seq[..monomer_index.unwrap_or(seq.len())];

        let rc = dna::revcomp(monomer);
        // debug!("monomer: {:?}", std::str::from_utf8(&rc).unwrap());
        let rc_monomer_index = self.first_monomer_end_index(&rc);
        // debug!("rc monomer index: {:?}", rc_monomer_index);
        match rc_monomer_index {
            None => monomer_index,
            Some(index) => Some(monomer_index.unwrap_or(seq.len()) - (monomer.len() - index)),
        }
    }

    /// A helper function to compute to get a slice of the monomer from a sequence.
    pub fn monomerize(self, seq: &[u8]) -> &[u8] {
        let end = self.last_monomer_end_index(seq);
        match end {
            None => seq,
            Some(end) => &seq[..end],
        }
    }

    pub fn monomerize_sensitive(self, seq: &[u8]) -> &[u8] {
        let end = self.last_monomer_end_index_sensitive(seq);
        match end {
            None => seq,
            Some(end) => &seq[..end],
        }
    }
}

#[cfg(test)]
mod test {
    use crate::monomerize::Monomerizer;
    use rstest::rstest;

    mod basic_tests {
        use super::*;
        use pretty_assertions::assert_eq;

        fn monomerize_with_seed_len_and_overlap_dist(
            seq: &str,
            seed_len: usize,
            overlap_dist: u64,
        ) -> &[u8] {
            let monomerizer = Monomerizer::builder()
                .seed_len(seed_len)
                .overlap_dist(overlap_dist)
                .build()
                .unwrap();
            monomerizer.monomerize(seq.as_bytes())
        }

        #[rstest]
        #[case("ATGCATGC", "ATGC", 4, 0, "identical repeat")]
        #[case("ATGCATGC", "ATGC", 4, 1, "identical repeat")]
        #[case("ATGCATGC", "ATGC", 2, 0, "identical repeat")]
        #[case("ATGCATGC", "ATGC", 2, 1, "identical repeat")]
        #[case("AAAAATTTTTAAAAATTTTT", "AAAAATTTTT", 10, 0, "repeat at beginning")]
        fn monomerize_with_seed_len_and_overlap_dist_works(
            #[case] seq: &str,
            #[case] expected: &str,
            #[case] seed_len: usize,
            #[case] overlap_dist: u64,
            #[case] name: &str,
        ) {
            let actual = monomerize_with_seed_len_and_overlap_dist(seq, seed_len, overlap_dist);
            assert_eq!(
                actual,
                expected.as_bytes(),
                "{}: seed_len = {}, overlap_dist = {}",
                name,
                seed_len,
                overlap_dist
            );
        }

        #[rstest]
        fn no_overlap(
            #[values(4, 5, 6, 7, 8, 9, 10)] seed_len: usize,
            #[values(0, 1, 2)] overlap_dist: u64,
        ) {
            assert_eq!(
                "TTTTTTTTTTTTAAAAAAAAAA".as_bytes(),
                monomerize_with_seed_len_and_overlap_dist(
                    "TTTTTTTTTTTTAAAAAAAAAA",
                    seed_len,
                    overlap_dist
                )
            );
        }

        #[rstest]
        fn mismatches_within_limit(
            #[values(4, 5, 6, 7, 8, 9, 10)] seed_len: usize,
            #[values(2, 3, 4)] overlap_dist: u64,
        ) {
            let x = b"TTAGCCCGTGTTTTATCGGAAGCTATCCTCAAAGCCCGTGTTTTATCGGAAGCTATCCTC";
            assert_eq!(
                Monomerizer::builder()
                    .seed_len(seed_len)
                    .overlap_dist(overlap_dist)
                    .build()
                    .unwrap()
                    .monomerize(x),
                b"TTAGCCCGTGTTTTATCGGAAGCTATCCTC"
            );
        }

        #[rstest]
        fn too_many_mismatches(
            #[values(4, 5, 6, 7, 8, 9, 10)] seed_len: usize,
            #[values(0, 1)] overlap_dist: u64,
        ) {
            let x = b"TTTTTGGTTTTTAAAAAAAAAATTTTTTTTTTTTAAAAAAAAAA";
            assert_eq!(
                Monomerizer::builder()
                    .seed_len(seed_len)
                    .overlap_dist(overlap_dist)
                    .build()
                    .unwrap()
                    .monomerize(x),
                x
            );
        }

        #[rstest]
        fn complete_identical_multimer(
            #[values(4, 5, 6, 7, 8, 9, 10)] seed_len: usize,
            #[values(0, 1, 2)] overlap_dist: u64,
        ) {
            assert_eq!(
                Monomerizer::builder()
                    .seed_len(seed_len)
                    .overlap_dist(overlap_dist)
                    .build()
                    .unwrap()
                    .monomerize(b"AAAAATTTTTAAAAATTTTTAAAAATTTTT"),
                b"AAAAATTTTT"
            );
        }

        #[rstest]
        /// This test checks to make sure that the overlap_dist method is computed per monomer, not globally.
        fn multimer_with_mismatch_in_each(
            #[values(4, 5, 6, 7)] seed_len: usize,
            #[values(1, 2, 3)] overlap_dist: u64, // 0 is not allowed since there's a mismatch in each
        ) {
            assert_eq!(
                Monomerizer::builder()
                    .seed_len(seed_len)
                    .overlap_dist(overlap_dist)
                    .build()
                    .unwrap()
                    .monomerize(b"AACAATTTTTAAGAATTTTTAAAAATTTTT"),
                //                11 111111122 22222223333333333
                //                     ^^^^^     ^^^^^    ^^^^^
                b"AACAATTTTT"
            );
        }

        #[test]
        fn multimer_with_mismatch_in_middle() {
            assert_eq!(
                Monomerizer::builder()
                    .seed_len(5)
                    .overlap_dist(1)
                    .build()
                    .unwrap()
                    .monomerize(b"AAAAATTTTTAAGAATTTTTAAAAATTTTT"),
                //                111111111122 22222223333333333
                //                     ^^^^^     ^^^^^    ^^^^^
                b"AAAAATTTTT"
            );
        }
        // #[test]
        // fn multimer_with_partial_mer() {
        //     assert_eq!(
        //         Monomerizer::builder()
        //             .seed_len(5)
        //             .overlap_dist(1)
        //             .build()
        //             .unwrap()
        //             .monomerize(b"TATTTTTAAGAATTTTTAAAAATTTTT"),
        //         //                111111122 22222223333333333
        //         //                  ^^^^^     ^^^^^    ^^^^^
        //         b"AAAAATTTTT"
        //     );
        // }

        #[test]
        fn multimer_with_seed_repeated() {
            assert_eq!(
                Monomerizer::builder()
                    .seed_len(4)
                    .overlap_dist(0)
                    .build()
                    .unwrap()
                    .monomerize(b"TGCCAATGCATGCCAATGC"),
                //                     ^^^^      ^^^^
                //                         ^^^^ <- seed repeated but not a valid overlap
                b"TGCCAATGCA"
            );
        }

        #[rstest]
        fn big_dimer(
            #[values(6, 7, 8, 9, 10, 11)] seed_len: usize,
            #[values(0, 1, 3)] overlap_dist: u64,
        ) {
            let input = b"ATGACAGGTACAGCATAATGACAGGTACAGCATA";
            let output = b"ATGACAGGTACAGCATA";

            let m = Monomerizer::builder()
                .seed_len(seed_len)
                .overlap_dist(overlap_dist)
                .build()
                .unwrap();
            assert_eq!(
                m.monomerize(input),
                output,
                "seed_len = {}, overlap_dist = {}",
                seed_len,
                overlap_dist
            );
        }

        #[rstest]
        fn ambivirus(
            #[values(10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 32, 33, 63)] seed_len: usize,
            #[values(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)] overlap_dist: u64,
        ) {
            // MT354566.2 Cryphonectria parasitica ambivirus 1 strain ACP34, complete genome
            let input = b"CCCCTAATGCCGGGGTACTGATGAAGCCATGGCAGGCCGAAAGCAAATCTCATATTATGGCATACGATACGTGTGCGTCAGTTCCCCACCTGACTCTCGTCGGGTGTGGCTGAATGCAACAGTATCTTATGCGCAGGCATGATGACTGTCTTTCGCGGCCATGTTCTTCCATCAGCGCTAGGGCATTAGGGAGTAGCGATTTCGTTCTACCATGAACTCATCAGGCACCCGGTATCAGGGCACTATCGCACTACAGGCATGTCAATTAACGAACGACGAACTACGAAAAGATACGACATCGATCCTAAGAAAGTACACTATGGAACCAGGTCATCAACCAATCAATGTTGCCAGGAGGTCATCAAGATCAACCTCACCCTCACCAGCATCCTCTTCACGACGAGCTCTCTTCTTCTTTGGACCAGCTTCTCCAGCAATCTTGACCAACGACAGCCAAATCCGCTTCCTGAATTCCTCACTCTCCACCACCATATTCCGATATTTCCCAGCACGTTCATGGCTGTCTACCTTGATTGCCTGGCGTTCGAGCATCTTCTCAAAATCCTGGATGGCGATACTGAGAGGCTTCGGAGTAATGAAGAGTTTATCAGGCATGTTGGCAAAAGCATCCATCTCCCCAAGCTTACGTTTCCCACCCGAAACAACAGGTTTCTCTTTCGCTTCACACCGAATCTCGACACCTCGTTCATTCCAAAATGAGATCGCCTCAGGACCGAACCGAGAAAGGATGATAAACGCTCGACTCGAAAGGTTGACACGAATGGATGGAATGTACGTGGGACGACTGAGCTTAAGGACAGTATCAGAAGCAAGTGTCTCGAGCGCATCCGCAATCATCTGAGGGTTCCGAACGATGTACCCGGCACCCATATAGTTGAGGCTTCGAAGACTCTGGTTGACCTCCTTTTCGACCACATCATCCAACTCAGTCAGACGCACACCAGAGAAGACCTCAATGAGATTCTCAGGTTCAGTAAACGTCTTCCGATCGACCGATGGACCCGTATAAGTGCCTTCAGCAGTCATCCTCCCAAAAGCAGAGACCAGATTCGAGACAGCTGATTCGTGGGGATCCAAGCTGAGGAGTGCAAGAGAAAGTTGTTCTTCAGTACTGGGGGCAATCCACTTGGTCTGGTCGAAGATTGAGAATCGTGCGCCAAGCAGAACGAATCCCTTATACTCAGCAGCCTCGAAGATCATAAAGCATCTGCTCTGTGTATCAAGTGCTAAATCGATCCCCTTAAGTAGGTGAGAGAGCATCATCCCTTCCTTGGTCGTTGACAGAGATCCACAGCCGCGTCGAATAGTCAGAAAGGTGTCTTGAGACTTCTGCACCGTCGATCCGAGAAGTTGAAAGAACCGTTGAAGGATCACATTCGACATAAACACCGCGTCAGGCTGGACAAGACCGTGGAAGTATGGAAACGCAATTCCAGGAAGGTTGGGGACGTCGCGAGATGAGCCAATATTGGTCTGTGGCTTGCGACACGGTTTGGCCTTGTGGACAGTATCGTTTCCTGACGTGGAGAATCGTCCCTGAGGGGTTTCAAGTCGTGGATATGTCGTGGCGTCGTCGTCAGTCATCTGGTCAGGTGTGTACATCTCCCGTCCTTGAAGGGTTGACGAAAAGGCAATAGGAGCGAACGCGGATACATCCTTTCCTTCCTTGTCCTTCGTGTTGTAAGGGATATACTCAGTGTCCTTTACGTCATGAAGTGTGTGGAACTCTTGGATCGACTCGTCAGATTCGTATGAGAAGGCGGGGTATGAATGAGTACGCAAGAAGGTGGCGCAGACTGAGGACAAGAACTTGTACGCAAAGGCACCAGCTGGGTTGACGGCATCTCTACCACAGTTGAAATGGGCTGGAGTCAAGGCGATGATCGTATCCAACATGCCGTCAGGAAGTGCAGGTGCATGATAGATAGTCGTTTTACCAGAGGTACTCGTGGAGGCTCCAGGACCGAGACTGTTGACTAGCCTGATGGATTTGTCGGAATCGACCGCTTTCTTGAGCGTAAACTGCACGGGAGGACAGTTCAGAGGCAGAACATCTGATGAGTCAAAGAGTCCGAGGGGGAGTCCATTGACTAGAGTCGTTGCCTTGGTCTTTGTGGCCCAGAGGAAACAGGAGATGTTCGTCGCATAAGCAGTGTAGGCAGTCTGTTTGGTTAGGGCGACTCCGGCGAGATTGAGTGTTGAGATAGGCGAGGGGGCGAAGTTGTCCAATCCATGTAGTTTGAAGAGGGGAACGGTGTTCTGGCTGCCGGCAACGTGCCGAGCTGTGTAGATAGAGGAGGTCATGATGTCTGGGGACTGGAGGTGTGTAGTGTTGGCTCTTTTGTACCATCCCACTACACAGGCGAGATCTGTCGATCGACCTACCCACCAGCCCACTCAAGGACACCCCAGATCACTGCGTCCTGTGATACGGGCTTGTTGTGTGTCGATGTTCAAAAGAATGAAGCACACAGCCGCAACGTACTGGACCTGGTCTTGGACCTTCTCTCACCACACTGCTCTATGCACTGTCGTACAGACCATCCCATCCCTATCCTTACTACAAAGGAAAGACATACTCGGCCTCTGCAGGTCTCAAACACCGATCATACCAGGTTCCCACCAACACCGACAACTAAAGTGCTTAACCGACATCATCCAAGTTCCATCTGAACCTCAAAGAAGGGAAAACGACATCCAGTACTTGGAATACTCCTTCCATCGCTTTCTCAGAGACTACCAGGATGTATACCTTGACCATCTCTCCGCTGCTGTTATCGCTTTATCAGTCCGCGTCGCGCCCGGTATATCCGATTCAGACTTTCAGAACAGTCTCTCGTCCTCGATTTCCAAAGAGATTTCCGCTGCTATAACGATCTCCAACCAAGCACGCTTCGATTTCCGTTCGGTCATGGATTCGTCATTATCGTCCCGATTCATCCCTTTCACTGATAGACTGCCGCCGTCTCCAGAAGTCGCCTCTGATACGTATATCGCTCTTTGTCGAGTCCGCTCCCGACTTGATAATCGATTCCTACTGCCAACTCTTGGGAATAAGACCGTCGGCCAGGTTGTTCGTCATAACTTCTATGACATATATGACTCTTTCGGTTTAACCTTCAAACATGAGGGGAGGAGGAATGATGATACCGTGTCGACTTCGGATTGTATGCGTCTGTACCTCGAGACCGGAGTGTATCCACACGGACCTGTTGAAATGCGACGAGCCTGGACATATAATCAACTTGACCCAAGGGTGTATTACGCACGTGGAGGAGACGTTATGCACACATCCCAATACGTGCAGTCAATCGCCAATATGCTTATCGACGCTTTTCCCGAAACACACCGGAAAGATCGCTTCATGCCACCACGGGATCCACTAGCTGATGACGATGTTGAAGTCATATATGATTATTCATCGTTCACTTCGACTCTCGATTCGGTGGTTCCGTTTCTCGACAATCTCGCTGAGTTCTTCCGTGGTACTGTGGTTCATCTTGTCGACTTCAGGAACGGAGTTGTCCCAACCGACTTAGGTGATCTTATCGCACGATACAACACAGAGTGTAATCTCTATGGTACCTTCGACGCATCTAGAGTTCTGGGTCAATCTGCCGGTACGACTCTTCTTCAACATACCTGTGGAATGCTCGGAGTAGAAGGGAACATCTTCTTTGCCACACTCTTGCATGGGATACACCTCCGATTCATTGCTGGCCTGAACCGATCTCGTTGCGTTGGTGATGACGCTAGAATGCATCATAGGGTACCTTTCGGCATCATGGACAACACTGAAACCGATTACCTCGCTTGGGTTCTTGCTGGGTGTGGCGATCTGAGCAAGGAAAAGATGGGGAAATTCGAATCGGGTGTCGACAGTGAGCTCCAGGCCTATCGCTATATTAAGAGGCCGATACATCGCGACGGTTCCATTATGATTGAGGGGATACTCCTAACACTTCCTTCCATCATACCACTCCTCGGGGCCCAAGATCGATTCCACACTGTCACTCCATCCGTCAGTCATCCCTCGAGAAGAACGTACTCACAGATTCTCAGGTTCATCCAGGAGCTTTTCATTCATGGGTTATGTTACGACAGTGACGACGTTTCGTGGAAGTCGATATTGAAACATTTGATGTTCCTTAGGAGGTTATGTATTGCTGAGGATCCGGACTTTGAGCACTCCATGTTCATGAACTCATCTTATCTGACGAAGTATCGGTTTCCGCCACCAGAGACTTGGGGGAAGATGCCGATCACAGACTGGGTAGTAGGCGACATCATGTACGACGAGGTTATACGGTTTCCGATGAAGGGAGAGAGAGAATCAGAAGGAGGGTGTGACGGTCGCGTTGGTTCGGAAATGCTGAGGATGGGTTCCAAGGCTAGAGGATGTCTCGTGAAGCTAGGGTATCTGGAGGAGGAGAAGATGTTCGATGACGTTTCTGTCAAACTTGTCGGTTTGGACCTGTTCCTAGAGTACCTTGGAGGAAGATACCGTTCTATCAGCAAGTTCGTTGTCGTTAAGGACATACCTGGGTGGGTAGCACAAGTACCGAGTAGTCTATGAATGATGATAGAGCGAGATTTGCTATACCCCTAATGCCGGGGTACTGATGAAGCCATGGCAGGCCGAAAGCAAATCTCATATTATGGCATACGATACGTGTGCGTCAGTTCCCCACCTGACTCTCGTCGGGTGTGGCTGAATGCAACAGTATCTTATGCGCAGGCATGATGACTGTCTTTCGCGGCCATGTTCTTCCATCAGCGCTAGGGCATTAGGGAGTAGCGATTTCGTTCTACCATGAACTCATCAGGCACCCGGTATCAGGGCAC";
            let monomer = b"CCCCTAATGCCGGGGTACTGATGAAGCCATGGCAGGCCGAAAGCAAATCTCATATTATGGCATACGATACGTGTGCGTCAGTTCCCCACCTGACTCTCGTCGGGTGTGGCTGAATGCAACAGTATCTTATGCGCAGGCATGATGACTGTCTTTCGCGGCCATGTTCTTCCATCAGCGCTAGGGCATTAGGGAGTAGCGATTTCGTTCTACCATGAACTCATCAGGCACCCGGTATCAGGGCACTATCGCACTACAGGCATGTCAATTAACGAACGACGAACTACGAAAAGATACGACATCGATCCTAAGAAAGTACACTATGGAACCAGGTCATCAACCAATCAATGTTGCCAGGAGGTCATCAAGATCAACCTCACCCTCACCAGCATCCTCTTCACGACGAGCTCTCTTCTTCTTTGGACCAGCTTCTCCAGCAATCTTGACCAACGACAGCCAAATCCGCTTCCTGAATTCCTCACTCTCCACCACCATATTCCGATATTTCCCAGCACGTTCATGGCTGTCTACCTTGATTGCCTGGCGTTCGAGCATCTTCTCAAAATCCTGGATGGCGATACTGAGAGGCTTCGGAGTAATGAAGAGTTTATCAGGCATGTTGGCAAAAGCATCCATCTCCCCAAGCTTACGTTTCCCACCCGAAACAACAGGTTTCTCTTTCGCTTCACACCGAATCTCGACACCTCGTTCATTCCAAAATGAGATCGCCTCAGGACCGAACCGAGAAAGGATGATAAACGCTCGACTCGAAAGGTTGACACGAATGGATGGAATGTACGTGGGACGACTGAGCTTAAGGACAGTATCAGAAGCAAGTGTCTCGAGCGCATCCGCAATCATCTGAGGGTTCCGAACGATGTACCCGGCACCCATATAGTTGAGGCTTCGAAGACTCTGGTTGACCTCCTTTTCGACCACATCATCCAACTCAGTCAGACGCACACCAGAGAAGACCTCAATGAGATTCTCAGGTTCAGTAAACGTCTTCCGATCGACCGATGGACCCGTATAAGTGCCTTCAGCAGTCATCCTCCCAAAAGCAGAGACCAGATTCGAGACAGCTGATTCGTGGGGATCCAAGCTGAGGAGTGCAAGAGAAAGTTGTTCTTCAGTACTGGGGGCAATCCACTTGGTCTGGTCGAAGATTGAGAATCGTGCGCCAAGCAGAACGAATCCCTTATACTCAGCAGCCTCGAAGATCATAAAGCATCTGCTCTGTGTATCAAGTGCTAAATCGATCCCCTTAAGTAGGTGAGAGAGCATCATCCCTTCCTTGGTCGTTGACAGAGATCCACAGCCGCGTCGAATAGTCAGAAAGGTGTCTTGAGACTTCTGCACCGTCGATCCGAGAAGTTGAAAGAACCGTTGAAGGATCACATTCGACATAAACACCGCGTCAGGCTGGACAAGACCGTGGAAGTATGGAAACGCAATTCCAGGAAGGTTGGGGACGTCGCGAGATGAGCCAATATTGGTCTGTGGCTTGCGACACGGTTTGGCCTTGTGGACAGTATCGTTTCCTGACGTGGAGAATCGTCCCTGAGGGGTTTCAAGTCGTGGATATGTCGTGGCGTCGTCGTCAGTCATCTGGTCAGGTGTGTACATCTCCCGTCCTTGAAGGGTTGACGAAAAGGCAATAGGAGCGAACGCGGATACATCCTTTCCTTCCTTGTCCTTCGTGTTGTAAGGGATATACTCAGTGTCCTTTACGTCATGAAGTGTGTGGAACTCTTGGATCGACTCGTCAGATTCGTATGAGAAGGCGGGGTATGAATGAGTACGCAAGAAGGTGGCGCAGACTGAGGACAAGAACTTGTACGCAAAGGCACCAGCTGGGTTGACGGCATCTCTACCACAGTTGAAATGGGCTGGAGTCAAGGCGATGATCGTATCCAACATGCCGTCAGGAAGTGCAGGTGCATGATAGATAGTCGTTTTACCAGAGGTACTCGTGGAGGCTCCAGGACCGAGACTGTTGACTAGCCTGATGGATTTGTCGGAATCGACCGCTTTCTTGAGCGTAAACTGCACGGGAGGACAGTTCAGAGGCAGAACATCTGATGAGTCAAAGAGTCCGAGGGGGAGTCCATTGACTAGAGTCGTTGCCTTGGTCTTTGTGGCCCAGAGGAAACAGGAGATGTTCGTCGCATAAGCAGTGTAGGCAGTCTGTTTGGTTAGGGCGACTCCGGCGAGATTGAGTGTTGAGATAGGCGAGGGGGCGAAGTTGTCCAATCCATGTAGTTTGAAGAGGGGAACGGTGTTCTGGCTGCCGGCAACGTGCCGAGCTGTGTAGATAGAGGAGGTCATGATGTCTGGGGACTGGAGGTGTGTAGTGTTGGCTCTTTTGTACCATCCCACTACACAGGCGAGATCTGTCGATCGACCTACCCACCAGCCCACTCAAGGACACCCCAGATCACTGCGTCCTGTGATACGGGCTTGTTGTGTGTCGATGTTCAAAAGAATGAAGCACACAGCCGCAACGTACTGGACCTGGTCTTGGACCTTCTCTCACCACACTGCTCTATGCACTGTCGTACAGACCATCCCATCCCTATCCTTACTACAAAGGAAAGACATACTCGGCCTCTGCAGGTCTCAAACACCGATCATACCAGGTTCCCACCAACACCGACAACTAAAGTGCTTAACCGACATCATCCAAGTTCCATCTGAACCTCAAAGAAGGGAAAACGACATCCAGTACTTGGAATACTCCTTCCATCGCTTTCTCAGAGACTACCAGGATGTATACCTTGACCATCTCTCCGCTGCTGTTATCGCTTTATCAGTCCGCGTCGCGCCCGGTATATCCGATTCAGACTTTCAGAACAGTCTCTCGTCCTCGATTTCCAAAGAGATTTCCGCTGCTATAACGATCTCCAACCAAGCACGCTTCGATTTCCGTTCGGTCATGGATTCGTCATTATCGTCCCGATTCATCCCTTTCACTGATAGACTGCCGCCGTCTCCAGAAGTCGCCTCTGATACGTATATCGCTCTTTGTCGAGTCCGCTCCCGACTTGATAATCGATTCCTACTGCCAACTCTTGGGAATAAGACCGTCGGCCAGGTTGTTCGTCATAACTTCTATGACATATATGACTCTTTCGGTTTAACCTTCAAACATGAGGGGAGGAGGAATGATGATACCGTGTCGACTTCGGATTGTATGCGTCTGTACCTCGAGACCGGAGTGTATCCACACGGACCTGTTGAAATGCGACGAGCCTGGACATATAATCAACTTGACCCAAGGGTGTATTACGCACGTGGAGGAGACGTTATGCACACATCCCAATACGTGCAGTCAATCGCCAATATGCTTATCGACGCTTTTCCCGAAACACACCGGAAAGATCGCTTCATGCCACCACGGGATCCACTAGCTGATGACGATGTTGAAGTCATATATGATTATTCATCGTTCACTTCGACTCTCGATTCGGTGGTTCCGTTTCTCGACAATCTCGCTGAGTTCTTCCGTGGTACTGTGGTTCATCTTGTCGACTTCAGGAACGGAGTTGTCCCAACCGACTTAGGTGATCTTATCGCACGATACAACACAGAGTGTAATCTCTATGGTACCTTCGACGCATCTAGAGTTCTGGGTCAATCTGCCGGTACGACTCTTCTTCAACATACCTGTGGAATGCTCGGAGTAGAAGGGAACATCTTCTTTGCCACACTCTTGCATGGGATACACCTCCGATTCATTGCTGGCCTGAACCGATCTCGTTGCGTTGGTGATGACGCTAGAATGCATCATAGGGTACCTTTCGGCATCATGGACAACACTGAAACCGATTACCTCGCTTGGGTTCTTGCTGGGTGTGGCGATCTGAGCAAGGAAAAGATGGGGAAATTCGAATCGGGTGTCGACAGTGAGCTCCAGGCCTATCGCTATATTAAGAGGCCGATACATCGCGACGGTTCCATTATGATTGAGGGGATACTCCTAACACTTCCTTCCATCATACCACTCCTCGGGGCCCAAGATCGATTCCACACTGTCACTCCATCCGTCAGTCATCCCTCGAGAAGAACGTACTCACAGATTCTCAGGTTCATCCAGGAGCTTTTCATTCATGGGTTATGTTACGACAGTGACGACGTTTCGTGGAAGTCGATATTGAAACATTTGATGTTCCTTAGGAGGTTATGTATTGCTGAGGATCCGGACTTTGAGCACTCCATGTTCATGAACTCATCTTATCTGACGAAGTATCGGTTTCCGCCACCAGAGACTTGGGGGAAGATGCCGATCACAGACTGGGTAGTAGGCGACATCATGTACGACGAGGTTATACGGTTTCCGATGAAGGGAGAGAGAGAATCAGAAGGAGGGTGTGACGGTCGCGTTGGTTCGGAAATGCTGAGGATGGGTTCCAAGGCTAGAGGATGTCTCGTGAAGCTAGGGTATCTGGAGGAGGAGAAGATGTTCGATGACGTTTCTGTCAAACTTGTCGGTTTGGACCTGTTCCTAGAGTACCTTGGAGGAAGATACCGTTCTATCAGCAAGTTCGTTGTCGTTAAGGACATACCTGGGTGGGTAGCACAAGTACCGAGTAGTCTATGAATGATGATAGAGCGAGATTTGCTATA";
            assert_eq!(
                Monomerizer::builder()
                    .seed_len(seed_len)
                    .overlap_dist(overlap_dist)
                    .build()
                    .unwrap()
                    .monomerize(input),
                monomer,
                "seed_len: {}, overlap_dist: {}",
                seed_len,
                overlap_dist
            );
        }

        #[rstest]
        #[case("TCCTCCATCACCTAGTTTATGTAGAAACGCTGCTAAATCAATTTCCTCCATCACCTAGTTTATGTAGAAAAGCTGCTAAATCAATTTCCTCCATCACCTAGTTTATGTAGAAACGCTGCTAAATCAATTTCCTCCATCACCTAGTTTATGTAGAAAAGCTGCTA", "TCCTCCATCACCTAGTTTATGTAGAAACGCTGCTAAATCAATT")]
        #[case("GCAGTTATAGAGAGAGTGGGTCAGTTCATTATTACACTGCAGTAATAGAGAGAGTGGGTCAGTTCATTATTACACTGCAGTGATAGAGAGAGTGGGTCAGTTCATTATTACACTGCAGTTATAGAGAGAGTGGGTCAGTTCATTATTACACTGCAGTAATAGAGAGAGTGGGTCAGTTCATTATTACACTGCAGTGATAG", "GCAGTTATAGAGAGAGTGGGTCAGTTCATTATTACACT")]
        #[case("CTGGCCCAGGGGCTTCTAGTCAAACAGGCCTCTCTTCCCCACTCCTTACCTCTTCTGGTCTCTGGCCCAGGGGCTTCTAGTCAAACAGGCCTCTCTTCCCCACTCCTTACCTCTTCTGGTCTCTGGCCCCTGGCCCAGGGGCTTCTAGTCAAACAGGCCTCTCTTCCCCACTCCTTACCTCTTCTGGTCTCTGGCCCAGGGGCTTCT", "CTGGCCCAGGGGCTTCTAGTCAAACAGGCCTCTCTTCCCCACTCCTTACCTCTTCTGGTCT")]
        fn single_pass_regressions(#[case] input: &str, #[case] known_monomer: &str) {
            // monomerize the input
            let m = Monomerizer::builder()
                .overlap_min_identity(0.95)
                .seed_len(10)
                .build()
                .unwrap();
            let monomer = m.monomerize(input.as_bytes());
        }
    }

    mod validation {
        use super::*;
        #[test]
        #[should_panic(expected = "seed_len")]
        fn seed_not_set() {
            Monomerizer::builder().build().unwrap().monomerize(b"");
        }

        #[test]
        #[should_panic(expected = "at least 1 and at most 63")]
        fn seed_too_long() {
            Monomerizer::builder()
                .seed_len(100)
                .build()
                .unwrap()
                .monomerize(b"");
        }

        #[test]
        #[should_panic(expected = "at least 1 and at most 63")]
        fn seed_too_short() {
            Monomerizer::builder()
                .seed_len(0)
                .build()
                .unwrap()
                .monomerize(b"");
        }
    }

    mod overlap_percent {
        use super::*;
        use pretty_assertions::assert_eq;

        #[test]
        fn dimer_with_overlap_percentage() {
            let input = b"ATGCCCATGCGCCAGCGCAGATGCGAATGCGCCAGCGCAG";
            let output = b"ATGCCCATGCGCCAGCGCAG";

            // ATGCCCATGCGCCAGCGCAG Overlap (20 nt)
            // ||||..||||||||||||||
            // ATGCGAATGCGCCAGCGCAG Monomer (20 nt)

            let m = Monomerizer::builder()
                .overlap_min_identity(0.95)
                .seed_len(4)
                .build()
                .unwrap();
            assert_eq!(m.monomerize(input), input);

            let m = Monomerizer::builder()
                .overlap_min_identity(0.90)
                .seed_len(4)
                .build()
                .unwrap();
            assert_eq!(m.monomerize(input), output);
        }

        #[test]
        fn overlap_percentage_rounds_down_to_nearest_nt() {
            let input = b"TGCCCATGCGCCAGCGCAGATGCGAATGCGCCAGCGCAG";
            let output = b"TGCCCATGCGCCAGCGCAGA";

            //  TGCCCATGCGCCAGCGCAG Overlap (19 nt)
            //  |||..||||||||||||||
            // ATGCGAATGCGCCAGCGCAG Monomer (20 nt)

            let m = Monomerizer::builder()
                .overlap_min_identity(0.95) // 18.05/19 nt identity required
                .seed_len(4)
                .build()
                .unwrap();
            assert_eq!(m.monomerize(input), input);

            let m = Monomerizer::builder()
                .overlap_min_identity(0.90) // 17.1/19 nt identity required
                .seed_len(4)
                .build()
                .unwrap();
            assert_eq!(m.monomerize(input), output);

            let m = Monomerizer::builder()
                .overlap_min_identity(0.94) // 17.86/19 nt identity required
                .seed_len(4)
                .build()
                .unwrap();
            assert_eq!(m.monomerize(input), output);
        }

        #[test]
        #[should_panic(expected = "overlap_dist and overlap_min_identity")]
        fn overlap_percentage_and_dist_panics() {
            let input = b"TGCCCATGCGCCAGCGCAGATGCGAATGCGCCAGCGCAG";

            let m = Monomerizer::builder()
                .overlap_dist(1)
                .overlap_min_identity(0.95)
                .seed_len(4)
                .build()
                .unwrap();

            m.monomerize(input);
        }
    }

    mod sensitive {
        use crate::canonicalize;

        use super::*;
        use pretty_assertions::{assert_eq, assert_ne};
    }

    mod fuzzing {

        use super::*;
        use proptest::prelude::*;

        proptest! {
                #[test]
                fn concatenated_always_monomerizes(input in "[ACGT]{12,100}") {
                    let concatenated = format!("{}{}{}", input, input, input);
                    let m = Monomerizer::builder()
                        .overlap_min_identity(0.95)
                        .seed_len(10)
                        .build()
                        .unwrap();
                    prop_assert_eq!(m.clone().monomerize(concatenated.as_bytes()), input.as_bytes());
                    // prop_assert_eq!(m.clone().monomerize_sensitive(concatenated.as_bytes()), input.as_bytes());
                }
                #[test]
                fn small_mutations_outside_seed_still_monomerize(input in "[ACGT]{100,200}", mutation_index in 10..90usize) {
                    let mut concat = format!("{}{}", input, input).into_bytes();

                    concat[mutation_index] = match concat[mutation_index] {
                        b'A' => b'C',
                        b'C' => b'G',
                        b'G' => b'T',
                        b'T' => b'A',
                        _ => unreachable!(),
                    };

                    let m = Monomerizer::builder()
                        .overlap_min_identity(0.95)
                        .seed_len(10)
                        .build()
                        .unwrap();
                    prop_assert_eq!(m.monomerize(&concat), &concat[..input.len()]);
                }
        }
    }
}
