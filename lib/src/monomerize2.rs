use bio::alignment::distance::simd::*;
use bio::pattern_matching::shift_and;
use log::warn;
use memchr::memmem;

#[derive(Builder, Default, Clone)]
#[builder(setter(strip_option), build_fn(validate = "Self::validate"))]
pub struct Monomerizer {
    /// The maximum number of mismatches allowed in an overlap.
    #[builder(default)]
    pub overlap_dist: Option<u64>,
    /// The minimum percent identity within an overlap that may be considered a match. Overrides `overlap_dist`.
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
                1..=64 => {}
                _ => {
                    return Err(format!(
                        "Seed length must be at least 1 and at most 64 but was set to {}.",
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
    pub fn monomerize(self, seq: &[u8]) -> &[u8] {
        // if the sequence is shorter than the seed, give up
        let seed_len = self.seed_len;
        if seq.len() < seed_len {
            warn!("Sequence is less than seed length");
            return seq;
        }

        // slice last n bases of the record
        let seed = &seq[seq.len() - seed_len..];

        // we start out assuming the seed isn't even present anywhere else
        let mut last_match: Option<usize> = None;

        let text = &seq[..seq.len() - seed_len];

        // create the matcher
        let matcher = shift_and::ShiftAnd::new(seed);
        // let matcher = memmem::find_iter(text, seed);

        // loop over the searcher (noting that the final n bases are not included)
        for occ in matcher.find_all(text) {
            // slice the potential overlap from the sequence
            let potential_overlap = &seq[last_match.unwrap_or(0)..occ + seed.len()];

            // extend backward from the end of the sequence to the start of the potential overlap
            let seed_extension = &seq[seq.len() - potential_overlap.len()..];

            // compare the potential overlap to the seed
            let dist = hamming(potential_overlap, seed_extension);

            // decide whether the overlap is good enough to be a monomer
            let max_dist = match self.overlap_min_identity {
                Some(identity) => {
                    potential_overlap.len() as u64
                        - (potential_overlap.len() as f64 * identity).floor() as u64
                }
                None => self.overlap_dist.unwrap_or(0),
            };

            if dist <= max_dist {
                last_match = Some(occ);
            }
        }

        match last_match {
            Some(last_match) => &seq[last_match + seed.len()..],
            None => seq,
        }
    }
}

#[cfg(test)]
mod test {
    use crate::monomerize2::Monomerizer;
    use pretty_assertions::assert_eq;

    #[test]
    fn no_overlap_at_all() {
        assert_eq!(
            Monomerizer::builder()
                .seed_len(10)
                .overlap_dist(1)
                .build()
                .unwrap()
                .monomerize(b"TTTTTTTTTTTTAAAAAAAAAA"),
            b"TTTTTTTTTTTTAAAAAAAAAA"
        );
    }
    #[test]
    fn identical_repeat() {
        assert_eq!(
            Monomerizer::builder()
                .seed_len(2)
                .overlap_dist(1)
                .build()
                .unwrap()
                .monomerize(b"ATGCATGC"),
            b"ATGC"
        );
    }
    #[test]
    fn big_with_mismatch_overlap() {
        let x = b"TTTTTGGTTTTTAAAAAAAAAATTTTTTTTTTTTAAAAAAAAAA";
        assert_eq!(
            Monomerizer::builder()
                .seed_len(10)
                .overlap_dist(3)
                .build()
                .unwrap()
                .monomerize(x),
            b"TTTTTTTTTTTTAAAAAAAAAA"
        );
        assert_eq!(
            Monomerizer::builder()
                .seed_len(10)
                .overlap_dist(2)
                .build()
                .unwrap()
                .monomerize(x),
            b"TTTTTTTTTTTTAAAAAAAAAA"
        );
        assert_eq!(
            Monomerizer::builder()
                .seed_len(10)
                .overlap_dist(1)
                .build()
                .unwrap()
                .monomerize(x),
            x
        );
        assert_eq!(
            Monomerizer::builder()
                .seed_len(10)
                .overlap_dist(0)
                .build()
                .unwrap()
                .monomerize(x),
            x
        );
    }
    #[test]
    fn repeat_at_beginning() {
        assert_eq!(
            Monomerizer::builder()
                .seed_len(10)
                .overlap_dist(0)
                .build()
                .unwrap()
                .monomerize(b"AAAAATTTTTAAAAATTTTT"),
            b"AAAAATTTTT"
        );
    }
    #[test]
    fn multimer() {
        assert_eq!(
            Monomerizer::builder()
                .seed_len(8)
                .overlap_dist(0)
                .build()
                .unwrap()
                .monomerize(b"AAAAATTTTTAAAAATTTTTAAAAATTTTT"),
            b"AAAAATTTTT"
        );
    }
    #[test]
    fn multimer_with_mismatch_in_each() {
        assert_eq!(
            Monomerizer::builder()
                .seed_len(5)
                .overlap_dist(1)
                .build()
                .unwrap()
                .monomerize(b"AACAATTTTTAAGAATTTTTAAAAATTTTT"),
            //                11 111111122 22222223333333333
            //                     ^^^^^     ^^^^^    ^^^^^
            b"AAAAATTTTT"
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
    #[test]
    fn multimer_with_partial_mer() {
        assert_eq!(
            Monomerizer::builder()
                .seed_len(5)
                .overlap_dist(1)
                .build()
                .unwrap()
                .monomerize(b"TATTTTTAAGAATTTTTAAAAATTTTT"),
            //                111111122 22222223333333333
            //                  ^^^^^     ^^^^^    ^^^^^
            b"AAAAATTTTT"
        );
    }

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
            b"ATGCCAATGC"
        );
    }

    #[test]
    fn dimer() {
        let input = b"TAAAAAAAAAAAAATAAAAAAAAAAAAA";
        let output = b"TAAAAAAAAAAAAA";

        // for seed_len in 6..=10 {
        let m = Monomerizer::builder()
            .overlap_dist(0)
            .seed_len(6)
            .build()
            .unwrap();
        assert_eq!(m.monomerize(input), output);
        println!("\n");
        // }
    }

    #[test]
    #[should_panic(expected = "seed_len")]
    fn seed_not_set() {
        Monomerizer::builder().build().unwrap().monomerize(b"");
    }

    #[test]
    #[should_panic(expected = "at least 1 and at most 64")]
    fn seed_too_long() {
        Monomerizer::builder()
            .seed_len(100)
            .build()
            .unwrap()
            .monomerize(b"");
    }

    #[test]
    #[should_panic(expected = "at least 1 and at most 64")]
    fn seed_too_short() {
        Monomerizer::builder()
            .seed_len(0)
            .build()
            .unwrap()
            .monomerize(b"");
    }

    #[test]
    fn ambivirus() {
        // MT354566.2 Cryphonectria parasitica ambivirus 1 strain ACP34, complete genome
        let input = b"CCCCTAATGCCGGGGTACTGATGAAGCCATGGCAGGCCGAAAGCAAATCTCATATTATGGCATACGATACGTGTGCGTCAGTTCCCCACCTGACTCTCGTCGGGTGTGGCTGAATGCAACAGTATCTTATGCGCAGGCATGATGACTGTCTTTCGCGGCCATGTTCTTCCATCAGCGCTAGGGCATTAGGGAGTAGCGATTTCGTTCTACCATGAACTCATCAGGCACCCGGTATCAGGGCACTATCGCACTACAGGCATGTCAATTAACGAACGACGAACTACGAAAAGATACGACATCGATCCTAAGAAAGTACACTATGGAACCAGGTCATCAACCAATCAATGTTGCCAGGAGGTCATCAAGATCAACCTCACCCTCACCAGCATCCTCTTCACGACGAGCTCTCTTCTTCTTTGGACCAGCTTCTCCAGCAATCTTGACCAACGACAGCCAAATCCGCTTCCTGAATTCCTCACTCTCCACCACCATATTCCGATATTTCCCAGCACGTTCATGGCTGTCTACCTTGATTGCCTGGCGTTCGAGCATCTTCTCAAAATCCTGGATGGCGATACTGAGAGGCTTCGGAGTAATGAAGAGTTTATCAGGCATGTTGGCAAAAGCATCCATCTCCCCAAGCTTACGTTTCCCACCCGAAACAACAGGTTTCTCTTTCGCTTCACACCGAATCTCGACACCTCGTTCATTCCAAAATGAGATCGCCTCAGGACCGAACCGAGAAAGGATGATAAACGCTCGACTCGAAAGGTTGACACGAATGGATGGAATGTACGTGGGACGACTGAGCTTAAGGACAGTATCAGAAGCAAGTGTCTCGAGCGCATCCGCAATCATCTGAGGGTTCCGAACGATGTACCCGGCACCCATATAGTTGAGGCTTCGAAGACTCTGGTTGACCTCCTTTTCGACCACATCATCCAACTCAGTCAGACGCACACCAGAGAAGACCTCAATGAGATTCTCAGGTTCAGTAAACGTCTTCCGATCGACCGATGGACCCGTATAAGTGCCTTCAGCAGTCATCCTCCCAAAAGCAGAGACCAGATTCGAGACAGCTGATTCGTGGGGATCCAAGCTGAGGAGTGCAAGAGAAAGTTGTTCTTCAGTACTGGGGGCAATCCACTTGGTCTGGTCGAAGATTGAGAATCGTGCGCCAAGCAGAACGAATCCCTTATACTCAGCAGCCTCGAAGATCATAAAGCATCTGCTCTGTGTATCAAGTGCTAAATCGATCCCCTTAAGTAGGTGAGAGAGCATCATCCCTTCCTTGGTCGTTGACAGAGATCCACAGCCGCGTCGAATAGTCAGAAAGGTGTCTTGAGACTTCTGCACCGTCGATCCGAGAAGTTGAAAGAACCGTTGAAGGATCACATTCGACATAAACACCGCGTCAGGCTGGACAAGACCGTGGAAGTATGGAAACGCAATTCCAGGAAGGTTGGGGACGTCGCGAGATGAGCCAATATTGGTCTGTGGCTTGCGACACGGTTTGGCCTTGTGGACAGTATCGTTTCCTGACGTGGAGAATCGTCCCTGAGGGGTTTCAAGTCGTGGATATGTCGTGGCGTCGTCGTCAGTCATCTGGTCAGGTGTGTACATCTCCCGTCCTTGAAGGGTTGACGAAAAGGCAATAGGAGCGAACGCGGATACATCCTTTCCTTCCTTGTCCTTCGTGTTGTAAGGGATATACTCAGTGTCCTTTACGTCATGAAGTGTGTGGAACTCTTGGATCGACTCGTCAGATTCGTATGAGAAGGCGGGGTATGAATGAGTACGCAAGAAGGTGGCGCAGACTGAGGACAAGAACTTGTACGCAAAGGCACCAGCTGGGTTGACGGCATCTCTACCACAGTTGAAATGGGCTGGAGTCAAGGCGATGATCGTATCCAACATGCCGTCAGGAAGTGCAGGTGCATGATAGATAGTCGTTTTACCAGAGGTACTCGTGGAGGCTCCAGGACCGAGACTGTTGACTAGCCTGATGGATTTGTCGGAATCGACCGCTTTCTTGAGCGTAAACTGCACGGGAGGACAGTTCAGAGGCAGAACATCTGATGAGTCAAAGAGTCCGAGGGGGAGTCCATTGACTAGAGTCGTTGCCTTGGTCTTTGTGGCCCAGAGGAAACAGGAGATGTTCGTCGCATAAGCAGTGTAGGCAGTCTGTTTGGTTAGGGCGACTCCGGCGAGATTGAGTGTTGAGATAGGCGAGGGGGCGAAGTTGTCCAATCCATGTAGTTTGAAGAGGGGAACGGTGTTCTGGCTGCCGGCAACGTGCCGAGCTGTGTAGATAGAGGAGGTCATGATGTCTGGGGACTGGAGGTGTGTAGTGTTGGCTCTTTTGTACCATCCCACTACACAGGCGAGATCTGTCGATCGACCTACCCACCAGCCCACTCAAGGACACCCCAGATCACTGCGTCCTGTGATACGGGCTTGTTGTGTGTCGATGTTCAAAAGAATGAAGCACACAGCCGCAACGTACTGGACCTGGTCTTGGACCTTCTCTCACCACACTGCTCTATGCACTGTCGTACAGACCATCCCATCCCTATCCTTACTACAAAGGAAAGACATACTCGGCCTCTGCAGGTCTCAAACACCGATCATACCAGGTTCCCACCAACACCGACAACTAAAGTGCTTAACCGACATCATCCAAGTTCCATCTGAACCTCAAAGAAGGGAAAACGACATCCAGTACTTGGAATACTCCTTCCATCGCTTTCTCAGAGACTACCAGGATGTATACCTTGACCATCTCTCCGCTGCTGTTATCGCTTTATCAGTCCGCGTCGCGCCCGGTATATCCGATTCAGACTTTCAGAACAGTCTCTCGTCCTCGATTTCCAAAGAGATTTCCGCTGCTATAACGATCTCCAACCAAGCACGCTTCGATTTCCGTTCGGTCATGGATTCGTCATTATCGTCCCGATTCATCCCTTTCACTGATAGACTGCCGCCGTCTCCAGAAGTCGCCTCTGATACGTATATCGCTCTTTGTCGAGTCCGCTCCCGACTTGATAATCGATTCCTACTGCCAACTCTTGGGAATAAGACCGTCGGCCAGGTTGTTCGTCATAACTTCTATGACATATATGACTCTTTCGGTTTAACCTTCAAACATGAGGGGAGGAGGAATGATGATACCGTGTCGACTTCGGATTGTATGCGTCTGTACCTCGAGACCGGAGTGTATCCACACGGACCTGTTGAAATGCGACGAGCCTGGACATATAATCAACTTGACCCAAGGGTGTATTACGCACGTGGAGGAGACGTTATGCACACATCCCAATACGTGCAGTCAATCGCCAATATGCTTATCGACGCTTTTCCCGAAACACACCGGAAAGATCGCTTCATGCCACCACGGGATCCACTAGCTGATGACGATGTTGAAGTCATATATGATTATTCATCGTTCACTTCGACTCTCGATTCGGTGGTTCCGTTTCTCGACAATCTCGCTGAGTTCTTCCGTGGTACTGTGGTTCATCTTGTCGACTTCAGGAACGGAGTTGTCCCAACCGACTTAGGTGATCTTATCGCACGATACAACACAGAGTGTAATCTCTATGGTACCTTCGACGCATCTAGAGTTCTGGGTCAATCTGCCGGTACGACTCTTCTTCAACATACCTGTGGAATGCTCGGAGTAGAAGGGAACATCTTCTTTGCCACACTCTTGCATGGGATACACCTCCGATTCATTGCTGGCCTGAACCGATCTCGTTGCGTTGGTGATGACGCTAGAATGCATCATAGGGTACCTTTCGGCATCATGGACAACACTGAAACCGATTACCTCGCTTGGGTTCTTGCTGGGTGTGGCGATCTGAGCAAGGAAAAGATGGGGAAATTCGAATCGGGTGTCGACAGTGAGCTCCAGGCCTATCGCTATATTAAGAGGCCGATACATCGCGACGGTTCCATTATGATTGAGGGGATACTCCTAACACTTCCTTCCATCATACCACTCCTCGGGGCCCAAGATCGATTCCACACTGTCACTCCATCCGTCAGTCATCCCTCGAGAAGAACGTACTCACAGATTCTCAGGTTCATCCAGGAGCTTTTCATTCATGGGTTATGTTACGACAGTGACGACGTTTCGTGGAAGTCGATATTGAAACATTTGATGTTCCTTAGGAGGTTATGTATTGCTGAGGATCCGGACTTTGAGCACTCCATGTTCATGAACTCATCTTATCTGACGAAGTATCGGTTTCCGCCACCAGAGACTTGGGGGAAGATGCCGATCACAGACTGGGTAGTAGGCGACATCATGTACGACGAGGTTATACGGTTTCCGATGAAGGGAGAGAGAGAATCAGAAGGAGGGTGTGACGGTCGCGTTGGTTCGGAAATGCTGAGGATGGGTTCCAAGGCTAGAGGATGTCTCGTGAAGCTAGGGTATCTGGAGGAGGAGAAGATGTTCGATGACGTTTCTGTCAAACTTGTCGGTTTGGACCTGTTCCTAGAGTACCTTGGAGGAAGATACCGTTCTATCAGCAAGTTCGTTGTCGTTAAGGACATACCTGGGTGGGTAGCACAAGTACCGAGTAGTCTATGAATGATGATAGAGCGAGATTTGCTATACCCCTAATGCCGGGGTACTGATGAAGCCATGGCAGGCCGAAAGCAAATCTCATATTATGGCATACGATACGTGTGCGTCAGTTCCCCACCTGACTCTCGTCGGGTGTGGCTGAATGCAACAGTATCTTATGCGCAGGCATGATGACTGTCTTTCGCGGCCATGTTCTTCCATCAGCGCTAGGGCATTAGGGAGTAGCGATTTCGTTCTACCATGAACTCATCAGGCACCCGGTATCAGGGCAC";
        let monomer = b"TATCGCACTACAGGCATGTCAATTAACGAACGACGAACTACGAAAAGATACGACATCGATCCTAAGAAAGTACACTATGGAACCAGGTCATCAACCAATCAATGTTGCCAGGAGGTCATCAAGATCAACCTCACCCTCACCAGCATCCTCTTCACGACGAGCTCTCTTCTTCTTTGGACCAGCTTCTCCAGCAATCTTGACCAACGACAGCCAAATCCGCTTCCTGAATTCCTCACTCTCCACCACCATATTCCGATATTTCCCAGCACGTTCATGGCTGTCTACCTTGATTGCCTGGCGTTCGAGCATCTTCTCAAAATCCTGGATGGCGATACTGAGAGGCTTCGGAGTAATGAAGAGTTTATCAGGCATGTTGGCAAAAGCATCCATCTCCCCAAGCTTACGTTTCCCACCCGAAACAACAGGTTTCTCTTTCGCTTCACACCGAATCTCGACACCTCGTTCATTCCAAAATGAGATCGCCTCAGGACCGAACCGAGAAAGGATGATAAACGCTCGACTCGAAAGGTTGACACGAATGGATGGAATGTACGTGGGACGACTGAGCTTAAGGACAGTATCAGAAGCAAGTGTCTCGAGCGCATCCGCAATCATCTGAGGGTTCCGAACGATGTACCCGGCACCCATATAGTTGAGGCTTCGAAGACTCTGGTTGACCTCCTTTTCGACCACATCATCCAACTCAGTCAGACGCACACCAGAGAAGACCTCAATGAGATTCTCAGGTTCAGTAAACGTCTTCCGATCGACCGATGGACCCGTATAAGTGCCTTCAGCAGTCATCCTCCCAAAAGCAGAGACCAGATTCGAGACAGCTGATTCGTGGGGATCCAAGCTGAGGAGTGCAAGAGAAAGTTGTTCTTCAGTACTGGGGGCAATCCACTTGGTCTGGTCGAAGATTGAGAATCGTGCGCCAAGCAGAACGAATCCCTTATACTCAGCAGCCTCGAAGATCATAAAGCATCTGCTCTGTGTATCAAGTGCTAAATCGATCCCCTTAAGTAGGTGAGAGAGCATCATCCCTTCCTTGGTCGTTGACAGAGATCCACAGCCGCGTCGAATAGTCAGAAAGGTGTCTTGAGACTTCTGCACCGTCGATCCGAGAAGTTGAAAGAACCGTTGAAGGATCACATTCGACATAAACACCGCGTCAGGCTGGACAAGACCGTGGAAGTATGGAAACGCAATTCCAGGAAGGTTGGGGACGTCGCGAGATGAGCCAATATTGGTCTGTGGCTTGCGACACGGTTTGGCCTTGTGGACAGTATCGTTTCCTGACGTGGAGAATCGTCCCTGAGGGGTTTCAAGTCGTGGATATGTCGTGGCGTCGTCGTCAGTCATCTGGTCAGGTGTGTACATCTCCCGTCCTTGAAGGGTTGACGAAAAGGCAATAGGAGCGAACGCGGATACATCCTTTCCTTCCTTGTCCTTCGTGTTGTAAGGGATATACTCAGTGTCCTTTACGTCATGAAGTGTGTGGAACTCTTGGATCGACTCGTCAGATTCGTATGAGAAGGCGGGGTATGAATGAGTACGCAAGAAGGTGGCGCAGACTGAGGACAAGAACTTGTACGCAAAGGCACCAGCTGGGTTGACGGCATCTCTACCACAGTTGAAATGGGCTGGAGTCAAGGCGATGATCGTATCCAACATGCCGTCAGGAAGTGCAGGTGCATGATAGATAGTCGTTTTACCAGAGGTACTCGTGGAGGCTCCAGGACCGAGACTGTTGACTAGCCTGATGGATTTGTCGGAATCGACCGCTTTCTTGAGCGTAAACTGCACGGGAGGACAGTTCAGAGGCAGAACATCTGATGAGTCAAAGAGTCCGAGGGGGAGTCCATTGACTAGAGTCGTTGCCTTGGTCTTTGTGGCCCAGAGGAAACAGGAGATGTTCGTCGCATAAGCAGTGTAGGCAGTCTGTTTGGTTAGGGCGACTCCGGCGAGATTGAGTGTTGAGATAGGCGAGGGGGCGAAGTTGTCCAATCCATGTAGTTTGAAGAGGGGAACGGTGTTCTGGCTGCCGGCAACGTGCCGAGCTGTGTAGATAGAGGAGGTCATGATGTCTGGGGACTGGAGGTGTGTAGTGTTGGCTCTTTTGTACCATCCCACTACACAGGCGAGATCTGTCGATCGACCTACCCACCAGCCCACTCAAGGACACCCCAGATCACTGCGTCCTGTGATACGGGCTTGTTGTGTGTCGATGTTCAAAAGAATGAAGCACACAGCCGCAACGTACTGGACCTGGTCTTGGACCTTCTCTCACCACACTGCTCTATGCACTGTCGTACAGACCATCCCATCCCTATCCTTACTACAAAGGAAAGACATACTCGGCCTCTGCAGGTCTCAAACACCGATCATACCAGGTTCCCACCAACACCGACAACTAAAGTGCTTAACCGACATCATCCAAGTTCCATCTGAACCTCAAAGAAGGGAAAACGACATCCAGTACTTGGAATACTCCTTCCATCGCTTTCTCAGAGACTACCAGGATGTATACCTTGACCATCTCTCCGCTGCTGTTATCGCTTTATCAGTCCGCGTCGCGCCCGGTATATCCGATTCAGACTTTCAGAACAGTCTCTCGTCCTCGATTTCCAAAGAGATTTCCGCTGCTATAACGATCTCCAACCAAGCACGCTTCGATTTCCGTTCGGTCATGGATTCGTCATTATCGTCCCGATTCATCCCTTTCACTGATAGACTGCCGCCGTCTCCAGAAGTCGCCTCTGATACGTATATCGCTCTTTGTCGAGTCCGCTCCCGACTTGATAATCGATTCCTACTGCCAACTCTTGGGAATAAGACCGTCGGCCAGGTTGTTCGTCATAACTTCTATGACATATATGACTCTTTCGGTTTAACCTTCAAACATGAGGGGAGGAGGAATGATGATACCGTGTCGACTTCGGATTGTATGCGTCTGTACCTCGAGACCGGAGTGTATCCACACGGACCTGTTGAAATGCGACGAGCCTGGACATATAATCAACTTGACCCAAGGGTGTATTACGCACGTGGAGGAGACGTTATGCACACATCCCAATACGTGCAGTCAATCGCCAATATGCTTATCGACGCTTTTCCCGAAACACACCGGAAAGATCGCTTCATGCCACCACGGGATCCACTAGCTGATGACGATGTTGAAGTCATATATGATTATTCATCGTTCACTTCGACTCTCGATTCGGTGGTTCCGTTTCTCGACAATCTCGCTGAGTTCTTCCGTGGTACTGTGGTTCATCTTGTCGACTTCAGGAACGGAGTTGTCCCAACCGACTTAGGTGATCTTATCGCACGATACAACACAGAGTGTAATCTCTATGGTACCTTCGACGCATCTAGAGTTCTGGGTCAATCTGCCGGTACGACTCTTCTTCAACATACCTGTGGAATGCTCGGAGTAGAAGGGAACATCTTCTTTGCCACACTCTTGCATGGGATACACCTCCGATTCATTGCTGGCCTGAACCGATCTCGTTGCGTTGGTGATGACGCTAGAATGCATCATAGGGTACCTTTCGGCATCATGGACAACACTGAAACCGATTACCTCGCTTGGGTTCTTGCTGGGTGTGGCGATCTGAGCAAGGAAAAGATGGGGAAATTCGAATCGGGTGTCGACAGTGAGCTCCAGGCCTATCGCTATATTAAGAGGCCGATACATCGCGACGGTTCCATTATGATTGAGGGGATACTCCTAACACTTCCTTCCATCATACCACTCCTCGGGGCCCAAGATCGATTCCACACTGTCACTCCATCCGTCAGTCATCCCTCGAGAAGAACGTACTCACAGATTCTCAGGTTCATCCAGGAGCTTTTCATTCATGGGTTATGTTACGACAGTGACGACGTTTCGTGGAAGTCGATATTGAAACATTTGATGTTCCTTAGGAGGTTATGTATTGCTGAGGATCCGGACTTTGAGCACTCCATGTTCATGAACTCATCTTATCTGACGAAGTATCGGTTTCCGCCACCAGAGACTTGGGGGAAGATGCCGATCACAGACTGGGTAGTAGGCGACATCATGTACGACGAGGTTATACGGTTTCCGATGAAGGGAGAGAGAGAATCAGAAGGAGGGTGTGACGGTCGCGTTGGTTCGGAAATGCTGAGGATGGGTTCCAAGGCTAGAGGATGTCTCGTGAAGCTAGGGTATCTGGAGGAGGAGAAGATGTTCGATGACGTTTCTGTCAAACTTGTCGGTTTGGACCTGTTCCTAGAGTACCTTGGAGGAAGATACCGTTCTATCAGCAAGTTCGTTGTCGTTAAGGACATACCTGGGTGGGTAGCACAAGTACCGAGTAGTCTATGAATGATGATAGAGCGAGATTTGCTATACCCCTAATGCCGGGGTACTGATGAAGCCATGGCAGGCCGAAAGCAAATCTCATATTATGGCATACGATACGTGTGCGTCAGTTCCCCACCTGACTCTCGTCGGGTGTGGCTGAATGCAACAGTATCTTATGCGCAGGCATGATGACTGTCTTTCGCGGCCATGTTCTTCCATCAGCGCTAGGGCATTAGGGAGTAGCGATTTCGTTCTACCATGAACTCATCAGGCACCCGGTATCAGGGCAC";
        assert_eq!(
            Monomerizer::builder()
                .seed_len(20)
                .overlap_dist(0)
                .build()
                .unwrap()
                .monomerize(input),
            monomer
        );
    }

    mod overlap_percent {
        use super::*;
        use pretty_assertions::assert_eq;

        #[test]
        fn dimer_with_overlap_percentage() {
            let input = b"ATGCCCATGCGCCAGCGCAGATGCGAATGCGCCAGCGCAG";
            let output = b"ATGCGAATGCGCCAGCGCAG";

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
            let output = b"ATGCGAATGCGCCAGCGCAG";

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
}
