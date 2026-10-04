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
