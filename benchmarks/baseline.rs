// Adapter compiled against the PR #2 snapshot by run.py. It uses only that API.
use super::Settings;
use circkit::Monomerizer;
use std::borrow::Cow;

pub fn index(sequence: &[u8], _: Settings) -> usize {
    circkit::canonicalize::lmsr_index(sequence)
}

pub fn canonicalize(sequence: &[u8], _: Settings) -> Vec<u8> {
    circkit::canonicalize(sequence)
}

pub fn monomerizer(settings: Settings) -> Monomerizer {
    Monomerizer::builder()
        .seed_len(10)
        .overlap_dist(settings.max_mismatch)
        .build()
        .unwrap()
}

pub fn normalize(sequence: &[u8]) -> Cow<'_, [u8]> {
    Cow::Owned(
        needletail::sequence::normalize(sequence, false).unwrap_or_else(|| sequence.to_vec()),
    )
}
