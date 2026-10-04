use super::Settings;
use circkit::{canonicalize::RotationOptions, Monomerizer};
use std::{borrow::Cow, num::NonZeroUsize};

pub fn index(sequence: &[u8], settings: Settings) -> usize {
    RotationOptions {
        duval_max_len: settings.rotation_cutoff,
    }
    .lmsr_index(sequence)
}

pub fn canonicalize(sequence: &[u8], settings: Settings) -> Vec<u8> {
    RotationOptions {
        duval_max_len: settings.rotation_cutoff,
    }
    .canonicalize(sequence)
}

pub fn monomerizer(settings: Settings) -> Monomerizer {
    Monomerizer::builder()
        .seed_len(10)
        .overlap_dist(settings.max_mismatch)
        .mismatch_chunk_size(NonZeroUsize::new(settings.chunk_size).unwrap())
        .build()
        .unwrap()
}

pub fn normalize(sequence: &[u8]) -> Cow<'_, [u8]> {
    circkit_cli::utils::normalized_sequence(sequence)
}
