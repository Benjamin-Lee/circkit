#[macro_use]
extern crate derive_builder;

#[allow(dead_code)]
#[path = "reference/canonicalize.rs"]
mod baseline_canonicalize;
#[allow(dead_code)]
#[path = "reference/monomerize.rs"]
mod baseline_monomerize;
#[allow(dead_code)]
#[path = "reference/orfs.rs"]
mod baseline_orfs;

use proptest::prelude::*;

proptest! {
    #![proptest_config(ProptestConfig::with_cases(512))]

    #[test]
    fn circular_orfs_match_baseline(sequence in "[ACGTNacgt-]{3,1000}", config in 0..5usize) {
        let starts: &[&str] = match config {
            0 => &["ATG"], 1 => &["ATG", "CTG", "TTG"],
            2 => &["NNN", "ATG"], 3 => &["ATG", "ATG"], _ => &["TAA", "ATG"],
        };
        let stops = ["TAA", "TAG", "TGA"];
        let old = baseline_orfs::start_stop_codon_indices_by_frame_naive(&sequence, starts, &stops);
        let matcher = circkit::orfs::CodonMatcher::new(starts, &stops);
        let new = matcher.indices(sequence.as_bytes());
        prop_assert_eq!(&old, &new);
        let old_orfs = baseline_orfs::find_orfs_with_indices(sequence.len(), old.0, old.1);
        let new_orfs = circkit::orfs::find_orfs_with_indices(sequence.len(), new.0, new.1);
        let old_fields: Vec<_> = old_orfs.iter().map(|o| (o.start,o.stop,o.wraps,o.length)).collect();
        let new_fields: Vec<_> = new_orfs.iter().map(|o| (o.start,o.stop,o.wraps,o.length)).collect();
        prop_assert_eq!(old_fields, new_fields);
    }

    #[test]
    fn long_and_short_rotations_preserve_byte_order(sequence in prop::collection::vec(any::<u8>(),0..80000), cutoff in prop_oneof![Just(0), Just(usize::MAX), 1..80000usize]) {
        let options = circkit::canonicalize::RotationOptions { duval_max_len: cutoff };
        prop_assert_eq!(baseline_canonicalize::lmsr_index(&sequence),options.lmsr_index(&sequence));
        prop_assert_eq!(baseline_canonicalize::canonicalize(&sequence),options.canonicalize(&sequence));
        let (mut output, mut reverse) = (Vec::new(), Vec::new());
        options.canonicalize_into(&sequence, &mut output, &mut reverse);
        prop_assert_eq!(baseline_canonicalize::canonicalize(&sequence), output);
    }

    #[test]
    fn arbitrary_index_order_preserves_behavior(
        length in 3..80usize,
        starts in prop::collection::vec(prop::collection::vec(0..80usize,0..20),3),
        stops in prop::collection::vec(prop::collection::vec(0..80usize,0..20),3),
    ) {
        let starts: Vec<Vec<_>> = starts.into_iter().map(|v| v.into_iter().filter(|&x| x < length).collect()).collect();
        let stops: Vec<Vec<_>> = stops.into_iter().map(|v| v.into_iter().filter(|&x| x < length).collect()).collect();
        let old = baseline_orfs::find_orfs_with_indices(length, starts.clone(), stops.clone());
        let new = circkit::orfs::find_orfs_with_indices(length, starts, stops);
        prop_assert_eq!(old.iter().map(|o| (o.start,o.stop,o.wraps,o.length)).collect::<Vec<_>>(),
                        new.iter().map(|o| (o.start,o.stop,o.wraps,o.length)).collect::<Vec<_>>());
    }

    #[test]
    fn monomerization_matches_baseline(sequence in "[ACGTN-]{1,500}", seed in 1..30usize, config in 0..8usize, chunk in 1..1024usize) {
        let mut old = baseline_monomerize::Monomerizer::builder();
        let mut new = circkit::Monomerizer::builder();
        old.seed_len(seed); new.seed_len(seed);
        new.mismatch_chunk_size(std::num::NonZeroUsize::new(chunk).unwrap());
        if config < 4 {
            old.overlap_dist(config as u64); new.overlap_dist(config as u64);
        } else {
            let identity = [0.0,0.5,0.95,1.0][config-4];
            old.overlap_min_identity(identity); new.overlap_min_identity(identity);
        }
        let old = old.build().unwrap(); let new = new.build().unwrap();
        let bytes = sequence.as_bytes();
        prop_assert_eq!(old.first_monomer_end_index(bytes), new.first_monomer_end_index(bytes));
        prop_assert_eq!(old.last_monomer_end_index(bytes), new.last_monomer_end_index(bytes));
        prop_assert_eq!(old.last_monomer_end_index_sensitive(bytes), new.last_monomer_end_index_sensitive(bytes));
    }

    #[test]
    fn circular_extraction_matches_baseline(sequence in "[ACGT]{1,100}", start in 0..400usize, length in 3..800usize, include_stop in any::<bool>()) {
        let old = baseline_orfs::Orf {start,stop:None,wraps:3,length};
        let new = circkit::orfs::Orf {start,stop:None,wraps:3,length};
        prop_assert_eq!(old.seq_with_opts(sequence.as_bytes(),include_stop),new.seq_with_opts(sequence.as_bytes(),include_stop));
        let mut written = Vec::new();
        new.write_seq_with_opts(sequence.as_bytes(),include_stop,&mut written).unwrap();
        let materialized = new.seq_with_opts(sequence.as_bytes(),include_stop);
        prop_assert_eq!(written,materialized.as_bytes());
    }
}

#[test]
fn long_periodic_rotations_preserve_earliest_ties_and_reusable_buffers() {
    let mut output = Vec::new();
    let mut reverse = Vec::new();
    for length in [0, 1, 4096, 4097, 8192, 32768, 32769, 65536] {
        for motif in [b"A".as_slice(), b"ACGT", b"AAACAAA", &[255, 0, 127]] {
            let sequence: Vec<_> = motif.iter().copied().cycle().take(length).collect();
            assert_eq!(
                baseline_canonicalize::lmsr_index(&sequence),
                circkit::canonicalize::lmsr_index(&sequence)
            );
            circkit::canonicalize::canonicalize_into(&sequence, &mut output, &mut reverse);
            assert_eq!(output, baseline_canonicalize::canonicalize(&sequence));
        }
    }
}

#[test]
fn repetitive_monomerization_preserves_overlapping_seeds() {
    for length in [11, 31, 100, 1000] {
        for seed in [1, 5, 10, 30, 63] {
            for motif in [b"A".as_slice(), b"ACGT".as_slice(), b"AAACAAA".as_slice()] {
                let sequence: Vec<u8> = motif.iter().copied().cycle().take(length).collect();
                let old = baseline_monomerize::Monomerizer::builder()
                    .seed_len(seed)
                    .build()
                    .unwrap();
                let new = circkit::Monomerizer::builder()
                    .seed_len(seed)
                    .build()
                    .unwrap();
                assert_eq!(
                    old.last_monomer_end_index(&sequence),
                    new.last_monomer_end_index(&sequence)
                );
            }
        }
    }
}

#[test]
fn longest_orf_output_and_mutation_order_are_preserved() {
    let mut old: Vec<_> = (0..100)
        .map(|i| baseline_orfs::Orf {
            start: i,
            stop: if i % 5 == 0 { None } else { Some(i % 11) },
            length: 3 * (i % 7 + 1),
            wraps: i % 4,
        })
        .collect();
    let mut new: Vec<_> = old
        .iter()
        .map(|o| circkit::orfs::Orf {
            start: o.start,
            stop: o.stop,
            length: o.length,
            wraps: o.wraps,
        })
        .collect();
    let a = baseline_orfs::longest_orfs(&mut old);
    let b = circkit::orfs::longest_orfs(&mut new);
    assert_eq!(
        a.iter().map(|o| o.start).collect::<Vec<_>>(),
        b.iter().map(|o| o.start).collect::<Vec<_>>()
    );
    assert_eq!(
        old.iter().map(|o| o.start).collect::<Vec<_>>(),
        new.iter().map(|o| o.start).collect::<Vec<_>>()
    );
}
