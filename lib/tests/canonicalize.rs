use circkit::canonicalize::{lmsr, lmsr_index};
use proptest::prelude::*;
use rstest::rstest;

#[rstest]
#[case(b"", 0, b"")]
#[case(b"A", 0, b"A")]
#[case(b"AAAA", 0, b"AAAA")]
#[case(b"TAA", 1, b"AAT")]
#[case(b"GATGAT", 1, b"ATGATG")]
#[case(b"ACACAC", 0, b"ACACAC")]
#[case(b"CACACA", 1, b"ACACAC")]
#[case(b"banana", 5, b"abanan")]
#[case(b"\xff\0\xff\0", 1, b"\0\xff\0\xff")]
fn public_rotation_api(
    #[case] sequence: &[u8],
    #[case] expected_index: usize,
    #[case] expected_rotation: &[u8],
) {
    assert_eq!(lmsr_index(sequence), expected_index);
    assert_eq!(lmsr(sequence), expected_rotation);
}

proptest! {
    #[test]
    fn public_rotation_matches_brute_force(sequence in prop::collection::vec(any::<u8>(), 0..128)) {
        let (expected_index, expected_rotation) = (0..sequence.len())
            .map(|index| {
                let rotation = sequence[index..]
                    .iter()
                    .chain(&sequence[..index])
                    .copied()
                    .collect::<Vec<_>>();
                (index, rotation)
            })
            .min_by(|(left_index, left), (right_index, right)| {
                left.cmp(right).then_with(|| left_index.cmp(right_index))
            })
            .unwrap_or_default();

        prop_assert_eq!(lmsr_index(&sequence), expected_index);
        prop_assert_eq!(lmsr(&sequence), expected_rotation);
    }
}
