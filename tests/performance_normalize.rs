use circkit_cli::utils::normalized_sequence;
use std::borrow::Cow;

#[test]
fn normalization_preserves_all_byte_mappings_and_whitespace() {
    for byte in 0..=255_u8 {
        for length in [1, 15, 31, 32, 33, 63, 64, 150] {
            let mut input = vec![b'A'; length];
            for index in 0..length {
                input[index] = byte;
                let expected =
                    needletail::sequence::normalize(&input, false).unwrap_or_else(|| input.clone());
                assert_eq!(normalized_sequence(&input).as_ref(), expected);
                input[index] = b'A';
            }
        }
    }
    assert!(matches!(normalized_sequence(b"ACGTN-"), Cow::Borrowed(_)));
}
