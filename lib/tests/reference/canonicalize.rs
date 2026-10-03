use bio::alphabets;

/// Return the zero-based byte offset of the lexicographically minimal circular rotation.
///
/// Bytes are compared directly, without normalizing the alphabet or changing strands.
/// Empty input returns zero. If multiple offsets produce the same minimum rotation,
/// the smallest offset is returned.
///
/// Uses Duval's Lyndon factorization algorithm in linear time.
/// https://cp-algorithms.com/string/lyndon_factorization.html#finding-the-smallest-cyclic-shift
///
/// ```
/// use circkit::canonicalize::lmsr_index;
///
/// assert_eq!(lmsr_index(b"TAA"), 1);
/// ```
pub fn lmsr_index(s: &[u8]) -> usize {
    let n = s.len();
    let doubled: Vec<u8> = s.iter().chain(s.iter()).copied().collect();
    let mut i = 0;
    let mut ans = 0;
    while i < n {
        ans = i;
        let mut j = i + 1;
        let mut k = i;
        while j < 2 * n && doubled[k] <= doubled[j] {
            if doubled[k] < doubled[j] {
                k = i;
            } else {
                k += 1;
            }
            j += 1;
        }

        while i <= k {
            i += j - k;
        }
    }
    ans
}

/// Return the lexicographically minimal circular rotation, preserving the supplied strand.
///
/// The rotation starts at the offset returned by [`lmsr_index`]. Empty input returns
/// an empty vector. Bytes are compared directly, without alphabet normalization.
///
/// ```
/// use circkit::canonicalize::lmsr;
///
/// assert_eq!(lmsr(b"TAA"), b"AAT");
/// ```
pub fn lmsr(s: &[u8]) -> Vec<u8> {
    let index = lmsr_index(s);
    let mut rotated = Vec::with_capacity(s.len());
    rotated.extend_from_slice(&s[index..]);
    rotated.extend_from_slice(&s[..index]);
    rotated
}

/// Canonicalize a circular DNA sequence.
///
/// This function computes the lexicographically minimal string rotation of a string and its reverse complement, and returns the smaller of the two.
/// It does not check that the input is a valid DNA sequence, so RNA sequences will have unexpected results.
/// Ensure that the input is a valid DNA sequence before calling this function.
/// Non-ATGC characters will be treated normally, meaning that they too will be used when sorting lexicographically.
pub fn canonicalize(s: &[u8]) -> Vec<u8> {
    let lmsr_s = lmsr(s);
    let lmsr_revcomp_s = lmsr(&alphabets::dna::revcomp(&lmsr_s));

    if lmsr_s < lmsr_revcomp_s {
        lmsr_s
    } else {
        lmsr_revcomp_s
    }
}
