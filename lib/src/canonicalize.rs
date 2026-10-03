use bio::alphabets;

/// Performance choices that do not change the rotation or its tie-breaking.
#[derive(Clone, Copy, Debug)]
pub struct RotationOptions {
    /// Use contiguous Duval through this length, and a constant-space search above it.
    /// The default (32,768 bytes) is empirical. Zero selects constant-space search
    /// for every nonempty input; `usize::MAX` selects Duval for every input.
    pub duval_max_len: usize,
}

impl Default for RotationOptions {
    fn default() -> Self {
        Self {
            duval_max_len: 32768,
        }
    }
}

/// Return the zero-based byte offset of the lexicographically minimal circular rotation.
///
/// Bytes are compared directly, without normalizing the alphabet or changing strands.
/// Empty input returns zero. If multiple offsets produce the same minimum rotation,
/// the smallest offset is returned.
///
/// Uses linear-time searches: contiguous Duval factorization for short inputs,
/// and a two-candidate search with constant auxiliary space for longer inputs.
///
/// ```
/// use circkit::canonicalize::lmsr_index;
///
/// assert_eq!(lmsr_index(b"TAA"), 1);
/// ```
pub fn lmsr_index(s: &[u8]) -> usize {
    RotationOptions::default().lmsr_index(s)
}

fn duval_allocating_index(s: &[u8]) -> usize {
    let n = s.len();
    let doubled: Vec<u8> = s.iter().chain(s.iter()).copied().collect();
    let (mut i, mut answer) = (0, 0);
    while i < n {
        answer = i;
        let (mut j, mut k) = (i + 1, i);
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
    answer
}

fn duval_index(s: &[u8], doubled: &mut Vec<u8>) -> usize {
    let n = s.len();
    doubled.clear();
    doubled.reserve(2 * n);
    doubled.extend_from_slice(s);
    doubled.extend_from_slice(s);
    let (mut i, mut answer) = (0, 0);
    while i < n {
        answer = i;
        let (mut j, mut k) = (i + 1, i);
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
    answer
}

fn two_candidate_index(s: &[u8]) -> usize {
    let n = s.len();
    if n < 2 {
        return 0;
    }
    let (mut i, mut j, mut matched) = (0, 1, 0);
    while i < n && j < n && matched < n {
        let left = i + matched;
        let right = j + matched;
        let a = s[if left >= n { left - n } else { left }];
        let b = s[if right >= n { right - n } else { right }];
        match a.cmp(&b) {
            std::cmp::Ordering::Equal => matched += 1,
            std::cmp::Ordering::Greater => {
                // Every candidate in i..=i+matched loses to j at the mismatch.
                i += matched + 1;
                if i <= j {
                    i = j + 1;
                }
                matched = 0;
            }
            std::cmp::Ordering::Less => {
                j += matched + 1;
                if j <= i {
                    j = i + 1;
                }
                matched = 0;
            }
        }
    }
    i.min(j)
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
    RotationOptions::default().lmsr(s)
}

/// Canonicalize a circular DNA sequence.
///
/// This function computes the lexicographically minimal string rotation of a string and its reverse complement, and returns the smaller of the two.
/// It does not check that the input is a valid DNA sequence, so RNA sequences will have unexpected results.
/// Ensure that the input is a valid DNA sequence before calling this function.
/// Non-ATGC characters will be treated normally, meaning that they too will be used when sorting lexicographically.
pub fn canonicalize(s: &[u8]) -> Vec<u8> {
    RotationOptions::default().canonicalize(s)
}

/// Canonicalize into reusable buffers, avoiding per-record allocations in pipelines.
pub fn canonicalize_into(s: &[u8], output: &mut Vec<u8>, reverse: &mut Vec<u8>) {
    RotationOptions::default().canonicalize_into(s, output, reverse);
}

impl RotationOptions {
    /// Compute [`lmsr_index`] with a configurable algorithm cutoff.
    pub fn lmsr_index(self, s: &[u8]) -> usize {
        if s.len() <= self.duval_max_len {
            duval_allocating_index(s)
        } else {
            two_candidate_index(s)
        }
    }

    /// Compute [`lmsr`] with a configurable algorithm cutoff.
    pub fn lmsr(self, s: &[u8]) -> Vec<u8> {
        let index = self.lmsr_index(s);
        let mut rotated = Vec::with_capacity(s.len());
        rotated.extend_from_slice(&s[index..]);
        rotated.extend_from_slice(&s[..index]);
        rotated
    }

    /// Compute [`canonicalize`] with a configurable algorithm cutoff.
    pub fn canonicalize(self, s: &[u8]) -> Vec<u8> {
        if s.len() <= self.duval_max_len {
            let forward = self.lmsr(s);
            let backward = self.lmsr(&alphabets::dna::revcomp(&forward));
            return if forward < backward {
                forward
            } else {
                backward
            };
        }
        let mut output = Vec::with_capacity(s.len());
        let mut reverse = Vec::with_capacity(s.len());
        self.canonicalize_into(s, &mut output, &mut reverse);
        output
    }

    /// Compute [`canonicalize_into`] with a configurable algorithm cutoff.
    pub fn canonicalize_into(self, s: &[u8], output: &mut Vec<u8>, reverse: &mut Vec<u8>) {
        reverse.clear();
        reverse.extend(s.iter().rev().map(|&base| alphabets::dna::complement(base)));
        // The output buffer doubles as scratch until the winning rotation is known.
        let (forward_index, reverse_index) = if s.len() <= self.duval_max_len {
            (duval_index(s, output), duval_index(reverse, output))
        } else {
            (two_candidate_index(s), two_candidate_index(reverse))
        };
        let (sequence, index) =
            if compare_rotations(s, forward_index, reverse, reverse_index).is_lt() {
                (s, forward_index)
            } else {
                (reverse.as_slice(), reverse_index)
            };
        output.clear();
        output.extend_from_slice(&sequence[index..]);
        output.extend_from_slice(&sequence[..index]);
    }
}

#[inline]
fn compare_rotations(a: &[u8], mut left: usize, b: &[u8], mut right: usize) -> std::cmp::Ordering {
    let mut remaining = a.len();
    while remaining > 0 {
        // At most three contiguous comparisons cover two equal-length rotations.
        // Slice comparison can use the platform's optimized byte comparison.
        let length = remaining.min(a.len() - left).min(b.len() - right);
        let order = a[left..left + length].cmp(&b[right..right + length]);
        if !order.is_eq() {
            return order;
        }
        remaining -= length;
        left += length;
        right += length;
        if left == a.len() {
            left = 0;
        }
        if right == b.len() {
            right = 0;
        }
    }
    std::cmp::Ordering::Equal
}

#[cfg(test)]
mod lmsr_index_test {
    use super::*;

    #[test]
    fn aaa() {
        assert_eq!(lmsr_index(b"AAA"), 0);
    }

    #[test]
    fn banana() {
        assert_eq!(lmsr_index(b"banana"), 5);
    }

    #[test]
    fn taa() {
        assert_eq!(lmsr_index(b"TAA"), 1);
    }
}

#[cfg(test)]
mod lmsr_test {
    use super::*;
    #[test]
    fn aaa() {
        assert_eq!(lmsr(b"AAA"), b"AAA");
    }
    #[test]
    fn banana() {
        assert_eq!(lmsr(b"banana"), b"abanan");
    }
    #[test]
    fn taa() {
        assert_eq!(lmsr(b"TAA"), b"AAT");
    }

    #[test]
    fn second_application_is_identical() {
        let tmp = lmsr(b"ATGCAGATACAGA");
        let tmp2 = lmsr(&tmp);
        assert_eq!(tmp, tmp2);
    }
}

#[cfg(test)]
mod canonicalize_test {
    use super::*;
    #[test]
    fn aaa() {
        assert_eq!(canonicalize(b"AAA"), b"AAA");
    }
    #[test]
    fn att() {
        assert_eq!(canonicalize(b"ATT"), b"AAT");
    }

    #[test]
    fn real_monomer() {
        // Drawn from 3300000336_thermBogB3DRAFT_128220 in cated_Soil_microbial_communities_from_permafrost_in_Bonanza_Creek__Alaska
        // AATCAATTTCCTCCATCACCTAGTTTATGTAGAAACGCTGCTA
        //         |||||||||||||||||||||||||||||||||||
        //         TCCTCCATCACCTAGTTTATGTAGAAACGCTGCTAAATCAATT

        let a = "AATCAATTTCCTCCATCACCTAGTTTATGTAGAAACGCTGCTA";
        let b = "TCCTCCATCACCTAGTTTATGTAGAAACGCTGCTAAATCAATT";
        assert_eq!(lmsr(a.as_bytes()), lmsr(b.as_bytes()));
        assert_eq!(canonicalize(a.as_bytes()), canonicalize(b.as_bytes()));
    }
}

#[cfg(test)]
/// We have multiple implementations of lmsr_index, so we can compare them against each other to make sure the optimized version is correct
mod fuzzing {
    use super::*;
    use proptest::prelude::*;

    /// An auxiliary function to rotate a string by n characters
    fn rotate(s: &str, n: usize) -> String {
        let s = s.chars().collect::<Vec<char>>();
        let mut res = String::new();
        for &character in &s[n..] {
            res.push(character);
        }
        for &character in &s[..n] {
            res.push(character);
        }
        res
    }

    fn lmsr_index_simple(s: &str) -> usize {
        let mut result = 0;

        for i in 0..s.len() {
            let rotated = rotate(s, i);
            if rotated < rotate(s, result) {
                result = i;
            }
        }
        result
    }

    /// Find starting position of minimum acyclic string in (s)
    /// https://codeforces.com/blog/entry/90035#duval
    fn lmsr_index_2(s: &str) -> usize {
        let n = s.len(); // the real size of the string
        let mut s = s.chars().collect::<Vec<char>>(); // convert string to char vector

        let s_extended = s.clone(); // Clone s to avoid borrowing conflicts
        s.extend(&s_extended); // for convention since we are dealing with acyclic

        let mut res = 0; // minimum acyclic string

        // while s2 is a lyndon word, try to add s2 with s[p]
        let mut l = 0;
        while l < n {
            res = l;

            // Extend as much as possible lyndon word s2 = s[l..r]
            let mut r = l;
            let mut p = l + 1;
            while p < s.len() {
                // (s2 + s[p]) is not a lyndon word
                if s[r] > s[p] {
                    break;
                }

                // (s2 + s[p]) is still a lyndon word, hence extend s2
                if s[r] == s[p] {
                    r += 1;
                    p += 1;
                    continue;
                }

                // (s2 + s[p]) is a lyndon word, but it may be a repeated string
                if s[r] < s[p] {
                    r = l;
                    p += 1;
                    continue;
                }
            }

            // The lyndon word may have the form of s2 = sx + sx + .. + sx like "12312123"
            while l <= r {
                l += p - r;
            }
        }

        // Don't forget to return the value ;)
        res
    }

    proptest! {
        #[test]
        fn lmsr_index_implementations_are_identical(s in "[ -~]{1, 100}") {
            prop_assert_eq!(lmsr_index_2(&s), lmsr_index_simple(&s));
            prop_assert_eq!(lmsr_index_2(&s), lmsr_index(s.as_bytes()));
        }

        #[test]
        fn lmsr_is_idempotent(s in "[ -~]{1, 100}") {
            prop_assert_eq!(lmsr(&lmsr(s.as_bytes())), lmsr(s.as_bytes()));
        }
        #[test]
        fn canonicalize_is_idempotent(s in "[ATGC]{1, 100}") {
            prop_assert_eq!(canonicalize(&canonicalize(s.as_bytes())), canonicalize(s.as_bytes()));
        }
    }
}
