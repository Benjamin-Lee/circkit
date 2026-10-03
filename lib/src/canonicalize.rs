use bio::alphabets;

/// Compute the lexicographically minimal string rotation of a string.
/// https://cp-algorithms.com/string/lyndon_factorization.html#finding-the-smallest-cyclic-shift
fn lmsr(s: &[u8]) -> Vec<u8> {
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
    doubled[ans..ans + n].to_vec()
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
        for i in n..s.len() {
            res.push(s[i]);
        }
        for i in 0..n {
            res.push(s[i]);
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
