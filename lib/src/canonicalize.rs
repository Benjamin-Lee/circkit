use bio::alphabets;

/// Compute index of the lexicographically minimal string rotation of a string.
pub fn lmsr_index(s: &[u8]) -> usize {
    let n = s.len();
    let mut i = 0;
    let mut ans = 0;
    while i < n {
        ans = i;
        let mut j = i + 1;
        let mut k = i;
        while j < n && s[k % n] <= s[j % n] {
            if s[k % n] < s[j % n] {
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

/// Compute the lexicographically minimal string rotation of a string.
///
/// Internally, this function computes the index of the LMSR and then converts it to a `Vec<u8>`.
pub fn lmsr(s: &[u8]) -> Vec<u8> {
    // TODO: compensate for RNA/DNA ambiguity
    let mut buf = Vec::<u8>::with_capacity(s.len());
    let i = lmsr_index(s);
    buf.extend_from_slice(&s[i..]);
    buf.extend_from_slice(&s[..i]);
    buf
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
    fn auu() {
        assert_eq!(canonicalize(b"AUU"), b"AAU");
    }
}
