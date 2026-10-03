//! Shared benchmark driver. run.py changes only the adapter path in the baseline snapshot.
#[path = "../benchmarks/optimized.rs"]
mod implementation;
use std::{fmt::Debug, hint::black_box, time::Instant};

#[derive(Clone, Copy)]
pub struct Settings {
    pub rotation_cutoff: usize,
    pub chunk_size: usize,
    pub max_mismatch: u64,
}

#[inline(always)]
fn ignore<T: Debug>(value: T) -> u64 {
    black_box(value);
    0
}

fn digest<T: Debug>(value: T) -> u64 {
    format!("{value:?}")
        .bytes()
        .fold(0xcbf29ce484222325, |hash, byte| {
            (hash ^ u64::from(byte)).wrapping_mul(0x100000001b3)
        })
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 4 {
        eprintln!("Usage: performance OPERATION LENGTH ITERATIONS [PATTERN [ROTATION_CUTOFF [CHUNK_SIZE [MAX_MISMATCH]]]]");
        std::process::exit(2);
    }
    let operation = args[1].as_str();
    let length: usize = args[2].parse().unwrap();
    let iterations: usize = args[3].parse().unwrap();
    assert!(length > 0 && iterations > 0);
    let pattern = args.get(4).map(String::as_str).unwrap_or("random");
    let settings = Settings {
        rotation_cutoff: args.get(5).map_or(32768, |s| s.parse().unwrap()),
        chunk_size: args.get(6).map_or(256, |s| s.parse().unwrap()),
        max_mismatch: args.get(7).map_or(0, |s| s.parse().unwrap()),
    };
    let mut state = 0x123456789abcdef_u64;
    let sequences: Vec<Vec<u8>> = (0..8)
        .map(|_| {
            let mut sequence: Vec<u8> = (0..length)
                .map(|i| {
                    state ^= state << 13;
                    state ^= state >> 7;
                    state ^= state << 17;
                    match pattern {
                        "dense" => b"ATGAAATAA"[i % 9],
                        "homopolymer" => b'A',
                        "false-seeds" => {
                            if i == 0 {
                                b'C'
                            } else {
                                b"AAAAAAAAAAG"[(i - 1) % 11]
                            }
                        }
                        _ => b"ACGT"[(state & 3) as usize],
                    }
                })
                .collect();
            if pattern == "dimer" {
                for i in length / 2..length {
                    sequence[i] = sequence[i - length / 2];
                }
            }
            if pattern == "false-seeds" && length > 10 {
                sequence[length - 10..].fill(b'A');
            }
            sequence
        })
        .collect();
    let monomerizer = implementation::monomerizer(settings);
    let indexed: Vec<_> = if operation == "orf-indices" {
        sequences
            .iter()
            .map(|s| {
                circkit::orfs::start_stop_codon_indices_by_frame_naive(
                    std::str::from_utf8(s).unwrap(),
                    &["ATG"],
                    &["TAA", "TAG", "TGA"],
                )
            })
            .collect()
    } else {
        Vec::new()
    };
    let orf = circkit::orfs::Orf {
        start: length / 2,
        stop: Some(0),
        wraps: 2,
        length: length * 2,
    };
    macro_rules! dispatch {
        ($consume:ident, $index:expr) => {{
            let index = $index;
            let sequence = black_box(sequences[index].as_slice());
            match operation {
                "lmsr-index" => $consume(implementation::index(sequence, settings)),
                "canonicalize" => $consume(implementation::canonicalize(sequence, settings)),
                "monomerize" => $consume(monomerizer.last_monomer_end_index(sequence)),
                "monomerize-sensitive" => {
                    $consume(monomerizer.last_monomer_end_index_sensitive(sequence))
                }
                "find-orfs" => $consume(circkit::orfs::find_orfs(
                    std::str::from_utf8(sequence).unwrap(),
                )),
                "orf-indices" => {
                    let (s, t) = &indexed[index];
                    $consume(circkit::orfs::find_orfs_with_indices(
                        length,
                        s.clone(),
                        t.clone(),
                    ))
                }
                "extract-orf" => $consume(orf.seq_with_opts(sequence, true)),
                "normalize" => $consume(implementation::normalize(sequence)),
                "hamming-scalar" => $consume(bio::alignment::distance::hamming(
                    sequence,
                    &sequences[(index + 1) % 8],
                )),
                "hamming-simd" => $consume(bio::alignment::distance::simd::hamming(
                    sequence,
                    &sequences[(index + 1) % 8],
                )),
                _ => panic!("Unknown operation: {operation}"),
            }
        }};
    }
    let start = Instant::now();
    for iteration in 0..iterations {
        dispatch!(ignore, iteration % sequences.len());
    }
    let elapsed = start.elapsed().as_nanos();
    // Validate all fixture outputs separately from the timed loop.
    let checksum = (0..sequences.len()).fold(0_u64, |hash, index| {
        hash.wrapping_add(dispatch!(digest, index))
    });
    println!("{{\"elapsed_ns\":{elapsed},\"iterations\":{iterations},\"checksum\":{checksum}}}");
}
