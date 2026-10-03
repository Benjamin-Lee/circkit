# Performance and portable benchmarking

Run on the machine whose performance matters. Python 3.9+, Git, a native Rust
toolchain with Cargo, and the repository's history are required. Linux and macOS
are supported; the runner uses only Python's standard library. Cargo may download
dependencies on the first build. In a shallow clone, run `git fetch --unshallow`
so the baseline commit is available.

```sh
python3 benchmarks/run.py --suite full --threads 1,2,4 --output target/benchmarks-arm \
  --note "ARM laptop, plugged in, default power mode"
```

Start with `--suite smoke --trials 2 --no-tuning` to check the build and output
validation. `quick` uses smaller CLI datasets; `full` reduces startup effects and
covers more sequence lengths. Builds can take several minutes. Benchmark with no
competing builds or workloads. The runner leaves CPU scheduling to the OS; optional
`--cpus 0` pins the entire run on Linux. macOS does not need an affinity utility.

Send back `target/benchmarks-arm/report.zip`. It contains:

- `results.csv` and `results.json`: medians, quartiles, standard deviations, raw
  ranges, settings, and speedups relative to the named reference stage.
- `raw.json`: each timed trial, its configuration, and library loop count.
- `environment.json`: OS, architecture, CPU description/features, logical CPUs,
  affinity and Linux CPU quota when available, Python version, seed, and run options.
- `builds.json`: source fingerprint, baseline/candidate commit IDs, dirty status,
  binary hashes, Rust/Cargo versions, dependency versions, features, and build flags.
- `inputs.json`: input hashes, sizes, record counts, and sequence-length statistics.
- `validation.json` and build logs: output comparisons and build evidence.

Input sequences and source/build directories are retained locally and excluded
from the report archive. Failed runs exit nonzero and retain partial results.

## What is compared

The default baseline is merged PR #2, commit
`deb0b407832cfc276ccf351f38f9201b90f90d34`. Use `--baseline REF` for another
compatible revision. The candidate is a snapshot of the current checkout,
including uncommitted files. Both builds use the same installed compiler, locked
dependencies, and release settings: fat LTO and one codegen unit. Environment
overrides are recorded. No `target-cpu=native` flag is added. A third build with
`scalar-normalization` measures the fallback for the new normalization fast path.

The shared Rust driver uses small adapters to call each revision's public API.
Library timings exclude preparation, output checksums, and process startup. CLI
timings include startup, parsing, processing, and output; plain output goes to the
OS null device, while gzip output goes to a real temporary file. Input caches are
warm. Stages are interleaved in a seeded random order. Library loop counts are
calibrated independently to avoid making slow stages dominate the run.

Before timing each configuration, the runner compares library checksums, FASTA
record multisets including full headers, and metadata row multisets. Multiple
workers can reorder FASTA batches. `uniq` stays at one worker in comparisons
because parallel deduplication can select a different representative ID.

Deterministic fixtures cover short reads, dimers, wrapped medium sequences, long
contigs, dense codons, and false seed matches. Dense codons and false seeds are
stress cases, not representative genomics estimates. Add local plain/gzip FASTA
with repeated `--input /path/to/sequences.fasta` options; these files are never
modified. Empty or very small outputs still count in equivalence checks.

`--skip-build` reuses binaries only when the source fingerprint, baseline,
toolchain, relevant build flags, and binary hashes match. Run with the same
`--output` directory. `--no-scalar` skips the fallback build; `--no-tuning` skips
parameter sweeps. `--offline` requires all Cargo dependencies to be cached.
Use `python3 benchmarks/run.py --help` for all options.

## Empirical tuning controls

These controls affect execution cost, not biological matching criteria:

| CLI control | Default | Meaning |
| --- | --- | --- |
| `--rotation-cutoff BYTES` | 32768 | `canonicalize`/`uniq`: contiguous allocating Duval at or below this length, constant-space two-candidate search above it. Zero forces constant-space search for nonempty inputs; a sufficiently large value forces Duval. |
| `--mismatch-chunk-size BYTES` | 256 | `monomerize`: check approximate overlaps in chunks and stop once mismatches exceed the permitted count. Must be positive. Exact matches use byte equality; debug logging computes the full distance. |
| `--queue-depth BUFFERS` | 64; ORFs: max(2, 2*threads) | Bound queued FASTA work in pipeline execution. Must be positive. Legacy `monomerize --batch-size` remains an alias. |
| `--execution auto\|serial\|pipeline` | auto | Auto uses serial processing with one worker and one available CPU, otherwise a reader/worker/writer pipeline. Serial requires `--threads 1`. |
| `--threads N` | logical CPUs | Existing worker count, now validated to be positive. |

For example:

```sh
circkit canonicalize sequences.fasta --threads 1 --execution serial --rotation-cutoff 65536
circkit monomerize sequences.fasta --threads 4 --mismatch-chunk-size 128 --queue-depth 8
```

The library exposes `canonicalize::RotationOptions { duval_max_len }` methods and
`Monomerizer::builder().mismatch_chunk_size(NonZeroUsize)`. Existing rotation
functions keep the default settings. Builder-based monomerizer callers keep
working; struct-literal callers must supply the new `mismatch_chunk_size` field.

Default benchmark sweeps compare rotation cutoffs 0/32768/1048576, mismatch chunks
64/256/1024, serial/pipeline execution, queue depths 2/8/64, and the selected
worker counts. Override them with `--rotation-cutoffs`, `--mismatch-chunks`,
`--queue-depths`, and `--threads`. Rotation sweeps time the public index API;
canonicalization and buffer behavior are also covered by equivalence tests.

The cutoff and chunk defaults were chosen from x86 measurements and remain
provisional for ARM. Algorithmic ORF improvements and allocation reductions are
portable. The normalization fast path checks AVX2 at runtime on x86 and falls back
to scalar code on other CPUs; there is no new ARM NEON kernel. memchr and zlib-rs
use their own platform support. Vector widths and the 64-entry codon table are
architecture/biology constants and are not tuning parameters.

Native ARM Linux and macOS CI runs check correctness and smoke-test this runner.
Compare medians and variability on your own workload before changing defaults.
The benchmark does not measure peak memory, cold-cache behavior, or energy use,
and it does not promise a speedup for every command or worker count.
