# Command-line guide

The CLI uses Clap 4 and requires Rust 1.85+. Existing biological commands and flags
remain available. `circkit --help` lists commands; `circkit COMMAND --help` explains
options, defaults, and accepted values.

## FASTA streams and files

```sh
circkit canonicalize reads.fasta -o canonical.fasta
circkit monomerize reads.fasta.gz --keep-all --threads 1 | circkit uniq --threads 1
circkit rotate - --bases -5 -o - < reads.fasta
```

Omitting the input or using `-` reads stdin. Omitting `-o` or using `-o -` writes
FASTA to stdout. Use `./-` for a file literally named `-`. Paths containing spaces
should be quoted. Compression is detected from input bytes, including stdin;
output filenames ending in `.gz`, `.bz2`, `.xz`, or `.zst` select compression.
Stdout is plain FASTA. Headers are preserved; sequence line wrapping can change.

Empty input and empty records are accepted. The FASTA parser reports invalid
record syntax, but it is not a validator for biological alphabets. Canonicalization,
deduplication, and ORF finding normalize sequences as before. `cat`, `decat`, and
`rotate` preserve sequence bytes. `decat` keeps the first half, rounding down.

Input and output files must be distinct, including symbolic and hard links. FASTA
and metadata outputs must also differ. These checks run before any output is
opened. Named output files are overwritten; output is streamed, so an error may
leave partial output. Check the exit status before accepting results.

Output is buffered and stdout is locked for processing. Compressed output is
explicitly finalized, including trailers, and write errors are propagated. Closing
stdout early (for example, piping into `head`) ends processing quietly with exit 0.
A broken named output or metadata pipe remains an error. Early stdout closure can
also leave a partial named output when metadata is being piped.

## Metadata

`monomerize`, `uniq`, and `orfs` support `--table FILE`:

```sh
circkit monomerize reads.fasta --keep-all --table monomers.csv -o monomers.fasta
circkit orfs reads.fasta --table orfs.jsonl.gz -o orfs.fasta
circkit uniq reads.fasta --table - --table-format jsonl -o unique.fasta
```

`.tsv` selects tab-separated values, `.jsonl` selects one JSON object per line,
and other suffixes select CSV. Compression suffixes are stripped before choosing
the table format, so `.tsv.gz` and `.jsonl.gz` work. `--table-format csv|tsv|jsonl`
overrides format selection and requires `--table`. `--table -` requires a named
FASTA output (`-o FILE`) so the two streams do not mix. No rows means an empty
metadata file (or a compressed representation of an empty file).

| Command | Fields | Meaning |
| --- | --- | --- |
| `monomerize` | `id`, `original_length`, `monomer_length` | One row per retained sequence. ID is the full header. |
| `uniq` | `id`, `duplicate_id` | One row per discarded duplicate, referring to the retained representative ID (first header token). |
| `orfs` | `orf_id`, `seq_id`, `start`, `stop`, `length`, `wraps`, `ratio` | One row per output ORF. IDs use the full FASTA header. |

Metadata identifiers must be UTF-8; FASTA-only processing accepts byte headers.
JSONL uses numbers for lengths, positions, wraps, and ratios. ORF positions are
zero-based; a missing stop is JSON `null` or an empty CSV/TSV cell. `length` follows
`--include-stop`; `ratio` includes the stop codon regardless of that flag.
`--min-length` for ORFs filters coding length **excluding** the stop codon.

## Discoverability and automation

```sh
circkit schema              # Full JSON command catalog
circkit schema orfs         # Options for one command
circkit describe canon      # Aliases also work
circkit --error-format json orfs missing.fasta -o orfs.fasta
```

`schema` returns a JSON object with `schema_version: 1`, the program version,
commands, stream contracts, metadata formats, ordering guidance, and exit codes.
Command descriptions come from the same Clap definitions used for parsing. Each
command includes its name, aliases, usage, arguments, groups, and metadata fields.
Arguments describe flag spellings, positions, value types, possible enum values,
defaults, arity, required/global status, conflicts, and help. Required exclusive
groups describe alternatives such as `rotate --bases` or `--percent`. Numeric
ranges and conditional requirements are explained in help and enforced by the
parser; the catalog is not a JSON Schema for constructing arbitrary argument lists.
Agents should read the catalog, choose a command, and supply individual argument
values through their subprocess API rather than interpolate paths into a shell.

FASTA remains the primary data format. Diagnostics always use stderr; normal runs
are quiet. Global `-v` increases logging and `-q` decreases it. `--error-format json`
selects a versioned JSON error object on stderr, including argument-parse failures:

```json
{"schema_version":1,"error":{"code":"io_error","message":"open input missing.fasta: No such file or directory (os error 2)","exit_code":1}}
```

| Exit code | Meaning |
| --- | --- |
| 0 | Success, help/version, or stdout closed by a downstream consumer. |
| 1 | Input, output, or processing failure. |
| 2 | Invalid arguments. |

Error codes are `invalid_arguments`, `invalid_fasta`, `invalid_header`, `io_error`,
and `processing_error`. Messages provide context and may vary by OS. Additional
verbose logging is textual, so omit `-v` when consuming stderr as a JSON error.
Successful help/version output stays human-readable regardless of error format.

For reproducible record order and duplicate representatives, use `--threads 1`.
Multiple workers may reorder batches and select different duplicate representatives.
Performance controls and the portable benchmark/report workflow are described in
[the performance guide](../benchmarks/README.md).

## Input validation and compatibility

`rotate` requires exactly one of nonzero `--bases INTEGER` or finite, nonzero
`--percent FRACTION`. Positive values rotate right, negative values left; rotation
wraps modulo the sequence length. A percentage is a decimal fraction (`0.5` means
50%); its base shift is rounded down. Values whose computed shift exceeds the signed
64-bit range are rejected.

`--threads` and queue/chunk sizes must be positive. Serial execution requires
`--threads 1`. Identities must be between 0 and 1; overlap/ORF ratios must be finite
and nonnegative. Minimum lengths/wrap counts cannot exceed their maxima; wrap
counts range from 0 to 3. Codon lists accept comma-separated DNA triplets using
A/C/G/T/N, normalize case and surrounding spaces, and reject malformed entries.
N matches a literal N, not a wildcard.

`orfs --strand reverse` searches and outputs only reverse-strand ORFs; earlier
versions also included forward-strand ORFs in that mode. Default `both` behavior
and the ORF matching algorithms are unchanged.

`canon` aliases `canonicalize`; `--min-overlap-ratio` aliases the existing
`--min-overlap-percent` flag. `--batch-size` remains a queue-depth alias, and
`uniq --norm`/`--canon` remain canonical-output aliases.

## Shell completion

```sh
circkit completions bash > circkit.bash
source circkit.bash
circkit completions fish > circkit.fish
circkit completions zsh > _circkit
```

Bash, Fish, Zsh, PowerShell, and Elvish scripts are generated from the parser.
Source or install the generated script according to your shell's completion setup.

## Design references

The stream and error conventions follow the [Command Line Interface Guidelines](https://clig.dev/).
Machine-readable discovery follows patterns used by the
[GitHub CLI's JSON interface](https://cli.github.com/manual/gh_help_formatting) and
[Google Workspace CLI's runtime schema discovery](https://github.com/googleworkspace/cli).
The latter project describes itself as not officially supported by Google.
