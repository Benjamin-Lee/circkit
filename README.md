# circKit

circKit is a library for manipulating circular biological sequences such as DNA and RNA.

## Features

- Easy to install
- Written in Rust for performance and safety
- Inputs and outputs can be gzip, bzip2, xz, or zstd compressed
- Streaming FASTA with CSV, TSV, or JSONL metadata
- Shell completions, JSON command discovery, and structured errors for automation

## Usage

Building the CLI requires Rust 1.85 or newer:

```sh
cargo build --release --locked
./target/release/circkit --help
```

See the [CLI guide](docs/cli.md) for streaming, metadata, shell completions,
and automation examples. Run `circkit COMMAND --help` for each command's options.

```text
$ circkit --help
A toolkit for working with circular sequences.

Usage: circkit [OPTIONS] <COMMAND>

Commands:
  monomerize    Find monomers of (potentially) circular or multimeric sequences
  cat           Concatenate sequences to themselves
  decat         Deconcatenate sequences to themselves
  canonicalize  Select a consistent rotation and strand for circular sequences [alias: canon]
  uniq          Deduplicate circular sequences
  rotate        Rotate circular sequences to the left or right
  orfs          Find ORFs in circular sequences
  schema        Print a versioned JSON command catalog for automation and agents [alias: describe]
  completions   Generate shell completion scripts on stdout
  help          Print this message or the help of the given subcommand(s)

Options:
  -v, --verbose...                   Increase logging verbosity
  -q, --quiet...                     Decrease logging verbosity
      --error-format <ERROR_FORMAT>  Diagnostic format on stderr; FASTA output is unaffected [default: text] [possible values: text, json]
  -h, --help                         Print help
  -V, --version                      Print version

Examples:
  circkit canonicalize reads.fasta -o canonical.fasta
  circkit monomerize - --keep-all --threads 1
  circkit schema orfs

Use '-' for stdin/stdout. Run 'circkit COMMAND --help' for options.
```

## Subcommands

### `cat` and `decat`

`cat` and `decat` are used to concatenate and deconcatenate sequences to themselves. These commands are useful for dealing with tools that assume that the input sequence is linear.

For example, let's say you had a tool that does some sort of filtering on linear sequences. You could use `cat` to convert your circular sequences to linear sequences, run the tool, and then use `decat` to convert the output linear sequences back to circular sequences.

As a note, circKit's `cat` is functionally equivalent to `seqkit concat file.fasta file.fasta`.

### `canonicalize`

`canonicalize` computes a single, canonical representation of a circular sequence. This is useful for comparing sequences that are the same but have different polarities or different starting positions. For example, the following two sequences are the same if they are circular:

```text
>seq1
TGCA
>seq2
GCAT
```

We define the canonical representation as the lexicographically smallest rotation of either polarity. In other words, we compute the sequence rotation that would come first in the alphabet for each polarity (known as the [lexicographically minimal string rotation](https://en.wikipedia.org/wiki/Lexicographically_minimal_string_rotation)). Then we simply compare the LMSRs for each polarity and return the one that comes first in the alphabet. So, in the example above, the normalized representation is `ATGC`.

## Roadmap

There's still a lot to do before an initial release.
Here's how it's going:

- [x] `rotate`
- [x] `cat`
- [x] `decat`
- [x] `canonicalize`
- [x] `monomerize`
- [x] `orfs`
- [ ] ~~`grep`~~ (use `seqkit grep` instead)
- [ ] `cluster` (future)
- [ ] `prealign` (future)

## Benchmarks

See [the performance guide](benchmarks/README.md) for portable benchmarks,
machine-readable results, and tuning controls. The runner compares the merged
PR #2 baseline with the current checkout on the same machine.

## See Also

- [vdsearch](https://github.com/Benjamin-Lee/vdsearch): A tool for searching for viroid-like sequences. Eventually, `circkit` be used for all of the data manipulation tasks.
