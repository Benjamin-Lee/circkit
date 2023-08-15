# circKit

circKit is a library for manipulating circular biological sequences such as DNA and RNA.

## Features

- Easy to install
- Written in Rust for performance and safety
- Inputs and outputs can be gzip, bzip2, xz, or zstd compressed

## Usage

```text
$ circkit --help
circkit 0.1.0
Benjamin D. Lee <benjamin.lee@chch.ox.ac.uk>

USAGE:
    circkit [OPTIONS] <SUBCOMMAND>

OPTIONS:
    -h, --help       Print help information
    -q, --quiet      Less output per occurrence
    -v, --verbose    More output per occurrence
    -V, --version    Print version information

SUBCOMMANDS:
    canonicalize    Normalize circular sequences [aliases: canon]
    cat             concatenate sequences to themselves [aliases: concat, concatenate]
    decat           deconcatenate sequences to themselves (this is the reverse of `cat`)
                        [aliases: deconcat, deconcatenate, uncat, unconcatenate]
    help            Print this message or the help of the given subcommand(s)
    monomerize      Find monomers of (potentially) circular or multimeric sequences
    orfs            Find ORFs in circular sequences
    rotate          Rotate circular sequences to the left or right
    uniq            Deduplicate circular sequences
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

To be done.

## See Also

- [vdsearch](https://github.com/Benjamin-Lee/vdsearch): A tool for searching for viroid-like sequences. Eventually, `circkit` be used for all of the data manipulation tasks.
