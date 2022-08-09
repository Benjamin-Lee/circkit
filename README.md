# circkit

## Introduction

Circkit is a library for manipulating circular biological sequences such as DNA and RNA.

## Features

- Easy to install
- Written in Rust for performance and safety
- Inputs can be gzip, bzip2, or xz compressed

## Usage

```text
$ circkit --help
circkit 0.1.0
Benjamin D. Lee <benjamin.lee@chch.ox.ac.uk>

USAGE:
    circkit <SUBCOMMAND>

OPTIONS:
    -h, --help       Print help information
    -V, --version    Print version information

SUBCOMMANDS:
    cat           concatenate sequences to themselves [aliases: concat, concatenate]
    decat         deconcatenate sequences to themselves (this is the reverse of `cat`) [aliases:
                      deconcat, deconcatenate, uncat, unconcatenate]
    help          Print this message or the help of the given subcommand(s)
    monomerize    Find monomers of (potentially) circular or multimeric sequences
    normalize     Normalize circular sequences [aliases: canonicalize, canon]
```

## Roadmap

There's still a lot to do before an initial release.
Here's how it's going:

- [ ] `rotate`
- [x] `cat`
- [x] `decat`
- [x] `normalize`
- [x] `monomerize`
- [ ] `grep`
- [ ] `cluster` (future)
- [ ] `prealign` (future)
