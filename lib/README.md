# circkit

Rust algorithms for circular DNA and RNA sequences: canonical rotation,
monomer detection, and open reading frames that can cross the sequence origin.

```rust
use circkit::canonicalize::{canonicalize, lmsr_index};

assert_eq!(canonicalize(b"TGCA"), b"ATGC");
assert_eq!(lmsr_index(b"TGCA"), 3);
```

Add the library with `cargo add circkit`. Rust 1.85 or newer is required.
See the [API documentation](https://docs.rs/circkit) for rotation options,
`Monomerizer`, and ORF extraction.

The command-line tool is a separate package:
`cargo install circkit-cli --locked` installs the `circkit` executable.

Source code and CLI documentation:
[Benjamin-Lee/circkit](https://github.com/Benjamin-Lee/circkit).
Licensed under the MIT license.
