Tools and Rust algorithms for circular DNA and RNA sequences: canonicalization,
deduplication, monomer detection, rotation, and ORFs crossing the sequence origin.

Download the archive matching your OS and CPU. Linux binaries are static musl
builds; macOS binaries require macOS 11 or newer. Every archive includes the
`circkit` executable, installation instructions, CLI guide, command schema,
shell completions, and build metadata.

Verify downloads against `SHA256SUMS`. GitHub artifact attestations are available
for binary archives and crate packages.

After crates.io publication, `cargo install circkit-cli --locked` installs the
CLI; `cargo add circkit` adds the library. Source builds require Rust 1.85+ and a
C compiler. See the included `docs/install.md` for exact installation commands.
