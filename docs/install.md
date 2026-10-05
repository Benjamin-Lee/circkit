# Installation

The CLI package is **circkit-cli**; the executable and Rust library are both
named **circkit**.

## Prebuilt binaries

Download an archive and its `.sha256` file from
[GitHub Releases](https://github.com/Benjamin-Lee/circkit/releases).
The initial release will appear there after publication.

| Computer | Target in the archive filename |
| --- | --- |
| Linux, Intel/AMD 64-bit | `x86_64-unknown-linux-musl` |
| Linux, ARM 64-bit | `aarch64-unknown-linux-musl` |
| macOS, Intel | `x86_64-apple-darwin` |
| macOS, Apple Silicon | `aarch64-apple-darwin` |

Run `uname -s` and `uname -m` if you are unsure. macOS binaries require macOS 11
or newer. Linux binaries include their C runtime and compression libraries;
they work on glibc and musl distributions without installing those libraries.
The CPU instruction set is the target's baseline, with AVX2 selected only when
available at runtime. Windows users can use the Linux binaries under WSL2.

For an Apple Silicon Mac, for example, after downloading the two files:

```sh
shasum -a 256 -c circkit-0.1.0-aarch64-apple-darwin.tar.gz.sha256
tar -xzf circkit-0.1.0-aarch64-apple-darwin.tar.gz
mkdir -p "$HOME/.local/bin"
install -m 755 circkit-0.1.0-aarch64-apple-darwin/bin/circkit "$HOME/.local/bin/circkit"
export PATH="$HOME/.local/bin:$PATH"
circkit --version
```

On Linux, `sha256sum -c ARCHIVE.tar.gz.sha256` verifies the download.
Replace the version and target above with those you downloaded. Add the PATH
line to your shell configuration to retain it in new terminals.
GitHub's `SHA256SUMS` covers all binary and crate archives. With the GitHub CLI,
you can also verify build provenance:

```sh
gh attestation verify ARCHIVE.tar.gz --repo Benjamin-Lee/circkit
```

Each archive includes a CLI guide, a JSON command catalog, shell completions,
and `build.json` recording the source commit, target, compiler, and binary hash.
Generate fresh discovery data with `circkit schema` or `circkit completions zsh`.
See [the CLI guide](cli.md) for completion installation and streaming examples.

macOS binaries have the normal ad hoc signature produced by the Rust linker;
they are not Developer ID signed or notarized. macOS may require you to allow
a downloaded executable under System Settings > Privacy & Security. Verify
the download first. Cargo installation builds locally and avoids that step.

## Cargo installation

After the first crates.io publication:

```sh
cargo install circkit-cli --version 0.1.0 --locked
circkit --help
```

Rust 1.85 or newer and a C compiler are required. A current Rust toolchain is
recommended. macOS: install Apple's command-line developer tools with
`xcode-select --install`. Debian/Ubuntu: install `build-essential` and
`pkg-config`. Compression libraries can be built from bundled source; use
`--features static-codecs` to force bundling of bzip2 and xz.

The package name is `circkit-cli`, not `circkit`. The latter is the library:

```sh
cargo add circkit
```

## Build from a checkout

Before publication, or for development:

```sh
cargo build --release --locked
cargo test --workspace --locked
cargo install --path . --locked
```

If a different Rust installation appears active, inspect `command -v rustc`,
`rustc --version`, and `rustup show`. With rustup, update using `rustup update
stable`; upgrading Homebrew's Rust does not update a rustup toolchain selected
earlier on PATH or by an override.

## Update and remove

To update a downloaded binary, verify and install the new archive in the same
location. Remove `$HOME/.local/bin/circkit` to uninstall it. For Cargo installs,
use `cargo install circkit-cli --locked --force` to update, or
`cargo uninstall circkit-cli` to remove it.
