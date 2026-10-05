# Maintainer release procedure

The initial version is **0.1.0** for both crates. `circkit` is the library;
`circkit-cli` depends on it and installs the `circkit` executable.

## What CI verifies

Every PR runs a publishing dry run for both crates and native builds/tests on
Linux and macOS, each on x86-64 and ARM64. The binary matrix and pinned release
compiler live in [dist/config.toml](../dist/config.toml). Release tooling needs
Python 3.12+; end users do not need Python or Rust for downloaded binaries.

GNU/Linux builds use digest-pinned manylinux 2.28 containers, with bundled
compression libraries. CI rejects non-system dependencies and glibc symbol
requirements above 2.28. Static musl builds provide an additional compatibility
option; CI rejects a dynamic loader or shared-library dependencies in those
builds. macOS targets version 11 and
bundles non-system codecs; CI rejects dependencies outside standard system
paths. No `target-cpu=native` flags are used. Each archive is extracted and
tested for canonicalization, monomers/JSONL metadata, command discovery,
completions, and round trips through gzip, bzip2, xz, and zstd.

Archives contain MIT licensing for circkit, a generated dependency license
report, bundled native codec notices, and the Rust toolchain's runtime notices.
GNU builds preserve the host glibc allocator, which can be much faster than
musl for multithreaded short-read processing. Choose GNU/Linux on supported
distributions and musl for Alpine/older Linux. This changes distribution
choices rather than algorithm tuning.

The license generator and its download checksum are pinned in the release
config. Update `dist/native-licenses` when changing bundled native codecs.

`cargo publish --workspace --dry-run --locked` verifies extracted crate sources
in dependency order, using Cargo's staging registry before the library exists
on crates.io. The library's packaged tests and README example are tested too.
The CLI source package omits the 77 MB integration fixtures and local benchmark
driver. The full checkout tests still run on every release target.

## First release

1. Merge prerequisite PRs and the packaging PR. Confirm CI is green on main.
2. Check both Cargo.toml versions and the CLI's `circkit` dependency agree.
   Run `python3 dist/release.py check --tag v0.1.0`.
3. Create the version tag on the tested main commit and push it:

   ```sh
   git switch main
   git pull --ff-only
   git tag -a v0.1.0 -m 'circkit 0.1.0'
   git push origin v0.1.0
   ```

4. Wait for **Release packages** to pass. It creates a **draft** GitHub release
   with six native binary archives, two `.crate` archives, checksums, and
   provenance attestations. Review its artifacts and notes. It does not publish
   crates or make the draft public. A failed job creates no draft release.
5. Publish the initial crates locally from the tagged checkout. crates.io needs
   this first publication before trusted publishing can be configured:

   ```sh
   git switch --detach v0.1.0
   rustup toolchain install 1.99.0 --profile minimal
   cargo +1.99.0 publish --workspace --dry-run --locked
   cargo login
   cargo +1.99.0 publish --workspace --locked
   ```

   Use a crates.io account you control and a short-lived token that can create
   these two crate names. Cargo publishes the library before the CLI and waits
   for dependencies to become available. The names were available when this
   workflow was prepared; publication reserves them. If publication stops
   after the library succeeds, finish with
   `cargo +1.99.0 publish -p circkit-cli --locked`.
6. Test the registry installation with
   `cargo install circkit-cli --version 0.1.0 --locked --root /tmp/circkit-install-check`
   and run its `bin/circkit --version` and `schema` commands.
7. Publish the reviewed GitHub draft. No automated workflow performs this step.
8. On each crate's crates.io settings page, add a GitHub trusted publisher for
   owner **Benjamin-Lee**, repository **circkit**, workflow **publish.yml**,
   with no environment restriction. This avoids long-lived repository secrets.
   See [crates.io trusted publishing](https://crates.io/docs/trusted-publishing).

The CLI minimum compiler remains Rust 1.85. Release packaging uses the compiler
pinned in `dist/config.toml`; Cargo's workspace publishing support needs Cargo
1.90 or newer. Recheck the pinned compiler and all platform jobs when upgrading.
The CLI/library versions need not increase for this packaging PR because
neither package has been published yet.

## Later releases

Update both crate versions and the CLI library dependency together, update
version examples in the installation guide, and regenerate Cargo.lock. Tag a
tested commit on main. The same workflow builds and creates a draft release.
Then run **Publish crates** manually on main with the tag as input. It requires
all release artifacts, reruns tests and a publishing dry run, obtains a
short-lived crates.io token, and publishes both crates in dependency order.
Publish the GitHub draft after the registry installation check.

Versions on crates.io are immutable. Completed public GitHub releases are
never overwritten by this workflow. Rerunning a tag build can update an
existing draft. If publishing fails partway through, use a local scoped token
to publish the missing crate; do not rerun a whole-workspace publication when
one crate already exists.

## Local artifact checks

To reproduce a Linux x86-64 archive, install a native musl compiler first
(`musl-tools` on Debian/Ubuntu):

```sh
rustup toolchain install 1.99.0 --profile minimal --target x86_64-unknown-linux-musl
export RUSTUP_TOOLCHAIN=1.99.0
cargo install cargo-about --version 0.9.2 --locked
cargo about generate --workspace --all-features --locked --fail --config dist/about.toml dist/licenses.hbs --output-file target/THIRD-PARTY-LICENSES.html
CC=musl-gcc cargo +1.99.0 build --release --locked --features static-codecs --target x86_64-unknown-linux-musl --bin circkit
python3 dist/release.py pack --binary target/x86_64-unknown-linux-musl/release/circkit --target x86_64-unknown-linux-musl --licenses target/THIRD-PARTY-LICENSES.html
python3 dist/release.py verify target/dist/circkit-0.1.0-x86_64-unknown-linux-musl.tar.gz
```

For macOS, use the native Apple target above and set
`MACOSX_DEPLOYMENT_TARGET=11.0` during the build. The packer runs on the target's
native architecture. `rustc` on PATH must be the compiler used to build the
binary so `build.json` records it accurately.

Archive timestamps, ordering, ownership, and modes are normalized. Builds use
a pinned compiler and lockfile; byte-identical binaries across different
machines are not guaranteed. Check each archive's build metadata and attestation
to identify its source and platform. macOS Developer ID signing/notarization,
Homebrew/Bioconda recipes, containers, and native Windows binaries are separate
distribution channels, not prerequisites for these Cargo/GitHub releases.
