# Bundled native code notices

These are the unmodified upstream notices for native libraries compiled by the
locked dependencies:

| Notice | Source in the crate archive |
| --- | --- |
| bzip2.txt | bzip2-sys 0.1.11+1.0.8, bzip2-1.0.8/LICENSE |
| zstd.txt | zstd-sys 2.0.8+zstd.1.5.5, zstd/LICENSE |
| xz.txt | lzma-sys 0.1.20, xz-5.2/COPYING |

Only liblzma and its public-domain helpers are compiled from XZ Utils. The
separate GPL command-line utilities and getopt implementation are not linked.
Zstd's BSD license is used for redistribution. These notices supplement the
generated Rust dependency license report. Update them when changing bundled
native dependencies. The installed Rust toolchain's COPYRIGHT-library.html
is also included, covering the standard library and its runtime components
(including musl for Linux builds).
