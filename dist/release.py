#!/usr/bin/env python3
"""Build-independent release packaging and checks (Python 3.12+, stdlib only)."""

import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import tarfile
import tempfile
import tomllib

ROOT = Path(__file__).resolve().parents[1]
SHELLS = {"bash": "circkit.bash", "fish": "circkit.fish", "zsh": "_circkit",
          "powershell": "circkit.ps1", "elvish": "circkit.elv"}


def run(*args, **kwargs):
    return subprocess.run(args, check=True, capture_output=True, **kwargs).stdout


def configuration(tag=None):
    cli = tomllib.loads((ROOT / "Cargo.toml").read_text())
    lib = tomllib.loads((ROOT / "lib/Cargo.toml").read_text())
    config = tomllib.loads((ROOT / "dist/config.toml").read_text())
    version = cli["package"]["version"]
    if not re.fullmatch(r"\d+\.\d+\.\d+", version):
        raise ValueError("release versions must have the form X.Y.Z")
    if lib["package"]["version"] != version:
        raise ValueError("library and CLI versions must match")
    if cli["dependencies"]["circkit"]["version"] != version:
        raise ValueError("CLI dependency must match the library version")
    if (ROOT / "LICENSE").read_bytes() != (ROOT / "lib/LICENSE").read_bytes():
        raise ValueError("crate licenses must match")
    if tag and tag != f"v{version}":
        raise ValueError(f"tag {tag!r} does not match v{version}")
    return version, config


def sha256(path):
    with path.open("rb") as source:
        return hashlib.file_digest(source, "sha256").hexdigest()


def check_linkage(binary, target):
    """Reject dynamic Linux binaries and non-system macOS dependencies."""
    architecture = target.split("-", 1)[0]
    native = {"arm64": "aarch64", "AMD64": "x86_64"}.get(platform.machine(), platform.machine())
    if architecture != native:
        raise ValueError(f"{target} must be smoke-tested on its native architecture")
    if target.endswith("linux-musl"):
        headers = run("readelf", "-l", str(binary)).decode()
        dynamic = run("readelf", "-d", str(binary)).decode()
        if "INTERP" in headers or "NEEDED" in dynamic:
            raise ValueError("Linux release binary must be fully static")
    else:
        libraries = run("otool", "-L", str(binary)).decode().splitlines()[1:]
        for library in libraries:
            name = library.strip().split(" (", 1)[0]
            if not name.startswith(("/usr/lib/", "/System/Library/")):
                raise ValueError(f"unbundled macOS dependency: {name}")


def smoke(binary, version):
    """Exercise the installed binary, including every compressed stream codec."""
    binary = str(binary.resolve())
    if run(binary, "--version").decode().strip() != f"circkit {version}":
        raise ValueError("binary version does not match the package")
    schema = json.loads(run(binary, "schema"))
    if schema["schema_version"] != 1 or not schema["commands"]:
        raise ValueError("invalid command catalog")
    fasta = b">circle\nTGCA\n"
    canonical = b">circle\nATGC\n"
    if run(binary, "canonicalize", "-", "--threads", "1", input=fasta) != canonical:
        raise ValueError("canonicalization smoke test failed")
    with tempfile.TemporaryDirectory(prefix="circkit-smoke-") as temporary:
        temporary = Path(temporary)
        for suffix in ("gz", "bz2", "xz", "zst"):
            compressed = temporary / f"sequences.fasta.{suffix}"
            run(binary, "canonicalize", "-", "--threads", "1", "-o", str(compressed), input=fasta)
            if run(binary, "canonicalize", str(compressed), "--threads", "1") != canonical:
                raise ValueError(f"{suffix} round trip failed")
        table = temporary / "monomers.jsonl"
        output = temporary / "monomers.fasta"
        run(binary, "monomerize", "-", "--threads", "1", "--keep-all", "--table", str(table),
            "-o", str(output), input=b">dimer\nATGCACTGGAATGCACTGGA\n")
        rows = [json.loads(line) for line in table.read_text().splitlines()]
        if rows != [{"id": "dimer", "original_length": 20, "monomer_length": 10}]:
            raise ValueError("JSONL metadata smoke test failed")
        if output.read_bytes() != b">dimer\nATGCACTGGA\n":
            raise ValueError("monomer smoke test failed")
    for shell in SHELLS:
        if not run(binary, "completions", shell):
            raise ValueError(f"empty {shell} completions")


def pack(binary, target, output, licenses):
    version, config = configuration()
    if target not in [item["target"] for item in config["targets"]]:
        raise ValueError(f"unsupported release target: {target}")
    binary = binary.resolve()
    check_linkage(binary, target)
    smoke(binary, version)
    output.mkdir(parents=True, exist_ok=True)
    name = f"circkit-{version}-{target}"
    archive = output / f"{name}.tar.gz"
    commit = run("git", "-C", str(ROOT), "rev-parse", "HEAD").decode().strip()
    epoch = int(run("git", "-C", str(ROOT), "show", "-s", "--format=%ct", "HEAD"))
    with tempfile.TemporaryDirectory(prefix="circkit-pack-") as temporary:
        stage = Path(temporary) / name
        (stage / "bin").mkdir(parents=True)
        (stage / "docs").mkdir()
        (stage / "completions").mkdir()
        (stage / "licenses").mkdir()
        shutil.copyfile(binary, stage / "bin/circkit")
        (stage / "bin/circkit").chmod(0o755)
        for source in ("README.md", "LICENSE", "docs/cli.md", "docs/install.md"):
            shutil.copyfile(ROOT / source, stage / source)
        if not licenses.is_file() or licenses.stat().st_size == 0:
            raise ValueError("generate the third-party license report before packaging")
        shutil.copyfile(licenses, stage / "THIRD-PARTY-LICENSES.html")
        shutil.copytree(ROOT / "dist/native-licenses", stage / "licenses", dirs_exist_ok=True)
        sysroot = Path(run("rustc", "--print", "sysroot").decode().strip())
        shutil.copyfile(sysroot / "share/doc/rust/COPYRIGHT-library.html", stage / "licenses/rust-runtime.html")
        for shell, filename in SHELLS.items():
            (stage / "completions" / filename).write_bytes(run(str(binary), "completions", shell))
        (stage / "command-schema.json").write_bytes(run(str(binary), "schema"))
        metadata = {"version": version, "target": target, "commit": commit,
                    "binary_sha256": sha256(binary), "rustc": run("rustc", "-vV").decode(),
                    "features": ["static-codecs"],
                    "macos_minimum": config["macos_deployment_target"] if "apple" in target else None}
        (stage / "build.json").write_text(json.dumps(metadata, indent=2) + "\n")
        # Normalize archive order, timestamps, ownership and permissions.
        with archive.open("wb") as destination, gzip.GzipFile(fileobj=destination, mode="wb", mtime=0, filename="") as compressed:
            with tarfile.open(fileobj=compressed, mode="w") as tar:
                for path in [stage, *sorted(stage.rglob("*"))]:
                    info = tar.gettarinfo(str(path), str(path.relative_to(stage.parent)))
                    info.uid = info.gid = 0
                    info.uname = info.gname = ""
                    info.mtime = epoch
                    info.mode = 0o755 if path.is_dir() or path == stage / "bin/circkit" else 0o644
                    if path.is_file():
                        with path.open("rb") as source:
                            tar.addfile(info, source)
                    else:
                        tar.addfile(info)
    archive.with_suffix(archive.suffix + ".sha256").write_text(f"{sha256(archive)}  {archive.name}\n")
    verify(archive)
    print(archive)


def verify(archive):
    expected = archive.with_suffix(archive.suffix + ".sha256").read_text().strip()
    if expected != f"{sha256(archive)}  {archive.name}":
        raise ValueError("archive checksum mismatch")
    with tempfile.TemporaryDirectory(prefix="circkit-unpack-") as temporary:
        temporary = Path(temporary)
        with tarfile.open(archive) as tar:
            tar.extractall(temporary, filter="data")
        roots = list(temporary.iterdir())
        if len(roots) != 1 or not roots[0].is_dir():
            raise ValueError("archive must contain one root directory")
        stage = roots[0]
        metadata = json.loads((stage / "build.json").read_text())
        binary = stage / "bin/circkit"
        if sha256(binary) != metadata["binary_sha256"]:
            raise ValueError("installed binary checksum mismatch")
        if not os.access(binary, os.X_OK):
            raise ValueError("installed binary is not executable")
        check_linkage(binary, metadata["target"])
        smoke(binary, metadata["version"])


def manifest(output):
    version, config = configuration()
    names = [f"circkit-{version}-{item['target']}.tar.gz" for item in config["targets"]]
    names += [f"{crate}-{version}.crate" for crate in ("circkit", "circkit-cli")]
    for name in names:
        path = output / name
        if not path.is_file():
            raise ValueError(f"missing release artifact: {name}")
        checksum = f"{sha256(path)}  {name}\n"
        sidecar = output / f"{name}.sha256"
        if sidecar.exists() and sidecar.read_text() != checksum:
            raise ValueError(f"checksum mismatch: {name}")
        sidecar.write_text(checksum)
    (output / "SHA256SUMS").write_text("".join((output / f"{name}.sha256").read_text() for name in names))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    check = commands.add_parser("check", help="check versions/licenses and print the build matrix")
    check.add_argument("--tag")
    package = commands.add_parser("pack", help="smoke-test and archive a native release binary")
    package.add_argument("--binary", type=Path, required=True)
    package.add_argument("--target", required=True)
    package.add_argument("--licenses", type=Path, required=True, help="cargo-about HTML report")
    package.add_argument("--output", type=Path, default=Path("target/dist"))
    verify_parser = commands.add_parser("verify", help="verify checksum and smoke-test an extracted archive")
    verify_parser.add_argument("archive", type=Path)
    collect = commands.add_parser("manifest", help="require all artifacts and write SHA256SUMS")
    collect.add_argument("output", type=Path)
    args = parser.parse_args()
    try:
        if args.command == "check":
            version, config = configuration(args.tag)
            print(json.dumps({"version": version, "rust": config["rust"],
                              "macos_minimum": config["macos_deployment_target"],
                              "include": config["targets"]}))
        elif args.command == "pack":
            pack(args.binary, args.target, args.output, args.licenses)
        elif args.command == "verify":
            verify(args.archive)
        else:
            manifest(args.output)
    except (ValueError, subprocess.CalledProcessError) as error:
        parser.exit(1, f"release check failed: {error}\n")


if __name__ == "__main__":
    main()
