#!/usr/bin/env python3
"""Native Linux/macOS benchmarks using only Python's standard library."""
import argparse
import csv
import datetime
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import platform
import random
import shutil
import statistics
import subprocess
import sys
import tarfile
import time
import zipfile

REPO = Path(__file__).resolve().parents[1]
BASELINE = "deb0b407832cfc276ccf351f38f9201b90f90d34"  # Merged PR #2.
SEED = 20261003
DEFAULT_CUTOFF = 32768
DEFAULT_CHUNK = 256


def command(args, cwd=REPO, **kwargs):
    return subprocess.check_output(args, cwd=cwd, text=True, **kwargs).strip()


def sha256(path):
    result = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def numbers(value, zero=False):
    try:
        result = list(dict.fromkeys(int(x) for x in value.split(",")))
        if not result or any(x < (0 if zero else 1) for x in result):
            raise ValueError()
        return result
    except ValueError:
        raise argparse.ArgumentTypeError("Expected comma-separated positive integers" if not zero
                                         else "Expected comma-separated nonnegative integers")


def positive(value):
    result = int(value)
    if result < 1:
        raise argparse.ArgumentTypeError("Must be at least one")
    return result


def optional_text(path):
    try:
        return Path(path).read_text().strip()
    except OSError:
        return None


def machine():
    cpu = platform.processor()
    features = None
    if sys.platform == "darwin":
        try:
            cpu = command(["sysctl", "-n", "machdep.cpu.brand_string"])
        except subprocess.CalledProcessError:
            cpu = command(["sysctl", "-n", "hw.model"])
    elif sys.platform.startswith("linux"):
        lines = (optional_text("/proc/cpuinfo") or "").splitlines()
        cpu = next((s.split(":", 1)[1].strip() for s in lines if s.startswith("model name")), cpu)
        if not cpu:
            cpu = "; ".join(s.strip() for s in lines if s.startswith(("CPU implementer", "CPU part")))
        features = next((s.split(":", 1)[1].split() for s in lines
                         if s.startswith(("flags", "Features"))), None)
    affinity = sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None
    return {
        "platform": platform.platform(), "architecture": platform.machine(), "cpu": cpu,
        "logical_cpus": os.cpu_count(), "affinity": affinity, "cpu_features": features,
        "cgroup_cpu_max": optional_text("/sys/fs/cgroup/cpu.max"),
        "cgroup_v1_quota": optional_text("/sys/fs/cgroup/cpu/cpu.cfs_quota_us"),
        "cgroup_v1_period": optional_text("/sys/fs/cgroup/cpu/cpu.cfs_period_us"),
        "python": platform.python_version(),
    }


def source_files():
    paths = subprocess.check_output(
        ["git", "ls-files", "-z", "--cached", "--others", "--exclude-standard"], cwd=REPO
    ).decode().split("\0")
    return sorted({p for p in paths if p and (REPO / p).is_file()})


def source_hash(files):
    result = hashlib.sha256()
    for name in files:
        result.update(name.encode() + b"\0" + (REPO / name).read_bytes() + b"\0")
    return result.hexdigest()


def dependency_versions(lock):
    versions = {}
    for block in lock.read_text().split("[[package]]")[1:]:
        fields = {}
        for line in block.splitlines():
            if line.startswith(("name = ", "version = ")):
                key, value = line.split(" = ", 1)
                fields[key] = value.strip('"')
        if fields.get("name") in ("bio", "memchr", "flate2", "zlib-rs", "triple_accel"):
            versions.setdefault(fields["name"], []).append(fields["version"])
    return versions


def build(args, output, files):
    baseline = command(["git", "rev-parse", args.baseline + "^{commit}"])
    rustc = os.environ.get("RUSTC") or (str(Path(args.cargo).with_name("rustc"))
                                      if Path(args.cargo).is_file() else "rustc")
    compiler = command([rustc, "-vV"])
    host = next(s.split(": ", 1)[1] for s in compiler.splitlines() if s.startswith("host:"))
    target = os.environ.get("CARGO_BUILD_TARGET")
    if target and target != host:
        raise ValueError("Benchmarks require a native build; CARGO_BUILD_TARGET differs from rustc's host")
    signature = {
        "baseline_commit": baseline, "candidate_head": command(["git", "rev-parse", "HEAD"]),
        "candidate_source_sha256": source_hash(files), "rustc": compiler,
        "cargo": command([args.cargo, "--version"]),
        "build_environment": {key: os.environ.get(key) for key in (
            "RUSTFLAGS", "CARGO_ENCODED_RUSTFLAGS", "CARGO_BUILD_TARGET",
            "RUSTC", "RUSTC_WRAPPER", "RUSTUP_TOOLCHAIN", "CARGO_BUILD_RUSTFLAGS",
            "CARGO_PROFILE_RELEASE_LTO", "CARGO_PROFILE_RELEASE_CODEGEN_UNITS",
            "CARGO_PROFILE_RELEASE_OPT_LEVEL", "CARGO_PROFILE_RELEASE_DEBUG")
            + tuple(k for k in os.environ if k.startswith("CARGO_TARGET_") and k.endswith("_RUSTFLAGS"))},
        "scalar_variant": not args.no_scalar,
    }
    info_path = output / "builds.json"
    if args.skip_build:
        info = json.loads(info_path.read_text())
        if info["signature"] != signature:
            raise ValueError("Source, baseline, or toolchain changed; rerun without --skip-build")
        for stage, data in info["stages"].items():
            for name, digest in data["binaries"].items():
                if sha256(output / "binaries" / name) != digest:
                    raise ValueError("Benchmark binary changed: " + stage)
        return info
    snapshots = output / "build" / "sources"
    snapshots.mkdir(parents=True, exist_ok=True)
    baseline_source = snapshots / "baseline"
    candidate_source = snapshots / "optimized"
    for path in (baseline_source, candidate_source):
        if path.exists():
            shutil.rmtree(path)
        path.mkdir()
    archive = subprocess.check_output(["git", "archive", baseline], cwd=REPO)
    with tarfile.open(fileobj=io.BytesIO(archive)) as source:
        # Python 3.9-compatible extraction, limited to ordinary repository files.
        for member in source.getmembers():
            destination = (baseline_source / member.name).resolve()
            if baseline_source.resolve() not in destination.parents or not (member.isfile() or member.isdir()):
                raise ValueError("Unsupported archive member: " + member.name)
        if sys.version_info >= (3, 12):
            source.extractall(baseline_source, filter="data")
        else:
            source.extractall(baseline_source)
    for name in files:
        destination = candidate_source / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(REPO / name, destination)
    (baseline_source / "examples").mkdir(exist_ok=True)
    shutil.copytree(REPO / "benchmarks", baseline_source / "benchmarks", dirs_exist_ok=True)
    driver = (REPO / "examples/performance.rs").read_text()
    adapter = '#[path = "../benchmarks/optimized.rs"]'
    if driver.count(adapter) != 1:
        raise ValueError("Expected a single benchmark adapter path")
    (baseline_source / "examples/performance.rs").write_text(
        driver.replace(adapter, '#[path = "../benchmarks/baseline.rs"]'))
    binary_dir = output / "binaries"
    binary_dir.mkdir(exist_ok=True)
    info = {"signature": signature, "default_profile": {"release": True, "lto": "fat", "codegen_units": 1},
            "candidate_dirty": command(["git", "status", "--porcelain"]),
            "instrumentation": "Identical driver; baseline adapter uses only the PR #2 public API",
            "stages": {}}
    stages = [("baseline", baseline_source, []), ("optimized", candidate_source, [])]
    if not args.no_scalar:
        stages.append(("optimized-scalar", candidate_source, ["--features", "scalar-normalization"]))
    for stage, source, features in stages:
        # Cargo keeps separate baseline/candidate targets; scalar shares candidate dependencies.
        target_dir = output / "build" / ("target-baseline" if stage == "baseline" else "target-optimized")
        env = dict(os.environ, CARGO_TARGET_DIR=str(target_dir))
        build_command = [args.cargo, "build", "--verbose", "--locked", "--release", "--bin", "circkit",
                         "--example", "performance", "--jobs", str(args.build_jobs)] + features
        if args.offline:
            build_command.append("--offline")
        print("Building " + stage + " (log: build-" + stage + ".log)", flush=True)
        with (output / ("build-" + stage + ".log")).open("w") as log:
            subprocess.run(build_command, cwd=source, env=env, stdout=log, stderr=log, check=True)
        release = target_dir / target / "release" if target else target_dir / "release"
        for name, path in ((stage, release / "circkit"), (stage + "-micro", release / "examples/performance")):
            shutil.copy2(path, binary_dir / name)
        info["stages"][stage] = {
            "command": build_command, "dependencies": dependency_versions(source / "Cargo.lock"),
            "lock_sha256": sha256(source / "Cargo.lock"),
            "binaries": {name: sha256(binary_dir / name) for name in (stage, stage + "-micro")},
        }
    write_json(info_path, info)
    return info


def fasta(path):
    opener = gzip.open if path.suffix == ".gz" else open
    records, header, sequence = [], None, []
    with opener(path, "rb") as source:
        for raw in source:
            line = raw.rstrip(b"\r\n")
            if line.startswith(b">"):
                if header is not None:
                    records.append((header, b"".join(sequence)))
                header, sequence = line[1:], []
            elif header is not None:
                sequence.append(line)
            elif line:
                raise ValueError("Expected FASTA input: " + str(path))
    if header is not None:
        records.append((header, b"".join(sequence)))
    return records


def fingerprint(records):
    digest = hashlib.sha256()
    for header, sequence in sorted(records):
        for field in (header, sequence):
            digest.update(len(field).to_bytes(8, "little"))
            digest.update(field)
    return digest.hexdigest()


def fixtures(output, suite):
    data = output / "datasets"
    data.mkdir(exist_ok=True)
    count, medium, dense, false_seeds = {
        "smoke": (200, 3, 3000, 100),
        "quick": (10000, 25, 30000, 1000),
        "full": (100000, 1000, 300000, 10000),
    }[suite]
    rng = random.Random(SEED)
    def write(name, n, length, kind="random", wrapped=False):
        path = data / (name + ".fasta")
        with path.open("w", newline="\n") as out:
            for index in range(n):
                if kind == "dense":
                    sequence = ("ATGAAATAA" * ((length + 8) // 9))[:length]
                elif kind == "dimer":
                    half = "".join(rng.choices("ACGT", k=length // 2))
                    sequence = half + half
                elif kind == "false-seeds":
                    sequence = "C" + "AAAAAAAAAAG" * false_seeds + "AAAAAAAAAA"
                else:
                    sequence = "".join(rng.choices("ACGT", k=length))
                out.write(">sequence_" + str(index) + " description\n")
                out.write("\n".join(sequence[i:i + 60] for i in range(0, len(sequence), 60))
                          if wrapped else sequence)
                out.write("\n")
        return path
    paths = {
        "short": write("short", count, 150),
        "dimers": write("dimers", count // 2, 300, "dimer"),
        "medium": write("medium", medium, 10000, wrapped=True),
        "long": write("long", 3, 208399),
        "dense": write("dense", 1, dense, "dense"),
        "false-seeds": write("false-seeds", 1, 0, "false-seeds"),
    }
    paths["gzip"] = data / "short.fasta.gz"
    with paths["gzip"].open("wb") as dest:
        with gzip.GzipFile(fileobj=dest, mode="wb", mtime=0) as compressed:
            compressed.write(paths["short"].read_bytes())
    return paths


def distribution(values):
    quartiles = statistics.quantiles(values, n=4, method="inclusive") if len(values) > 1 else values * 3
    return {"median_seconds": statistics.median(values), "q1_seconds": quartiles[0],
            "q3_seconds": quartiles[2], "min_seconds": min(values), "max_seconds": max(values),
            "stdev_seconds": statistics.stdev(values) if len(values) > 1 else 0,
            "runs": len(values)}


class Runner:
    def __init__(self, args, output):
        self.args, self.output = args, output
        self.rng = random.Random(SEED)
        self.rows, self.raw, self.validation = [], [], []
        self.scratch = output / "scratch"
        self.scratch.mkdir(exist_ok=True)

    def save(self):
        write_json(self.output / "results.json", self.rows)
        write_json(self.output / "raw.json", self.raw)
        write_json(self.output / "validation.json", self.validation)
        columns = ["category", "case", "stage", "reference", "speedup", "median_seconds",
                   "q1_seconds", "q3_seconds", "stdev_seconds", "min_seconds", "max_seconds",
                   "runs", "settings"]
        with (self.output / "results.csv").open("w", newline="") as out:
            writer = csv.DictWriter(out, fieldnames=columns)
            writer.writeheader()
            for row in self.rows:
                writer.writerow(dict(row, settings=json.dumps(row["settings"], sort_keys=True)))

    def measure(self, category, case, stages, reference, settings, sample):
        samples = {stage: [] for stage in stages}
        for _ in range(self.args.warmups):
            for stage in stages:
                sample(stage)
        for trial in range(self.args.trials):
            order = list(stages)
            self.rng.shuffle(order)
            for stage in order:
                seconds, extra = sample(stage)
                samples[stage].append(seconds)
                self.raw.append({"category": category, "case": case, "stage": stage,
                                 "trial": trial, "seconds": seconds, "settings": settings, **extra})
        ref = statistics.median(samples[reference])
        for stage, values in samples.items():
            stats = distribution(values)
            self.rows.append({"category": category, "case": case, "stage": stage,
                              "reference": reference, "speedup": ref / stats["median_seconds"],
                              "settings": settings, **stats})
        self.save()
        ratios = ", ".join(stage + "=" + format(ref / statistics.median(v), ".2f") + "x"
                           for stage, v in samples.items() if stage != reference)
        print(case + ": " + ratios, flush=True)

    def cli(self, name, operation, dataset, threads, candidate_extra=(), scalar=False, compressed=False):
        stages = ["baseline", "optimized"]
        if scalar and not self.args.no_scalar:
            stages.append("optimized-scalar")
        settings = {"operation": operation, "input": str(dataset), "threads": threads,
                    "candidate_args": list(candidate_extra), "gzip_output": compressed}
        common = []
        if operation == "monomerize":
            common = ["--keep-all"]
            if "approx" in name:
                common += ["--max-mismatch", "2"]
        if operation == "orfs" and "dense" in name:
            common = ["--min-length", "0", "--include-stop"]
        if "sensitive" in name:
            common += ["--sensitive"]
        settings["common_args"] = common
        def cmd(stage, output):
            args = [str(self.output / "binaries" / stage), operation, str(dataset),
                    "--threads", str(threads), "-o", str(output)] + common
            return args + (list(candidate_extra) if stage != "baseline" else [])
        digests, tables = {}, {}
        for stage in stages:
            path = self.scratch / (stage + (".fasta.gz" if compressed else ".fasta"))
            table = self.scratch / (stage + ".tsv")
            args = cmd(stage, path)
            if operation != "canonicalize":
                args += ["--table", str(table)]
            subprocess.run(args, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
                           timeout=self.args.timeout, check=True)
            digests[stage] = fingerprint(fasta(path))
            if operation != "canonicalize":
                with table.open(newline="") as source:
                    rows = list(csv.reader(source, delimiter="\t"))
                tables[stage] = {"header": rows[:1], "rows": sorted(rows[1:])}
        if len(set(digests.values())) != 1 or any(t != tables[stages[0]] for t in tables.values()):
            raise ValueError("FASTA or metadata differs for " + name)
        self.validation.append({"case": name, "type": "cli", "settings": settings,
                                "fasta_record_hashes": digests, "metadata_identical": True})
        def sample(stage):
            dest = self.scratch / (stage + ".timed.gz") if compressed else os.devnull
            start = time.perf_counter_ns()
            subprocess.run(cmd(stage, dest), stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
                           timeout=self.args.timeout, check=True)
            return (time.perf_counter_ns() - start) / 1e9, {}
        self.measure("cli", name, stages, "baseline", settings, sample)

    def micro(self, name, operation, length, pattern="random", cutoffs=None, chunks=None, mismatch=0):
        configurations = {"baseline": ("baseline", DEFAULT_CUTOFF, DEFAULT_CHUNK)}
        if cutoffs is not None:
            configurations.update({str(c): ("optimized", c, DEFAULT_CHUNK) for c in cutoffs})
        elif chunks is not None:
            configurations.update({str(c): ("optimized", DEFAULT_CUTOFF, c) for c in chunks})
        else:
            configurations["optimized"] = ("optimized", DEFAULT_CUTOFF, DEFAULT_CHUNK)
        iterations, checksums = {}, {}
        def run(stage, loops):
            binary, cutoff, chunk = configurations[stage]
            result = subprocess.check_output(
                [str(self.output / "binaries" / (binary + "-micro")), operation,
                 str(length), str(loops), pattern, str(cutoff), str(chunk), str(mismatch)],
                timeout=self.args.timeout, text=True)
            return json.loads(result)
        for stage in configurations:
            result = run(stage, 1)
            checksums[stage] = result["checksum"]
            target = 0.02 if self.args.suite == "smoke" else 0.10
            iterations[stage] = max(8, min(2000000, int(target * 1e9 / max(1, result["elapsed_ns"]))))
        if len(set(checksums.values())) != 1:
            raise ValueError("Library output differs for " + name)
        settings = {"operation": operation, "length": length, "pattern": pattern, "max_mismatch": mismatch,
                    "configurations": configurations, "iterations": iterations}
        self.validation.append({"case": name, "type": "library", "settings": settings,
                                "checksums": checksums, "identical": True})
        def sample(stage):
            result = run(stage, iterations[stage])
            return result["elapsed_ns"] / result["iterations"] / 1e9, {"iterations": result["iterations"]}
        self.measure("library", name, configurations, "baseline", settings, sample)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=REPO / "target/benchmarks")
    parser.add_argument("--baseline", default=BASELINE)
    parser.add_argument("--suite", choices=["smoke", "quick", "full"], default="quick")
    parser.add_argument("--trials", type=positive, default=5)
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--threads", type=numbers, default=[1, min(4, os.cpu_count() or 1)])
    parser.add_argument("--rotation-cutoffs", type=lambda s: numbers(s, zero=True), default=[0, 32768, 1048576])
    parser.add_argument("--mismatch-chunks", type=numbers, default=[64, 256, 1024])
    parser.add_argument("--queue-depths", type=numbers, default=[2, 8, 64])
    parser.add_argument("--cpus", type=lambda s: numbers(s, zero=True), help="Linux affinity; leave unset on macOS")
    parser.add_argument("--input", type=Path, action="append", default=[], help="Additional local FASTA (plain or gzip)")
    parser.add_argument("--note", action="append", default=[], help="Machine context, e.g. plugged in or power mode")
    parser.add_argument("--cargo", default=shutil.which("cargo") or "cargo")
    parser.add_argument("--build-jobs", type=positive, default=2)
    parser.add_argument("--timeout", type=positive, default=120)
    parser.add_argument("--offline", action="store_true")
    parser.add_argument("--skip-build", action="store_true")
    parser.add_argument("--no-scalar", action="store_true")
    parser.add_argument("--no-tuning", action="store_true")
    args = parser.parse_args()
    if args.warmups < 0:
        parser.error("--warmups must be nonnegative")
    args.output = args.output.expanduser().resolve()
    # Output must be ignored by Git or outside the checkout to keep snapshots bounded.
    if args.output == REPO or REPO in args.output.parents:
        relative = str(args.output.relative_to(REPO))
        if subprocess.run(["git", "check-ignore", "-q", relative + "/probe"], cwd=REPO).returncode:
            parser.error("Choose an output under target/ or outside the checkout")
    if args.cpus:
        if not hasattr(os, "sched_setaffinity"):
            parser.error("--cpus is unavailable on this OS; omit it to use the OS scheduler")
        os.sched_setaffinity(0, args.cpus)
    args.output.mkdir(parents=True, exist_ok=True)
    files = source_files()
    environment = {"schema_version": 1, "started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
                   "machine": machine(), "settings": vars(args).copy(), "seed": SEED,
                   "candidate_defaults": {"rotation_cutoff": DEFAULT_CUTOFF,
                                          "mismatch_chunk_size": DEFAULT_CHUNK,
                                          "execution": "auto", "queue_depth": "64; orfs max(2, 2*threads)",
                                          "gzip_level": 6},
                   "method": {"cli": "Wall time including startup, parsing and writing; warm caches",
                              "library": "Time per call excluding fixture setup, checksum and startup",
                              "validation": "FASTA record multiset including headers, plus metadata row multiset",
                              "order": "Seeded interleaving of stages; independent calibration for library stages",
                              "memory": "No memory measurements"}}
    environment["settings"] = {k: str(v) if isinstance(v, Path) else
                               [str(x) if isinstance(x, Path) else x for x in v] if isinstance(v, list) else v
                               for k, v in environment["settings"].items()}
    write_json(args.output / "environment.json", environment)
    build(args, args.output, files)
    paths = fixtures(args.output, args.suite)
    for index, path in enumerate(args.input):
        paths["custom-" + str(index)] = path.expanduser().resolve()
    manifest = []
    for name, path in paths.items():
        records = fasta(path)
        lengths = [len(s) for _, s in records]
        manifest.append({"name": name, "path": str(path), "sha256": sha256(path),
                         "bytes": path.stat().st_size, "records": len(records),
                         "length_min": min(lengths, default=0), "length_max": max(lengths, default=0),
                         "length_mean": statistics.mean(lengths) if lengths else 0})
    write_json(args.output / "inputs.json", manifest)
    runner = Runner(args, args.output)
    for op, dataset, suffix in [
        ("canonicalize", "short", ""), ("canonicalize", "long", ""), ("canonicalize", "medium", ""),
        ("uniq", "short", ""), ("monomerize", "short", ""), ("monomerize", "dimers", ""),
        ("monomerize", "medium", "-sensitive"), ("monomerize", "false-seeds", "-approx-stress"),
        ("orfs", "short", ""), ("orfs", "medium", ""), ("orfs", "dense", "-stress"),
        ("canonicalize", "gzip", ""),
    ]:
        runner.cli(op + "-" + dataset + suffix, op, paths[dataset], 1,
                   scalar=dataset == "short", compressed=dataset == "gzip")
    for key in paths:
        if key.startswith("custom-"):
            for op in ["canonicalize", "monomerize", "orfs"]:
                runner.cli(op + "-" + key, op, paths[key], 1)
    for threads in sorted(set(args.threads) - {1}):
        for op, key in [("canonicalize", "short"), ("monomerize", "dimers"), ("orfs", "short")]:
            runner.cli(op + "-" + key + "-threads-" + str(threads), op, paths[key], threads)
    for op, length, pattern in [
        ("lmsr-index", 150, "random"), ("canonicalize", 150, "random"),
        ("lmsr-index", 208399, "random"), ("monomerize", 300, "dimer"),
        ("find-orfs", 150, "random"), ("find-orfs", 3000 if args.suite == "smoke" else 30000, "dense"),
        ("orf-indices", 3000 if args.suite == "smoke" else 30000, "dense"),
        ("extract-orf", 10000, "random"), ("normalize", 150, "random"),
    ]:
        runner.micro(op + "-" + str(length) + "-" + pattern, op, length, pattern)
    if not args.no_tuning:
        lengths = {"smoke": [150, 32769], "quick": [150, 4096, 32768, 65536],
                   "full": [150, 1000, 4096, 8192, 32768, 32769, 65536, 208399]}[args.suite]
        for length in lengths:
            for pattern in ["random", "homopolymer"]:
                runner.micro("rotation-cutoff-" + str(length) + "-" + pattern, "lmsr-index",
                             length, pattern, cutoffs=args.rotation_cutoffs)
        for pattern, length in [("dimer", 300), ("false-seeds", 11011)]:
            runner.micro("mismatch-chunk-" + pattern, "monomerize", length, pattern,
                         chunks=args.mismatch_chunks, mismatch=2)
        for op, key in [("canonicalize", "short"), ("monomerize", "dimers"), ("orfs", "short")]:
            for mode in ["serial", "pipeline"]:
                runner.cli(op + "-execution-" + mode, op, paths[key], 1, ["--execution", mode])
        queue_threads = max(args.threads)
        for depth in args.queue_depths:
            runner.cli("orfs-queue-" + str(depth), "orfs", paths["medium"], queue_threads,
                       ["--execution", "pipeline", "--queue-depth", str(depth)])
    environment["completed_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    write_json(args.output / "environment.json", environment)
    report_files = ["environment.json", "builds.json", "inputs.json", "results.csv",
                    "results.json", "raw.json", "validation.json"]
    report_files += [p.name for p in args.output.glob("build-*.log")]
    with zipfile.ZipFile(args.output / "report.zip", "w", zipfile.ZIP_DEFLATED) as archive:
        for name in report_files:
            archive.write(args.output / name, name)
    print("All outputs matched. Share " + str(args.output / "report.zip"), flush=True)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        print("Benchmark failed: " + str(error), file=sys.stderr)
        if isinstance(error, subprocess.CalledProcessError) and error.stderr:
            print(error.stderr.decode(errors="replace") if isinstance(error.stderr, bytes) else error.stderr,
                  file=sys.stderr)
        sys.exit(1)
