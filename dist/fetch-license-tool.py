"""Fetch the pinned, checksum-verified cargo-about binary for the Linux plan job."""
import hashlib
import io
from pathlib import Path
import tarfile
import tomllib
import urllib.request

config = tomllib.loads((Path(__file__).parent / "config.toml").read_text())
version = config["license_generator"]
url = (f"https://github.com/EmbarkStudios/cargo-about/releases/download/{version}/"
       f"cargo-about-{version}-x86_64-unknown-linux-musl.tar.gz")
with urllib.request.urlopen(url, timeout=60) as response:
    data = response.read()
if hashlib.sha256(data).hexdigest() != config["license_generator_sha256"]:
    raise ValueError("cargo-about download does not match the pinned checksum")
with tarfile.open(fileobj=io.BytesIO(data)) as archive:
    entry = next(item for item in archive.getmembers()
                 if Path(item.name).name == "cargo-about" and item.isfile())
    destination = Path("target/release-tools/cargo-about")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(archive.extractfile(entry).read())
    destination.chmod(0o755)
