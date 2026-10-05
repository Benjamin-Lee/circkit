"""Release failures that must be caught before uploading artifacts."""
import hashlib
import io
from pathlib import Path
import tarfile
import tempfile
import unittest
from unittest.mock import patch

import release


class ReleaseTests(unittest.TestCase):
    def test_tag_and_dependency_versions_must_match(self):
        version, _ = release.configuration()
        self.assertEqual(release.configuration(f"v{version}")[0], version)
        with self.assertRaisesRegex(ValueError, "does not match"):
            release.configuration("v999.0.0")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "lib").mkdir()
            (root / "dist").mkdir()
            for filename in ("Cargo.toml", "lib/Cargo.toml", "dist/config.toml", "LICENSE", "lib/LICENSE"):
                (root / filename).write_bytes((release.ROOT / filename).read_bytes())
            cli = root / "Cargo.toml"
            cli.write_text(cli.read_text().replace(f'path = "lib", version = "{version}"', 'path = "lib", version = "999.0.0"'))
            with patch.object(release, "ROOT", root), self.assertRaisesRegex(ValueError, "dependency"):
                release.configuration()

    def test_incomplete_release_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaisesRegex(ValueError, "missing release artifact"):
                release.manifest(Path(temporary))
            self.assertFalse((Path(temporary) / "SHA256SUMS").exists())

    def test_manifest_verifies_existing_checksums(self):
        version, config = release.configuration()
        names = [f"circkit-{version}-{item['target']}.tar.gz" for item in config["targets"]]
        names += [f"{crate}-{version}.crate" for crate in ("circkit", "circkit-cli")]
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            for name in names:
                (output / name).write_bytes(name.encode())
            release.manifest(output)
            recorded = {line.split("  ", 1)[1] for line in (output / "SHA256SUMS").read_text().splitlines()}
            self.assertEqual(recorded, set(names))
            (output / names[-1]).write_bytes(b"changed")
            with self.assertRaisesRegex(ValueError, "checksum mismatch"):
                release.manifest(output)

    def test_corrupt_archive_is_rejected_before_execution(self):
        with tempfile.TemporaryDirectory() as temporary:
            archive = Path(temporary) / "circkit.tar.gz"
            archive.write_bytes(b"changed")
            archive.with_suffix(".gz.sha256").write_text(f"{'0' * 64}  {archive.name}\n")
            with self.assertRaisesRegex(ValueError, "checksum mismatch"):
                release.verify(archive)

    def test_archive_cannot_extract_outside_temporary_directory(self):
        with tempfile.TemporaryDirectory() as temporary:
            archive = Path(temporary) / "circkit.tar.gz"
            with tarfile.open(archive, "w:gz") as tar:
                info = tarfile.TarInfo("../escape")
                info.size = 1
                tar.addfile(info, io.BytesIO(b"x"))
            checksum = hashlib.sha256(archive.read_bytes()).hexdigest()
            archive.with_suffix(".gz.sha256").write_text(f"{checksum}  {archive.name}\n")
            with self.assertRaises(tarfile.OutsideDestinationError):
                release.verify(archive)


if __name__ == "__main__":
    unittest.main()
