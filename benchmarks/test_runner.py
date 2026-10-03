import gzip
import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location("benchmark_runner", Path(__file__).with_name("run.py"))
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


class ReportTests(unittest.TestCase):
    def test_fasta_roundtrip_and_multiline_input(self):
        with tempfile.TemporaryDirectory() as directory:
            plain = Path(directory) / "input.fasta"
            plain.write_bytes(b">a description\r\nACG\r\nTN\r\n>b\nAAA\n")
            compressed = plain.with_suffix(".gz")
            with gzip.open(compressed, "wb") as output:
                output.write(plain.read_bytes())
            expected = [(b"a description", b"ACGTN"), (b"b", b"AAA")]
            self.assertEqual(runner.fasta(plain), expected)
            self.assertEqual(runner.fasta(compressed), expected)

    def test_equivalence_preserves_headers_and_duplicate_counts(self):
        records = [(b"a", b"ACGT"), (b"b description", b"TT")]
        self.assertEqual(runner.fingerprint(records), runner.fingerprint(records[::-1]))
        self.assertNotEqual(runner.fingerprint(records), runner.fingerprint(records + records[:1]))
        self.assertNotEqual(runner.fingerprint(records), runner.fingerprint([(b"changed", b"ACGT"), records[1]]))

    def test_empty_output_and_one_trial_are_supported(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "empty.fasta"
            path.write_bytes(b"")
            self.assertEqual(runner.fasta(path), [])
        self.assertEqual(runner.distribution([0.5])["median_seconds"], 0.5)
        self.assertEqual(runner.distribution([0.5])["stdev_seconds"], 0)

    def test_macos_metadata_does_not_require_linux_files(self):
        with patch.object(runner.sys, "platform", "darwin"), \
             patch.object(runner, "command", return_value="Apple M2"), \
             patch.object(runner, "optional_text", return_value=None):
            metadata = runner.machine()
        self.assertEqual(metadata["cpu"], "Apple M2")
        self.assertIsNone(metadata["cgroup_cpu_max"])

    def test_failed_rerun_cannot_leave_a_stale_success_report(self):
        with tempfile.TemporaryDirectory() as directory:
            report = Path(directory) / "report.zip"
            report.write_bytes(b"previous successful report")
            with patch.object(runner.sys, "argv", ["run.py", "--output", directory]), \
                 patch.object(runner, "source_files", return_value=[]), \
                 patch.object(runner, "machine", return_value={}), \
                 patch.object(runner, "build", side_effect=ValueError("changed baseline")):
                with self.assertRaisesRegex(ValueError, "changed baseline"):
                    runner.main()
            self.assertFalse(report.exists())


if __name__ == "__main__":
    unittest.main()
