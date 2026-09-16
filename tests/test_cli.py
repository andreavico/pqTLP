"""Exercise generation and solving in separate Sage processes."""

import subprocess
import sys
import unittest
from pathlib import Path


class CommandLineTests(unittest.TestCase):
    def test_hex_secret_survives_serialization(self):
        cli = Path(__file__).resolve().parents[1] / "TLP.sage"
        secret = "0001020300"
        generated = subprocess.run(
            [sys.executable, str(cli), "--generate", "0.01", "--secret-hex", secret],
            capture_output=True,
            text=True,
            timeout=120,
        )
        self.assertEqual(generated.returncode, 0, generated.stderr)
        self.assertTrue(generated.stdout.startswith("PQTLP\n"))
        solved = subprocess.run(
            [sys.executable, str(cli), "--solve", "-"],
            input=generated.stdout,
            capture_output=True,
            text=True,
            timeout=120,
        )
        self.assertEqual(solved.returncode, 0, solved.stderr)
        self.assertEqual(solved.stdout.strip(), secret)


if __name__ == "__main__":
    unittest.main()
