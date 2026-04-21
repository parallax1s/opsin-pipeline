"""End-to-end tests for the `cli plm` subcommand (spec §7)."""
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _run_cli(*args: str, cwd: Path | None = None) -> subprocess.CompletedProcess:
    env = {**os.environ, "PYTHONPATH": str(ROOT)}
    return subprocess.run(
        [sys.executable, "-m", "opsin_pipeline.cli", *args],
        cwd=cwd or ROOT,
        capture_output=True,
        text=True,
        env=env,
    )


def _write_scaffold(tmp_path: Path, name: str = "demo", sequence: str = "AFYKWVSG") -> Path:
    payload = {
        "scaffolds": [
            {
                "name": name,
                "family": "Synthetic",
                "sequence": sequence,
                "mutable_positions": [
                    {"position": 1, "allowed": ["F", "Y"]},
                    {"position": 4, "allowed": ["A"]},
                ],
            }
        ]
    }
    path = tmp_path / "scaffolds.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


class PLMCLIMockTests(unittest.TestCase):
    def test_mock_scorer_writes_predictions_json(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            scaffolds_path = _write_scaffold(tmp_path)
            out_path = tmp_path / "plm.json"
            result = _run_cli(
                "plm",
                "--scaffold", "demo",
                "--scaffolds", str(scaffolds_path),
                "--out", str(out_path),
                "--mock",
            )

            self.assertEqual(result.returncode, 0, msg=result.stderr)
            self.assertTrue(out_path.exists())
            data = json.loads(out_path.read_text(encoding="utf-8"))

        self.assertEqual(data["scaffold_name"], "demo")
        self.assertEqual(data["model_id"], "mock")
        # 2 allowed AAs at pos 1 (F, Y — F is silent? sequence is AFYK so pos 1=A, both F and Y real) + 1 at pos 4 (K->A is silent on 4? sequence pos 4=K so A is real)
        # pos 1: A->F, A->Y (both real); pos 4: K->A (real). So 3 predictions.
        self.assertEqual(len(data["predictions"]), 3)
        for p in data["predictions"]:
            self.assertIn(p["to_aa"], ("F", "Y", "A"))
            self.assertGreaterEqual(p["log_likelihood_delta"], -5.0)
            self.assertLessEqual(p["log_likelihood_delta"], 2.0)

    def test_unknown_scaffold_name_exits_nonzero(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            scaffolds_path = _write_scaffold(tmp_path)
            result = _run_cli(
                "plm",
                "--scaffold", "nonexistent",
                "--scaffolds", str(scaffolds_path),
                "--out", str(tmp_path / "plm.json"),
                "--mock",
            )

            self.assertNotEqual(result.returncode, 0)
            self.assertIn("not found", result.stderr + result.stdout)

    def test_non_mock_without_torch_exits_with_hint(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            scaffolds_path = _write_scaffold(tmp_path)
            result = _run_cli(
                "plm",
                "--scaffold", "demo",
                "--scaffolds", str(scaffolds_path),
                "--out", str(tmp_path / "plm.json"),
                # no --mock; also no real backend wired yet
            )

            self.assertNotEqual(result.returncode, 0)
            # Message should mention --mock, esm, or "follow-up" — anything actionable
            msg = (result.stderr + result.stdout).lower()
            self.assertTrue("mock" in msg or "follow-up" in msg or "esm" in msg)

    def test_non_standard_residue_hard_errors(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            # X at position 1 is non-standard -> plm must hard-error
            scaffolds_path = _write_scaffold(tmp_path, sequence="XFYK")
            result = _run_cli(
                "plm",
                "--scaffold", "demo",
                "--scaffolds", str(scaffolds_path),
                "--out", str(tmp_path / "plm.json"),
                "--mock",
            )

            self.assertNotEqual(result.returncode, 0)
            self.assertIn("non-standard", result.stderr + result.stdout)

    def test_sequence_over_1024_exits_with_context_window_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            big_seq = "A" * 1025
            scaffolds_path = _write_scaffold(tmp_path, sequence=big_seq)
            result = _run_cli(
                "plm",
                "--scaffold", "demo",
                "--scaffolds", str(scaffolds_path),
                "--out", str(tmp_path / "plm.json"),
                "--mock",
            )

            self.assertNotEqual(result.returncode, 0)
            self.assertIn("1024", result.stderr + result.stdout)


if __name__ == "__main__":
    unittest.main()
