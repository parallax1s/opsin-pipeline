import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

FIXTURE_DIR = Path(__file__).parent / "fixtures" / "structures"
SYNTHETIC = FIXTURE_DIR / "synthetic_mini.pdb"
SYNTHETIC_LYR = FIXTURE_DIR / "synthetic_lyr.pdb"


def _run_cli(*args: str) -> subprocess.CompletedProcess:
    env = {"PYTHONPATH": str(ROOT)}
    return subprocess.run(
        [sys.executable, "-m", "opsin_pipeline.cli", *args],
        cwd=ROOT,
        capture_output=True,
        text=True,
        env={**__import__("os").environ, **env},
    )


class PocketCLITests(unittest.TestCase):
    def test_pocket_subcommand_writes_json_with_expected_bands(self):
        with tempfile.TemporaryDirectory() as tmp:
            out_path = Path(tmp) / "synthetic_pocket.json"
            result = _run_cli(
                "pocket",
                "--scaffold", "synthetic",
                "--pdb", str(SYNTHETIC),
                "--out", str(out_path),
            )

            self.assertEqual(result.returncode, 0, msg=result.stderr)
            self.assertTrue(out_path.exists())
            data = json.loads(out_path.read_text(encoding="utf-8"))

        self.assertEqual(data["scaffold_name"], "synthetic")
        self.assertEqual(len(data["ligand_matches"]), 1)
        # 3 strong + 2 medium + 1 none within the default 6.0 A cutoff
        bands = [r["band"] for r in data["pocket_residues"]]
        self.assertEqual(bands.count("strong"), 3)
        self.assertEqual(bands.count("medium"), 2)
        self.assertEqual(bands.count("none"), 1)

    def test_pocket_with_offset_sets_seq_index(self):
        with tempfile.TemporaryDirectory() as tmp:
            out_path = Path(tmp) / "pocket.json"
            result = _run_cli(
                "pocket",
                "--scaffold", "synthetic",
                "--pdb", str(SYNTHETIC),
                "--pdb-offset", "100",
                "--out", str(out_path),
            )

            self.assertEqual(result.returncode, 0, msg=result.stderr)
            data = json.loads(out_path.read_text(encoding="utf-8"))

        for residue in data["pocket_residues"]:
            self.assertEqual(residue["seq_index"], residue["pdb_resnum"] + 100)

    def test_pocket_lyr_schiff_base_emitted(self):
        with tempfile.TemporaryDirectory() as tmp:
            out_path = Path(tmp) / "pocket.json"
            result = _run_cli(
                "pocket",
                "--scaffold", "br_synth",
                "--pdb", str(SYNTHETIC_LYR),
                "--out", str(out_path),
            )

            self.assertEqual(result.returncode, 0, msg=result.stderr)
            data = json.loads(out_path.read_text(encoding="utf-8"))

        schiff = [r for r in data["pocket_residues"] if r.get("role") == "schiff_base_linkage"]
        self.assertEqual(len(schiff), 1)
        self.assertEqual(schiff[0]["pdb_resnum"], 216)
        self.assertEqual(schiff[0]["band"], "strong")

    def test_pocket_missing_ligand_exits_nonzero(self):
        with tempfile.TemporaryDirectory() as tmp:
            # a PDB that has no whitelisted retinal-like HETATM
            bad = Path(tmp) / "no_retinal.pdb"
            bad.write_text(
                "HETATM    1  O1  GOL A 201      50.000  50.000  50.000  1.00 40.00           O\n",
                encoding="utf-8",
            )
            out_path = Path(tmp) / "pocket.json"
            result = _run_cli(
                "pocket",
                "--scaffold", "demo",
                "--pdb", str(bad),
                "--out", str(out_path),
            )

            self.assertNotEqual(result.returncode, 0)
            self.assertFalse(out_path.exists())


if __name__ == "__main__":
    unittest.main()
