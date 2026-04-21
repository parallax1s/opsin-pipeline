import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from opsin_pipeline.ingest import load_scaffolds
from opsin_pipeline.position_map import (
    UnmappedPocketMapError,
    annotate_with_pocket,
    write_draft_position_map,
)
from opsin_pipeline.structure.ligands import identify_ligands
from opsin_pipeline.structure.pdb import parse_pdb
from opsin_pipeline.structure.pocket import (
    apply_offset_mapping,
    compute_pocket,
    write_pocket_map,
)


FIXTURE_DIR = Path(__file__).parent / "fixtures" / "structures"
SYNTHETIC = FIXTURE_DIR / "synthetic_mini.pdb"


def _make_scaffold_json(tmp_path: Path, name: str = "demo") -> Path:
    payload = {
        "scaffolds": [
            {
                "name": name,
                "family": "Synthetic",
                "sequence": "AFYKWVSG",  # 8 residues, positions 1..8 match fixture
                "target_phenotypes": ["spectral_tuning"],
                "assay_architectures": ["biochemical"],
                "protected_positions": [],
                "mutable_positions": [],
            }
        ]
    }
    path = tmp_path / f"{name}_scaffolds.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def _make_mapped_pocket(tmp_path: Path, scaffold_name: str = "demo") -> Path:
    atoms = parse_pdb(SYNTHETIC)
    pocket = compute_pocket(
        atoms,
        identify_ligands(atoms),
        scaffold_name=scaffold_name,
        pdb_path=str(SYNTHETIC),
    )
    mapped = apply_offset_mapping(pocket, offset=0)
    return write_pocket_map(mapped, tmp_path / "pocket.json")


def _make_unmapped_pocket(tmp_path: Path, scaffold_name: str = "demo") -> Path:
    atoms = parse_pdb(SYNTHETIC)
    pocket = compute_pocket(
        atoms,
        identify_ligands(atoms),
        scaffold_name=scaffold_name,
        pdb_path=str(SYNTHETIC),
    )
    return write_pocket_map(pocket, tmp_path / "pocket.json")


class PocketAnnotateTests(unittest.TestCase):
    def test_mapped_pocket_annotates_scaffold_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            scaffolds_path = _make_scaffold_json(tmp_path)
            scaffolds = load_scaffolds(scaffolds_path)
            draft_path = write_draft_position_map(scaffolds, tmp_path / "draft.csv")
            pocket_path = _make_mapped_pocket(tmp_path)

            out_path = annotate_with_pocket(draft_path, pocket_path, tmp_path / "annotated.csv")

            with out_path.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle))

        by_index = {int(row["seq_index"]): row for row in rows if row["scaffold"] == "demo"}
        # position 1 (ALA): distance 2.0 A, band "strong", region gains retinal_pocket
        self.assertEqual(by_index[1]["distance_to_retinal_A"], "2.00")
        self.assertEqual(by_index[1]["pocket_band"], "strong")
        self.assertIn("retinal_pocket", by_index[1]["region"])
        # position 6 (VAL): 6.0 A, band "none", region is NOT tagged retinal_pocket
        self.assertEqual(by_index[6]["pocket_band"], "none")
        self.assertNotIn("retinal_pocket", by_index[6]["region"])
        # position 7 (SER): outside cutoff entirely, no annotation
        self.assertEqual(by_index[7]["distance_to_retinal_A"], "")
        self.assertEqual(by_index[7]["pocket_band"], "")

    def test_mutable_and_protected_are_not_touched(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            scaffolds_path = _make_scaffold_json(tmp_path)
            scaffolds = load_scaffolds(scaffolds_path)
            draft_path = write_draft_position_map(scaffolds, tmp_path / "draft.csv")
            pocket_path = _make_mapped_pocket(tmp_path)

            out_path = annotate_with_pocket(draft_path, pocket_path, tmp_path / "annotated.csv")

            with out_path.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle))

        self.assertTrue(all(row["mutable"] == "false" for row in rows))
        self.assertTrue(all(row["protected"] == "false" for row in rows))

    def test_unmapped_pocket_raises_and_cli_exits_nonzero(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            scaffolds_path = _make_scaffold_json(tmp_path)
            scaffolds = load_scaffolds(scaffolds_path)
            draft_path = write_draft_position_map(scaffolds, tmp_path / "draft.csv")
            pocket_path = _make_unmapped_pocket(tmp_path)

            with self.assertRaises(UnmappedPocketMapError):
                annotate_with_pocket(draft_path, pocket_path, tmp_path / "annotated.csv")

            # CLI variant of the same guard
            out_path = tmp_path / "cli_annotated.csv"
            result = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "opsin_pipeline.cli",
                    "pocket-annotate",
                    "--position-map",
                    str(draft_path),
                    "--pocket-map",
                    str(pocket_path),
                    "--out",
                    str(out_path),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
                env={**__import__("os").environ, "PYTHONPATH": str(ROOT)},
            )

            self.assertNotEqual(result.returncode, 0)
            self.assertFalse(out_path.exists())
            self.assertIn("Re-run", result.stderr + result.stdout)


if __name__ == "__main__":
    unittest.main()
