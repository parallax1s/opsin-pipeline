import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from opsin_pipeline.generate import generate_single_mutants
from opsin_pipeline.ingest import load_scaffolds
from opsin_pipeline.position_map import (
    apply_position_map_to_scaffolds,
    write_draft_position_map,
)


class PositionMapTests(unittest.TestCase):
    def write_scaffold_json(self, directory: Path) -> Path:
        payload = {
            "scaffolds": [
                {
                    "name": "miniRhGC",
                    "family": "RhGC",
                    "sequence": "MKTAYK",
                    "target_phenotypes": ["spectral_tuning", "cGMP_activity"],
                    "assay_architectures": ["biochemical", "growth_selection"],
                    "protected_positions": [],
                    "mutable_positions": [],
                },
                {
                    "name": "miniChR",
                    "family": "ChR",
                    "sequence": "MKTAFK",
                    "target_phenotypes": ["spectral_tuning", "photocurrent"],
                    "assay_architectures": ["growth_selection"],
                    "protected_positions": [],
                    "mutable_positions": [],
                },
            ]
        }
        path = directory / "scaffolds.json"
        path.write_text(json.dumps(payload), encoding="utf-8")
        return path

    def test_write_draft_position_map_contains_one_row_per_residue(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            scaffolds = load_scaffolds(self.write_scaffold_json(tmp_path))
            out_path = write_draft_position_map(scaffolds, tmp_path / "position_map.csv")

            with out_path.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle))

        self.assertEqual(len(rows), 12)
        self.assertEqual(rows[0]["scaffold"], "miniRhGC")
        self.assertEqual(rows[0]["seq_index"], "1")
        self.assertEqual(rows[0]["aa"], "M")
        self.assertEqual(rows[0]["alignment_col"], "1")
        self.assertEqual(rows[0]["review_status"], "draft")

    def test_apply_position_map_sets_protected_and_mutable_positions(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            scaffolds_path = self.write_scaffold_json(tmp_path)
            scaffolds = load_scaffolds(scaffolds_path)
            map_path = write_draft_position_map(scaffolds, tmp_path / "position_map.csv")

            with map_path.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle))
            for row in rows:
                if row["scaffold"] == "miniRhGC" and row["seq_index"] == "5":
                    row["role"] = "retinal_binding_lys"
                    row["protected"] = "true"
                    row["review_status"] = "reviewed"
                if row["scaffold"] == "miniRhGC" and row["seq_index"] == "4":
                    row["role"] = "retinal_pocket spectral_tuning"
                    row["mutable"] = "true"
                    row["allowed_mutations"] = "F;Y"
                    row["review_status"] = "reviewed"
            with map_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
                writer.writeheader()
                writer.writerows(rows)

            output_path = apply_position_map_to_scaffolds(
                scaffolds_path,
                map_path,
                tmp_path / "reviewed_scaffolds.json",
            )
            reviewed = load_scaffolds(output_path)
            candidates = generate_single_mutants(reviewed)

        rhgc = reviewed[0]
        self.assertEqual(rhgc.protected_positions, {5})
        self.assertEqual(rhgc.mutable_positions[0].position, 4)
        self.assertEqual(rhgc.mutable_positions[0].allowed, ["F", "Y"])
        self.assertEqual(
            [candidate.candidate_id for candidate in candidates],
            ["miniRhGC_p4AtoF", "miniRhGC_p4AtoY"],
        )


if __name__ == "__main__":
    unittest.main()
