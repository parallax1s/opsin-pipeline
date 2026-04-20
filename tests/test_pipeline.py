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
from opsin_pipeline.report import write_candidate_csv, write_decision_report
from opsin_pipeline.score import rank_candidates


class OpsinPipelineTests(unittest.TestCase):
    def write_scaffold_json(self, directory: Path) -> Path:
        payload = {
            "scaffolds": [
                {
                    "name": "CaRhGC",
                    "family": "RhGC",
                    "sequence": "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAN",
                    "starting_lambda_nm": 527,
                    "target_phenotypes": ["spectral_tuning", "cGMP_activity"],
                    "assay_architectures": ["biochemical", "growth_selection"],
                    "protected_positions": [5, 12],
                    "mutable_positions": [
                        {
                            "position": 8,
                            "allowed": ["F", "Y"],
                            "reason": "retinal_pocket spectral_tuning",
                        },
                        {
                            "position": 12,
                            "allowed": ["A"],
                            "reason": "protected_control",
                        },
                    ],
                },
                {
                    "name": "ChR2",
                    "family": "ChR",
                    "sequence": "MDYGGALSAVGRELLFVTNPVVVNGSVLVPEDQCYCAGWIESR",
                    "starting_lambda_nm": 470,
                    "target_phenotypes": ["spectral_tuning", "photocurrent"],
                    "assay_architectures": ["growth_selection", "patch_clamp"],
                    "protected_positions": [10],
                    "mutable_positions": [
                        {
                            "position": 15,
                            "allowed": ["F"],
                            "reason": "retinal_pocket spectral_tuning",
                        }
                    ],
                },
            ]
        }
        path = directory / "scaffolds.json"
        path.write_text(json.dumps(payload), encoding="utf-8")
        return path

    def test_load_scaffolds_from_json(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self.write_scaffold_json(Path(tmp))

            scaffolds = load_scaffolds(path)

        self.assertEqual([s.name for s in scaffolds], ["CaRhGC", "ChR2"])
        self.assertEqual(scaffolds[0].family, "RhGC")
        self.assertEqual(scaffolds[0].mutable_positions[0].allowed, ["F", "Y"])
        self.assertEqual(scaffolds[0].starting_lambda_nm, 527.0)

    def test_generate_single_mutants_skips_protected_positions(self):
        with tempfile.TemporaryDirectory() as tmp:
            scaffolds = load_scaffolds(self.write_scaffold_json(Path(tmp)))

        candidates = generate_single_mutants(scaffolds)

        candidate_ids = [candidate.candidate_id for candidate in candidates]
        self.assertIn("CaRhGC_p8KtoF", candidate_ids)
        self.assertIn("CaRhGC_p8KtoY", candidate_ids)
        self.assertIn("ChR2_p15LtoF", candidate_ids)
        self.assertNotIn("CaRhGC_p12ItoA", candidate_ids)
        self.assertTrue(all(not c.has_protected_violation for c in candidates))

    def test_rank_candidates_prefers_target_family_and_phenotype(self):
        with tempfile.TemporaryDirectory() as tmp:
            scaffolds = load_scaffolds(self.write_scaffold_json(Path(tmp)))
        candidates = generate_single_mutants(scaffolds)

        ranked = rank_candidates(
            candidates,
            scaffolds,
            target_family="RhGC",
            target_phenotype="spectral_tuning",
        )

        self.assertEqual(ranked[0].scaffold_name, "CaRhGC")
        self.assertGreater(ranked[0].scores["total"], ranked[-1].scores["total"])
        self.assertIn("retinal_pocket", ranked[0].tags)

    def test_report_writes_csv_and_markdown(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            scaffolds = load_scaffolds(self.write_scaffold_json(tmp_path))
            ranked = rank_candidates(
                generate_single_mutants(scaffolds),
                scaffolds,
                target_family="RhGC",
                target_phenotype="spectral_tuning",
            )

            csv_path = write_candidate_csv(ranked, scaffolds, tmp_path / "ranked_candidates.csv")
            report_path = write_decision_report(
                ranked,
                scaffolds,
                tmp_path / "decision_report.md",
                target_family="RhGC",
                target_phenotype="spectral_tuning",
            )

            with csv_path.open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            report = report_path.read_text(encoding="utf-8")

        self.assertEqual(rows[0]["scaffold_name"], "CaRhGC")
        self.assertIn("starting_lambda_nm", rows[0].keys())
        self.assertIn("hamming", rows[0].keys())
        self.assertIn("# Opsin Pipeline Decision Report", report)
        self.assertIn("CaRhGC_p8KtoF", report)


if __name__ == "__main__":
    unittest.main()
