import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from opsin_pipeline.calibration import evaluate_calibration
from opsin_pipeline.generate import generate_candidates
from opsin_pipeline.ingest import load_scaffolds
from opsin_pipeline.score import rank_candidates


class GenerationAndRankingTests(unittest.TestCase):
    def write_scaffold_json(self, directory: Path) -> Path:
        payload = {
            "scaffolds": [
                {
                    "name": "RhGC_A",
                    "family": "RhGC",
                    "sequence": "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAN",
                    "target_phenotypes": ["spectral_tuning", "cGMP_activity"],
                    "assay_architectures": ["biochemical", "growth_selection"],
                    "protected_positions": [5],
                    "mutable_positions": [
                        {
                            "position": 8,
                            "allowed": ["F", "Y"],
                            "reason": "retinal_pocket spectral_tuning",
                        },
                        {
                            "position": 10,
                            "allowed": ["W"],
                            "reason": "retinal_pocket spectral_tuning",
                        },
                    ],
                },
                {
                    "name": "RhGC_B",
                    "family": "RhGC",
                    "sequence": "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAN",
                    "target_phenotypes": ["spectral_tuning", "cGMP_activity"],
                    "assay_architectures": ["biochemical"],
                    "protected_positions": [],
                    "mutable_positions": [
                        {
                            "position": 9,
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

    def test_generate_candidates_includes_combinations_up_to_max_mutations(self):
        with tempfile.TemporaryDirectory() as tmp:
            scaffolds = load_scaffolds(self.write_scaffold_json(Path(tmp)))

        candidates = generate_candidates(scaffolds, max_mutations=2)
        ids = [candidate.candidate_id for candidate in candidates]

        self.assertIn("RhGC_A_p8KtoF", ids)
        self.assertIn("RhGC_A_p8KtoF__p10RtoW", ids)
        self.assertIn("RhGC_A_p8KtoY__p10RtoW", ids)
        self.assertNotIn("RhGC_A_p5YtoF", ids)

    def test_rank_candidates_can_cap_per_scaffold(self):
        with tempfile.TemporaryDirectory() as tmp:
            scaffolds = load_scaffolds(self.write_scaffold_json(Path(tmp)))
        candidates = generate_candidates(scaffolds, max_mutations=2)

        ranked = rank_candidates(
            candidates,
            scaffolds,
            target_family="RhGC",
            target_phenotype="spectral_tuning",
            per_scaffold_cap=1,
        )

        scaffolds_seen = [candidate.scaffold_name for candidate in ranked]
        self.assertEqual(scaffolds_seen.count("RhGC_A"), 1)
        self.assertEqual(scaffolds_seen.count("RhGC_B"), 1)

    def test_rank_candidates_can_cap_per_position(self):
        with tempfile.TemporaryDirectory() as tmp:
            scaffolds = load_scaffolds(self.write_scaffold_json(Path(tmp)))
        candidates = generate_candidates(scaffolds, max_mutations=1)

        ranked = rank_candidates(
            candidates,
            scaffolds,
            target_family="RhGC",
            target_phenotype="spectral_tuning",
            per_position_cap=1,
        )

        position_keys = [candidate.position_key for candidate in ranked]
        self.assertEqual(len(position_keys), len(set(position_keys)))

    def test_evaluate_calibration_reports_known_positive_recall(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            scaffolds = load_scaffolds(self.write_scaffold_json(tmp_path))
            candidates = rank_candidates(
                generate_candidates(scaffolds, max_mutations=2),
                scaffolds,
                target_family="RhGC",
                target_phenotype="spectral_tuning",
            )
            calibration_path = tmp_path / "known_mutants.csv"
            with calibration_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=["candidate_id", "label", "evidence"],
                )
                writer.writeheader()
                writer.writerow(
                    {
                        "candidate_id": "RhGC_A_p8KtoF",
                        "label": "positive",
                        "evidence": "toy known spectral tuning mutation",
                    }
                )

            result = evaluate_calibration(candidates, calibration_path, top_n=3)

        self.assertEqual(result.known_positive_count, 1)
        self.assertEqual(result.top_n, 3)
        self.assertGreaterEqual(result.positive_recall_at_n, 0.0)
        self.assertLessEqual(result.positive_recall_at_n, 1.0)


if __name__ == "__main__":
    unittest.main()
