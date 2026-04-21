"""Ensure the `min_plm_log_likelihood_delta` CSV column is populated when a
scaffold has a ``plm_predictions_path`` and blank otherwise (spec §8.4)."""
import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from opsin_pipeline.generate import generate_candidates
from opsin_pipeline.ingest import load_scaffolds
from opsin_pipeline.plm.predictions import PLMPrediction, PLMPredictionSet, write_predictions
from opsin_pipeline.report import write_candidate_csv
from opsin_pipeline.score import rank_candidates


def _fixture_with_plm(tmp_path: Path) -> Path:
    pset = PLMPredictionSet(
        scaffold_name="demo",
        sequence_sha256="fixture",
        model_id="mock",
        predictions=[
            PLMPrediction(position=1, from_aa="A", to_aa="F", log_likelihood_delta=-2.5, model_id="mock"),
            PLMPrediction(position=4, from_aa="K", to_aa="F", log_likelihood_delta=0.3, model_id="mock"),
        ],
    )
    plm_path = write_predictions(pset, tmp_path / "plm.json")
    payload = {
        "scaffolds": [
            {
                "name": "demo",
                "family": "Synthetic",
                "sequence": "AFYKWVSG",
                "plm_predictions_path": plm_path.name,
                "mutable_positions": [
                    {"position": 1, "allowed": ["F"]},
                    {"position": 4, "allowed": ["F"]},
                ],
            }
        ]
    }
    path = tmp_path / "scaffolds.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


class MinPLMDeltaColumnTests(unittest.TestCase):
    def test_single_mutant_emits_the_one_delta(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            scaffolds_path = _fixture_with_plm(tmp_path)
            scaffolds = load_scaffolds(scaffolds_path)
            candidates, _ = generate_candidates(scaffolds, max_mutations=1)
            ranked = rank_candidates(candidates, scaffolds)

            csv_path = write_candidate_csv(ranked, scaffolds, tmp_path / "ranked.csv")
            with csv_path.open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))

        by_id = {r["candidate_id"]: r for r in rows}
        self.assertEqual(by_id["demo_p1AtoF"]["min_plm_log_likelihood_delta"], "-2.500")
        self.assertEqual(by_id["demo_p4KtoF"]["min_plm_log_likelihood_delta"], "0.300")

    def test_multi_mutant_reports_min_across_mutations(self):
        """Per spec §12.3: aggregation is min (worst-case plausibility)."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            scaffolds_path = _fixture_with_plm(tmp_path)
            scaffolds = load_scaffolds(scaffolds_path)
            candidates, _ = generate_candidates(scaffolds, max_mutations=2)
            ranked = rank_candidates(candidates, scaffolds)

            csv_path = write_candidate_csv(ranked, scaffolds, tmp_path / "ranked.csv")
            with csv_path.open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))

        double = next(r for r in rows if r["candidate_id"] == "demo_p1AtoF__p4KtoF")
        # individual deltas are -2.5 and 0.3; min = -2.5
        self.assertEqual(double["min_plm_log_likelihood_delta"], "-2.500")

    def test_scaffold_without_plm_leaves_column_blank(self):
        """Backward-compat: no plm_predictions_path => blank cell, no error."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            payload = {
                "scaffolds": [
                    {
                        "name": "legacy",
                        "family": "Synthetic",
                        "sequence": "AFYK",
                        "mutable_positions": [{"position": 1, "allowed": ["F"]}],
                    }
                ]
            }
            scaffolds_path = tmp_path / "scaffolds.json"
            scaffolds_path.write_text(json.dumps(payload), encoding="utf-8")
            scaffolds = load_scaffolds(scaffolds_path)
            candidates, _ = generate_candidates(scaffolds, max_mutations=1)
            ranked = rank_candidates(candidates, scaffolds)

            csv_path = write_candidate_csv(ranked, scaffolds, tmp_path / "ranked.csv")
            with csv_path.open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))

        self.assertEqual(rows[0]["min_plm_log_likelihood_delta"], "")


if __name__ == "__main__":
    unittest.main()
