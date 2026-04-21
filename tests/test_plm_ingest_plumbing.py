"""Ingest + generate plumbing for PLM predictions (spec §8.2-§8.3)."""
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


def _write_plm(tmp_path: Path, scaffold_name: str, preds: list[tuple[int, str, str, float]]) -> Path:
    pset = PLMPredictionSet(
        scaffold_name=scaffold_name,
        sequence_sha256="fixture-hash",
        model_id="mock",
        predictions=[
            PLMPrediction(
                position=pos,
                from_aa=frm.upper(),
                to_aa=to.upper(),
                log_likelihood_delta=delta,
                model_id="mock",
            )
            for pos, frm, to, delta in preds
        ],
    )
    return write_predictions(pset, tmp_path / "plm.json")


def _write_scaffold(tmp_path: Path, payload: dict) -> Path:
    path = tmp_path / "scaffolds.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


class PLMIngestMergeTests(unittest.TestCase):
    def test_ingest_merges_plm_deltas_into_mutable_positions(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            plm_path = _write_plm(
                tmp_path,
                "demo",
                [
                    (1, "A", "F", -1.25),
                    (1, "A", "Y", 0.5),
                    (4, "K", "F", -0.7),
                ],
            )
            scaffold_payload = {
                "scaffolds": [
                    {
                        "name": "demo",
                        "family": "Synthetic",
                        "sequence": "AFYKWVSG",
                        "plm_predictions_path": plm_path.name,
                        "mutable_positions": [
                            {"position": 1, "allowed": ["F", "Y"]},
                            {"position": 4, "allowed": ["F"]},
                        ],
                    }
                ]
            }
            scaffolds_path = _write_scaffold(tmp_path, scaffold_payload)
            scaffolds = load_scaffolds(scaffolds_path)

        by_pos = {p.position: p for p in scaffolds[0].mutable_positions}
        self.assertEqual(by_pos[1].plm_log_likelihood_deltas, {"F": -1.25, "Y": 0.5})
        self.assertEqual(by_pos[4].plm_log_likelihood_deltas, {"F": -0.7})

    def test_unlisted_target_aas_are_dropped_from_merge(self):
        """If PLM JSON has a delta for a to_aa that isn't in the scaffold's `allowed`
        list, the merge quietly drops it so the generator never sees a stale entry."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            plm_path = _write_plm(
                tmp_path,
                "demo",
                [(1, "A", "F", -1.0), (1, "A", "W", 3.0)],  # W is NOT in allowed
            )
            scaffolds_path = _write_scaffold(
                tmp_path,
                {
                    "scaffolds": [
                        {
                            "name": "demo",
                            "family": "Synthetic",
                            "sequence": "AFYK",
                            "plm_predictions_path": plm_path.name,
                            "mutable_positions": [{"position": 1, "allowed": ["F"]}],
                        }
                    ]
                },
            )
            scaffolds = load_scaffolds(scaffolds_path)

        self.assertEqual(scaffolds[0].mutable_positions[0].plm_log_likelihood_deltas, {"F": -1.0})

    def test_positions_without_plm_predictions_stay_none(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            plm_path = _write_plm(tmp_path, "demo", [(1, "A", "F", -1.0)])
            scaffolds_path = _write_scaffold(
                tmp_path,
                {
                    "scaffolds": [
                        {
                            "name": "demo",
                            "family": "Synthetic",
                            "sequence": "AFYK",
                            "plm_predictions_path": plm_path.name,
                            "mutable_positions": [
                                {"position": 1, "allowed": ["F"]},
                                {"position": 4, "allowed": ["F"]},  # no PLM entry
                            ],
                        }
                    ]
                },
            )
            scaffolds = load_scaffolds(scaffolds_path)

        by_pos = {p.position: p for p in scaffolds[0].mutable_positions}
        self.assertEqual(by_pos[1].plm_log_likelihood_deltas, {"F": -1.0})
        self.assertIsNone(by_pos[4].plm_log_likelihood_deltas)


class PLMGeneratePropagationTests(unittest.TestCase):
    def test_generated_mutations_carry_plm_delta(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            plm_path = _write_plm(
                tmp_path,
                "demo",
                [(1, "A", "F", -1.25), (1, "A", "Y", 0.5), (4, "K", "F", -0.7)],
            )
            scaffolds_path = _write_scaffold(
                tmp_path,
                {
                    "scaffolds": [
                        {
                            "name": "demo",
                            "family": "Synthetic",
                            "sequence": "AFYKWVSG",
                            "plm_predictions_path": plm_path.name,
                            "mutable_positions": [
                                {"position": 1, "allowed": ["F", "Y"]},
                                {"position": 4, "allowed": ["F"]},
                            ],
                        }
                    ]
                },
            )
            scaffolds = load_scaffolds(scaffolds_path)

        candidates, _ = generate_candidates(scaffolds, max_mutations=1)
        by_id = {c.candidate_id: c for c in candidates}

        self.assertEqual(by_id["demo_p1AtoF"].mutations[0].plm_log_likelihood_delta, -1.25)
        self.assertEqual(by_id["demo_p1AtoY"].mutations[0].plm_log_likelihood_delta, 0.5)
        self.assertEqual(by_id["demo_p4KtoF"].mutations[0].plm_log_likelihood_delta, -0.7)

    def test_scaffold_without_plm_still_generates_candidates_cleanly(self):
        """Backward-compat: scaffolds without a plm_predictions_path produce
        candidates whose mutations have plm_log_likelihood_delta=None."""
        with tempfile.TemporaryDirectory() as tmp:
            scaffolds_path = _write_scaffold(
                Path(tmp),
                {
                    "scaffolds": [
                        {
                            "name": "legacy",
                            "family": "Synthetic",
                            "sequence": "AFYK",
                            "mutable_positions": [{"position": 1, "allowed": ["F"]}],
                        }
                    ]
                },
            )
            scaffolds = load_scaffolds(scaffolds_path)

        candidates, _ = generate_candidates(scaffolds, max_mutations=1)

        for c in candidates:
            for m in c.mutations:
                self.assertIsNone(m.plm_log_likelihood_delta)


if __name__ == "__main__":
    unittest.main()
