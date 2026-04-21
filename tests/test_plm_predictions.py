"""Tests for MockPLMScorer + PLMPredictionSet JSON round-trip (spec §10)."""
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from opsin_pipeline.plm.predictions import (
    PLMPrediction,
    PLMPredictionSet,
    deltas_by_position,
    predictions_from_dict,
    predictions_to_dict,
    read_predictions,
    sequence_sha256,
    write_predictions,
)
from opsin_pipeline.plm.scorer import (
    MockPLMScorer,
    MutationRequest,
    PositionOutOfRangeError,
    UnsupportedResidueError,
    WildTypeMismatchError,
)


BR_SEQUENCE_HEAD = "MLELLPTAVEGVSQAQITGRPEWIWLALGT"  # first 30 aa of BR_Hsal_BACR


class SequenceHashTests(unittest.TestCase):
    def test_sha256_is_stable_and_sensitive_to_changes(self):
        a = sequence_sha256("ACDE")
        b = sequence_sha256("ACDE")
        c = sequence_sha256("ACDF")  # last char differs
        self.assertEqual(a, b)
        self.assertNotEqual(a, c)


class MockScorerDeterminismTests(unittest.TestCase):
    def test_same_inputs_same_delta_across_runs(self):
        scorer = MockPLMScorer(seed=42)
        req = MutationRequest(position=5, from_aa="L", to_aa="A")
        pset1 = scorer.score(
            scaffold_name="BR_test", sequence=BR_SEQUENCE_HEAD, mutations=[req]
        )
        pset2 = scorer.score(
            scaffold_name="BR_test", sequence=BR_SEQUENCE_HEAD, mutations=[req]
        )
        self.assertEqual(
            pset1.predictions[0].log_likelihood_delta,
            pset2.predictions[0].log_likelihood_delta,
        )

    def test_different_seed_gives_different_delta(self):
        req = MutationRequest(position=5, from_aa="L", to_aa="A")
        d1 = MockPLMScorer(seed=0).score(
            scaffold_name="BR_test", sequence=BR_SEQUENCE_HEAD, mutations=[req]
        ).predictions[0].log_likelihood_delta
        d2 = MockPLMScorer(seed=1).score(
            scaffold_name="BR_test", sequence=BR_SEQUENCE_HEAD, mutations=[req]
        ).predictions[0].log_likelihood_delta
        self.assertNotEqual(d1, d2)

    def test_deltas_are_in_reasonable_range(self):
        # Spec §5: Mock deltas roughly in [-5.0, +2.0]
        scorer = MockPLMScorer()
        reqs = [
            MutationRequest(position=i, from_aa=BR_SEQUENCE_HEAD[i - 1], to_aa="A")
            for i in range(2, 11)
            if BR_SEQUENCE_HEAD[i - 1] != "A"
        ]
        pset = scorer.score(
            scaffold_name="BR_test", sequence=BR_SEQUENCE_HEAD, mutations=reqs
        )
        for p in pset.predictions:
            self.assertGreaterEqual(p.log_likelihood_delta, -5.0)
            self.assertLessEqual(p.log_likelihood_delta, 2.0)

    def test_silent_mutation_delta_is_zero(self):
        scorer = MockPLMScorer()
        # position 1 in BR head is M
        pset = scorer.score(
            scaffold_name="BR_test",
            sequence=BR_SEQUENCE_HEAD,
            mutations=[MutationRequest(position=1, from_aa="M", to_aa="M")],
        )
        self.assertEqual(pset.predictions[0].log_likelihood_delta, 0.0)


class MockScorerValidationTests(unittest.TestCase):
    def test_out_of_range_position_raises(self):
        with self.assertRaises(PositionOutOfRangeError):
            MockPLMScorer().score(
                scaffold_name="BR_test",
                sequence=BR_SEQUENCE_HEAD,
                mutations=[MutationRequest(position=999, from_aa="M", to_aa="A")],
            )

    def test_wild_type_mismatch_raises(self):
        # position 5 in BR head is L, not K
        with self.assertRaises(WildTypeMismatchError) as ctx:
            MockPLMScorer().score(
                scaffold_name="BR_test",
                sequence=BR_SEQUENCE_HEAD,
                mutations=[MutationRequest(position=5, from_aa="K", to_aa="A")],
            )
        self.assertIn("L", str(ctx.exception))

    def test_non_standard_target_aa_raises(self):
        with self.assertRaises(UnsupportedResidueError):
            MockPLMScorer().score(
                scaffold_name="BR_test",
                sequence=BR_SEQUENCE_HEAD,
                mutations=[MutationRequest(position=5, from_aa="L", to_aa="X")],
            )

    def test_non_standard_scaffold_residue_raises(self):
        # inject an X into the sequence at the requested position
        bad_sequence = "MABCDX"  # X at position 6
        with self.assertRaises(UnsupportedResidueError):
            MockPLMScorer().score(
                scaffold_name="bad_test",
                sequence=bad_sequence,
                mutations=[MutationRequest(position=6, from_aa="X", to_aa="A")],
            )


class JSONRoundTripTests(unittest.TestCase):
    def test_to_and_from_dict_preserves_content(self):
        pset = MockPLMScorer(seed=7).score(
            scaffold_name="BR_test",
            sequence=BR_SEQUENCE_HEAD,
            mutations=[
                MutationRequest(position=3, from_aa="E", to_aa="D"),
                MutationRequest(position=7, from_aa="T", to_aa="S"),
            ],
        )
        restored = predictions_from_dict(predictions_to_dict(pset))

        self.assertEqual(restored.scaffold_name, pset.scaffold_name)
        self.assertEqual(restored.sequence_sha256, pset.sequence_sha256)
        self.assertEqual(restored.model_id, pset.model_id)
        self.assertEqual(len(restored.predictions), len(pset.predictions))
        for a, b in zip(restored.predictions, pset.predictions):
            self.assertEqual(a.position, b.position)
            self.assertEqual(a.to_aa, b.to_aa)
            self.assertAlmostEqual(a.log_likelihood_delta, b.log_likelihood_delta, places=4)

    def test_write_and_read_file(self):
        pset = MockPLMScorer(seed=0).score(
            scaffold_name="BR_test",
            sequence=BR_SEQUENCE_HEAD,
            mutations=[MutationRequest(position=3, from_aa="E", to_aa="D")],
        )
        with tempfile.TemporaryDirectory() as tmp:
            path = write_predictions(pset, Path(tmp) / "plm.json")
            restored = read_predictions(path)

        self.assertEqual(restored.predictions[0].to_aa, "D")
        self.assertEqual(restored.model_id, "mock")


class DeltasByPositionTests(unittest.TestCase):
    def test_grouping_by_position_yields_to_aa_dict(self):
        pset = PLMPredictionSet(
            scaffold_name="BR_test",
            sequence_sha256="x",
            model_id="mock",
            predictions=[
                PLMPrediction(position=5, from_aa="L", to_aa="A", log_likelihood_delta=-1.0, model_id="mock"),
                PLMPrediction(position=5, from_aa="L", to_aa="V", log_likelihood_delta=0.5, model_id="mock"),
                PLMPrediction(position=8, from_aa="K", to_aa="R", log_likelihood_delta=-0.2, model_id="mock"),
            ],
        )

        grouped = deltas_by_position(pset)

        self.assertEqual(grouped[5], {"A": -1.0, "V": 0.5})
        self.assertEqual(grouped[8], {"R": -0.2})


if __name__ == "__main__":
    unittest.main()
