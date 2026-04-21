"""ESM2 backend tests (spec §5 + §9).

The module must always be importable. Failure modes:

- torch / transformers missing -> PLMBackendUnavailableError at score() time
- unknown --model alias -> ModelUnknownError at construction time
- sequence > max_length -> SequenceTooLongError at score() time

A real end-to-end test that actually downloads weights is opt-in via
``OPSIN_ENABLE_ESM_TESTS=1`` so the default 117-test suite never touches the
network.
"""
import importlib
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from opsin_pipeline.plm.esm import (
    ESM2Scorer,
    ModelUnknownError,
    PLMBackendUnavailableError,
    SequenceTooLongError,
    resolve_model_id,
)
from opsin_pipeline.plm.scorer import MutationRequest


class ModelResolutionTests(unittest.TestCase):
    def test_known_alias_resolves_to_hf_id(self):
        self.assertEqual(
            resolve_model_id("esm2_t12_35M_UR50D"),
            "facebook/esm2_t12_35M_UR50D",
        )

    def test_raw_hf_repo_id_passes_through(self):
        self.assertEqual(
            resolve_model_id("my-org/custom-esm2"),
            "my-org/custom-esm2",
        )

    def test_unknown_alias_raises(self):
        with self.assertRaises(ModelUnknownError) as ctx:
            resolve_model_id("bogus_model")
        self.assertIn("bogus_model", str(ctx.exception))


class LazyImportTests(unittest.TestCase):
    def test_module_imports_without_torch_or_transformers(self):
        """The plm.esm module must be importable even if torch / transformers
        are missing; only score() should trip the lazy-import path."""
        mod = importlib.import_module("opsin_pipeline.plm.esm")
        self.assertTrue(hasattr(mod, "ESM2Scorer"))
        self.assertTrue(hasattr(mod, "PLMBackendUnavailableError"))

    def test_construction_without_torch_is_ok(self):
        """Creating an ESM2Scorer instance must not attempt to import torch."""
        scorer = ESM2Scorer(model_alias="esm2_t12_35M_UR50D")
        self.assertEqual(scorer.model_id, "esm2_t12_35M_UR50D")

    def test_score_without_transformers_raises_backend_error(self):
        """transformers isn't installed in the default env; score() should raise
        PLMBackendUnavailableError with an actionable hint."""
        scorer = ESM2Scorer(model_alias="esm2_t12_35M_UR50D")
        with self.assertRaises(PLMBackendUnavailableError) as ctx:
            scorer.score(
                scaffold_name="demo",
                sequence="MKTAYK",
                mutations=[MutationRequest(position=2, from_aa="K", to_aa="A")],
            )
        msg = str(ctx.exception)
        self.assertIn("torch", msg.lower())
        self.assertIn("--mock", msg)


class ConstructionValidationTests(unittest.TestCase):
    def test_invalid_alias_fails_at_construction_time(self):
        with self.assertRaises(ModelUnknownError):
            ESM2Scorer(model_alias="not_a_real_model")


class SequenceLengthTests(unittest.TestCase):
    def test_sequence_over_max_length_raises(self):
        """Checked before the model is loaded, so this test works without torch."""
        scorer = ESM2Scorer(model_alias="esm2_t12_35M_UR50D", max_length=10)
        with self.assertRaises(SequenceTooLongError):
            scorer.score(
                scaffold_name="demo",
                sequence="A" * 11,
                mutations=[],
            )


class CLIESM2RoutingTests(unittest.TestCase):
    """With no --mock, `cli plm` must invoke ESM2Scorer and fail cleanly when
    transformers isn't installed."""

    def test_cli_non_mock_reports_backend_missing(self):
        env = {**os.environ, "PYTHONPATH": str(ROOT)}
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            scaffolds_path = tmp_path / "s.json"
            scaffolds_path.write_text(
                '{"scaffolds":[{"name":"demo","family":"S","sequence":"AFYK",'
                '"mutable_positions":[{"position":1,"allowed":["F"]}]}]}',
                encoding="utf-8",
            )
            result = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "opsin_pipeline.cli",
                    "plm",
                    "--scaffold", "demo",
                    "--scaffolds", str(scaffolds_path),
                    "--out", str(tmp_path / "p.json"),
                    # no --mock; ESM2Scorer should load, then fail on transformers import
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
                env=env,
            )

        self.assertNotEqual(result.returncode, 0)
        combined = (result.stderr + result.stdout).lower()
        # Either the backend-missing message or a torch/transformers mention
        self.assertTrue(
            "torch" in combined or "transformers" in combined or "--mock" in combined,
            msg=f"CLI error wasn't actionable: {combined}",
        )

    def test_cli_non_mock_unknown_model_exits_with_actionable_error(self):
        env = {**os.environ, "PYTHONPATH": str(ROOT)}
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            scaffolds_path = tmp_path / "s.json"
            scaffolds_path.write_text(
                '{"scaffolds":[{"name":"demo","family":"S","sequence":"AFYK",'
                '"mutable_positions":[{"position":1,"allowed":["F"]}]}]}',
                encoding="utf-8",
            )
            result = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "opsin_pipeline.cli",
                    "plm",
                    "--scaffold", "demo",
                    "--scaffolds", str(scaffolds_path),
                    "--out", str(tmp_path / "p.json"),
                    "--model", "not_a_real_model",
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
                env=env,
            )

        self.assertNotEqual(result.returncode, 0)
        self.assertIn("not_a_real_model", result.stderr + result.stdout)


@unittest.skipUnless(
    os.environ.get("OPSIN_ENABLE_ESM_TESTS") == "1",
    "Opt-in: set OPSIN_ENABLE_ESM_TESTS=1 to run the real ESM2 smoke test "
    "(downloads ~35M-param model from HuggingFace on first run).",
)
class RealESM2SmokeTest(unittest.TestCase):
    """End-to-end sanity check against a tiny ESM2 model. Runs only when the
    env var is set; default CI must never download weights."""

    def test_real_scorer_returns_finite_deltas(self):
        scorer = ESM2Scorer(model_alias="esm2_t6_8M_UR50D")  # smallest (~8M)
        pset = scorer.score(
            scaffold_name="demo",
            sequence="MKTAYIAKQRQISFVKSHFSR",  # 21 aa, well under any context window
            mutations=[
                MutationRequest(position=2, from_aa="K", to_aa="A"),
                MutationRequest(position=8, from_aa="K", to_aa="R"),
            ],
        )
        self.assertEqual(len(pset.predictions), 2)
        for p in pset.predictions:
            self.assertTrue(
                p.log_likelihood_delta == p.log_likelihood_delta,  # not NaN
                msg=f"delta is NaN for {p}",
            )
            # Loose bound — real ESM2 deltas typically within [-20, +10].
            self.assertGreater(p.log_likelihood_delta, -50.0)
            self.assertLess(p.log_likelihood_delta, 50.0)


if __name__ == "__main__":
    unittest.main()
