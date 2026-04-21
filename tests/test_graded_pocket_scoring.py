"""Graded pocket scoring (spec §7.4–§7.5).

Default is OFF: with ``use_graded_pocket=False`` the existing reason-string rule
produces identical numbers to main. With ``use_graded_pocket=True``, mutations
with ``distance_to_retinal`` populated switch to graded bands (strong <=4 A,
medium <=5.5 A, none) using ``max`` aggregation across a multi-mutant candidate.
Scaffolds without pocket data silently fall back to the reason-string rule.
"""
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from opsin_pipeline.generate import generate_candidates, generate_single_mutants
from opsin_pipeline.ingest import load_scaffolds
from opsin_pipeline.schemas import Candidate, Mutation, MutablePosition, Scaffold
from opsin_pipeline.score import (
    DEFAULT_MEDIUM_POINTS,
    DEFAULT_STRONG_POINTS,
    rank_candidates,
    score_candidate,
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


def _scaffold_with_pocket_map(tmp_path: Path) -> tuple[Path, Path]:
    atoms = parse_pdb(SYNTHETIC)
    pocket = compute_pocket(
        atoms,
        identify_ligands(atoms),
        scaffold_name="demo",
        pdb_path=str(SYNTHETIC),
    )
    mapped = apply_offset_mapping(pocket, offset=0)
    pocket_path = write_pocket_map(mapped, tmp_path / "pocket.json")
    scaffold_payload = {
        "scaffolds": [
            {
                "name": "demo",
                "family": "Synthetic",
                "sequence": "AFYKWVSG",
                "target_phenotypes": ["spectral_tuning"],
                "assay_architectures": ["biochemical"],
                "protected_positions": [],
                "pocket_map_path": pocket_path.name,
                "mutable_positions": [
                    # position 1 -> distance 2.0 A (strong)
                    {"position": 1, "allowed": ["F"], "reason": "pocket_test"},
                    # position 4 -> distance 4.5 A (medium)
                    {"position": 4, "allowed": ["F"], "reason": "pocket_test"},
                    # position 8 -> outside cutoff -> no distance
                    {"position": 8, "allowed": ["F"], "reason": "far_from_pocket"},
                ],
            }
        ]
    }
    scaffolds_path = tmp_path / "scaffolds.json"
    scaffolds_path.write_text(json.dumps(scaffold_payload), encoding="utf-8")
    return scaffolds_path, pocket_path


class FlagOffPreservesLegacyBehaviorTests(unittest.TestCase):
    def test_no_pocket_data_no_flag_gives_reason_string_score(self):
        scaffold = Scaffold(
            name="legacy",
            family="RhGC",
            sequence="MKTAYK",
            mutable_positions=[
                MutablePosition(position=4, allowed=["F"], reason="retinal_pocket test"),
            ],
            target_phenotypes=["spectral_tuning"],
        )
        candidate = Candidate(
            candidate_id="legacy_p4AtoF",
            scaffold_name="legacy",
            family="RhGC",
            mutations=[Mutation(position=4, from_aa="A", to_aa="F", reason="retinal_pocket test")],
            target_phenotypes=["spectral_tuning"],
            assay_architectures=[],
        )

        scored = score_candidate(candidate, scaffold)

        self.assertEqual(scored.scores.get("retinal_pocket"), 2.0)
        self.assertNotIn("retinal_pocket_strong", scored.scores)
        self.assertNotIn("retinal_pocket_medium", scored.scores)

    def test_pocket_data_present_but_flag_off_still_uses_reason_string(self):
        with tempfile.TemporaryDirectory() as tmp:
            scaffolds_path, _ = _scaffold_with_pocket_map(Path(tmp))
            scaffolds = load_scaffolds(scaffolds_path)
        candidates = generate_single_mutants(scaffolds)

        ranked = rank_candidates(candidates, scaffolds)  # use_graded_pocket defaults False

        for candidate in ranked:
            self.assertIn("retinal_pocket", candidate.scores)
            self.assertNotIn("retinal_pocket_strong", candidate.scores)
            self.assertNotIn("retinal_pocket_medium", candidate.scores)


class GradedPocketTests(unittest.TestCase):
    def test_strong_band_pays_strong_points(self):
        scaffold = Scaffold(
            name="demo",
            family="Synthetic",
            sequence="AFYK",
            mutable_positions=[
                MutablePosition(position=1, allowed=["F"], reason="pocket", distance_to_retinal=2.0),
            ],
        )
        candidate = Candidate(
            candidate_id="demo_p1AtoF",
            scaffold_name="demo",
            family="Synthetic",
            mutations=[
                Mutation(position=1, from_aa="A", to_aa="F", reason="pocket", distance_to_retinal=2.0),
            ],
            target_phenotypes=[],
            assay_architectures=[],
        )

        scored = score_candidate(candidate, scaffold, use_graded_pocket=True)

        self.assertEqual(scored.scores["retinal_pocket_strong"], DEFAULT_STRONG_POINTS)
        self.assertEqual(scored.scores["retinal_pocket_medium"], 0.0)
        self.assertNotIn("retinal_pocket", scored.scores)
        self.assertIn("retinal_pocket_strong", scored.tags)

    def test_medium_band_pays_medium_points(self):
        scaffold = Scaffold(
            name="demo",
            family="Synthetic",
            sequence="AFYK",
            mutable_positions=[
                MutablePosition(position=1, allowed=["F"], reason="pocket", distance_to_retinal=4.8),
            ],
        )
        candidate = Candidate(
            candidate_id="demo_p1AtoF",
            scaffold_name="demo",
            family="Synthetic",
            mutations=[
                Mutation(position=1, from_aa="A", to_aa="F", reason="pocket", distance_to_retinal=4.8),
            ],
            target_phenotypes=[],
            assay_architectures=[],
        )

        scored = score_candidate(candidate, scaffold, use_graded_pocket=True)

        self.assertEqual(scored.scores["retinal_pocket_medium"], DEFAULT_MEDIUM_POINTS)
        self.assertEqual(scored.scores["retinal_pocket_strong"], 0.0)

    def test_none_band_pays_zero(self):
        scaffold = Scaffold(
            name="demo",
            family="Synthetic",
            sequence="AFYK",
            mutable_positions=[
                MutablePosition(position=1, allowed=["F"], reason="pocket", distance_to_retinal=8.5),
            ],
        )
        candidate = Candidate(
            candidate_id="demo_p1AtoF",
            scaffold_name="demo",
            family="Synthetic",
            mutations=[
                Mutation(position=1, from_aa="A", to_aa="F", reason="pocket", distance_to_retinal=8.5),
            ],
            target_phenotypes=[],
            assay_architectures=[],
        )

        scored = score_candidate(candidate, scaffold, use_graded_pocket=True)

        self.assertEqual(scored.scores["retinal_pocket_strong"], 0.0)
        self.assertEqual(scored.scores["retinal_pocket_medium"], 0.0)

    def test_max_aggregation_across_multi_mutant(self):
        """Per spec §7.2: multi-mutant pocket score is max across mutations, not sum."""
        scaffold = Scaffold(
            name="demo",
            family="Synthetic",
            sequence="AFYKW",
            mutable_positions=[
                MutablePosition(position=1, allowed=["F"], reason="pocket", distance_to_retinal=2.0),
                MutablePosition(position=4, allowed=["F"], reason="pocket", distance_to_retinal=4.8),
            ],
        )
        candidate = Candidate(
            candidate_id="demo_p1AtoF__p4KtoF",
            scaffold_name="demo",
            family="Synthetic",
            mutations=[
                Mutation(position=1, from_aa="A", to_aa="F", distance_to_retinal=2.0),
                Mutation(position=4, from_aa="K", to_aa="F", distance_to_retinal=4.8),
            ],
            target_phenotypes=[],
            assay_architectures=[],
        )

        scored = score_candidate(candidate, scaffold, use_graded_pocket=True)

        # max(strong, medium) -> strong; medium key must be 0 even though one mutation is medium
        self.assertEqual(scored.scores["retinal_pocket_strong"], DEFAULT_STRONG_POINTS)
        self.assertEqual(scored.scores["retinal_pocket_medium"], 0.0)

    def test_graded_fallback_when_candidate_has_no_distance_data(self):
        """Spec §7.3 backward-compat: scaffolds without pocket data fall back to the
        reason-string rule even when the flag is on, so untouched scaffolds don't
        silently lose their score."""
        scaffold = Scaffold(
            name="legacy_in_mixed_run",
            family="RhGC",
            sequence="MKTAYK",
            mutable_positions=[
                MutablePosition(position=4, allowed=["F"], reason="retinal_pocket"),
            ],
        )
        candidate = Candidate(
            candidate_id="legacy_p4AtoF",
            scaffold_name="legacy_in_mixed_run",
            family="RhGC",
            mutations=[Mutation(position=4, from_aa="A", to_aa="F", reason="retinal_pocket")],
            target_phenotypes=[],
            assay_architectures=[],
        )

        scored = score_candidate(candidate, scaffold, use_graded_pocket=True)

        # falls back to legacy reason-string
        self.assertEqual(scored.scores.get("retinal_pocket"), 2.0)
        self.assertNotIn("retinal_pocket_strong", scored.scores)
        self.assertIn("pocket_graded_fallback", scored.tags)

    def test_custom_band_points_honored(self):
        scaffold = Scaffold(
            name="demo",
            family="Synthetic",
            sequence="AFYK",
            mutable_positions=[
                MutablePosition(position=1, allowed=["F"], reason="pocket", distance_to_retinal=2.0),
            ],
        )
        candidate = Candidate(
            candidate_id="demo_p1AtoF",
            scaffold_name="demo",
            family="Synthetic",
            mutations=[Mutation(position=1, from_aa="A", to_aa="F", distance_to_retinal=2.0)],
            target_phenotypes=[],
            assay_architectures=[],
        )

        scored = score_candidate(
            candidate,
            scaffold,
            use_graded_pocket=True,
            strong_points=1.5,
        )

        self.assertEqual(scored.scores["retinal_pocket_strong"], 1.5)


class EndToEndGradedRunTests(unittest.TestCase):
    """Integration: pocket data on a scaffold, --graded-pocket on, full rank pipeline."""

    def test_graded_run_reorders_candidates_by_band(self):
        with tempfile.TemporaryDirectory() as tmp:
            scaffolds_path, _ = _scaffold_with_pocket_map(Path(tmp))
            scaffolds = load_scaffolds(scaffolds_path)
        candidates = generate_single_mutants(scaffolds)

        graded = rank_candidates(candidates, scaffolds, use_graded_pocket=True)
        legacy = rank_candidates(candidates, scaffolds, use_graded_pocket=False)

        # graded run: position 1 (strong) > position 4 (medium) > position 8 (none, falls to reason)
        graded_by_pos = {c.mutations[0].position: c for c in graded}
        self.assertEqual(
            graded_by_pos[1].scores.get("retinal_pocket_strong"), DEFAULT_STRONG_POINTS
        )
        self.assertEqual(
            graded_by_pos[4].scores.get("retinal_pocket_medium"), DEFAULT_MEDIUM_POINTS
        )
        # position 8 in graded mode: no distance -> graded_fallback_reason; reason="far_from_pocket"
        # doesn't contain retinal_pocket -> fallback score is 0
        self.assertEqual(graded_by_pos[8].scores.get("retinal_pocket"), 0.0)

        # sanity: legacy run uses a single 'retinal_pocket' key for everyone; reason strings
        # here don't contain the tag so all score 0 on the pocket dimension
        for candidate in legacy:
            self.assertEqual(candidate.scores.get("retinal_pocket"), 0.0)

    def test_graded_run_produces_different_total_scores_than_legacy(self):
        with tempfile.TemporaryDirectory() as tmp:
            scaffolds_path, _ = _scaffold_with_pocket_map(Path(tmp))
            scaffolds = load_scaffolds(scaffolds_path)
        candidates = generate_single_mutants(scaffolds)

        graded_totals = {
            c.candidate_id: c.scores["total"]
            for c in rank_candidates(candidates, scaffolds, use_graded_pocket=True)
        }
        legacy_totals = {
            c.candidate_id: c.scores["total"]
            for c in rank_candidates(candidates, scaffolds, use_graded_pocket=False)
        }

        # at least one candidate must score differently under graded vs legacy, otherwise
        # the flag did nothing
        self.assertTrue(
            any(graded_totals[k] != legacy_totals[k] for k in graded_totals),
            "graded vs legacy totals are identical — flag had no effect",
        )


class ReportSignalLabelTests(unittest.TestCase):
    def test_pocket_signal_label_in_decision_report(self):
        from opsin_pipeline.report import write_decision_report

        with tempfile.TemporaryDirectory() as tmp:
            scaffolds_path, _ = _scaffold_with_pocket_map(Path(tmp))
            scaffolds = load_scaffolds(scaffolds_path)
            candidates = generate_single_mutants(scaffolds)
            ranked = rank_candidates(candidates, scaffolds, use_graded_pocket=True)

            report_path = write_decision_report(
                ranked,
                scaffolds,
                Path(tmp) / "report.md",
                use_graded_pocket=True,
            )
            text = report_path.read_text(encoding="utf-8")

        self.assertIn("Pocket signal:", text)
        self.assertIn("graded", text)

    def test_pocket_signal_label_reports_legacy_when_flag_off(self):
        from opsin_pipeline.report import write_decision_report

        with tempfile.TemporaryDirectory() as tmp:
            scaffolds_path, _ = _scaffold_with_pocket_map(Path(tmp))
            scaffolds = load_scaffolds(scaffolds_path)
            candidates = generate_single_mutants(scaffolds)
            ranked = rank_candidates(candidates, scaffolds, use_graded_pocket=False)

            report_path = write_decision_report(
                ranked, scaffolds, Path(tmp) / "report.md"
            )
            text = report_path.read_text(encoding="utf-8")

        self.assertIn("legacy", text)


if __name__ == "__main__":
    unittest.main()
