import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from opsin_pipeline.calibration import (
    CalibrationEntry,
    evaluate_ranking,
    load_calibration,
)
from opsin_pipeline.diversify import diversify_ranked
from opsin_pipeline.generate import generate_candidates, generate_single_mutants
from opsin_pipeline.ingest import load_scaffolds
from opsin_pipeline.report import write_candidate_csv, write_decision_report
from opsin_pipeline.schemas import Candidate, Mutation
from opsin_pipeline.score import rank_candidates, score_candidate


def _scaffold_payload() -> dict:
    return {
        "scaffolds": [
            {
                "name": "RhGC_A",
                "family": "RhGC",
                "sequence": "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAN",
                "starting_lambda_nm": 527,
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
                "starting_lambda_nm": 515,
                "target_phenotypes": ["spectral_tuning", "cGMP_activity"],
                "assay_architectures": ["biochemical"],
                "protected_positions": [],
                "mutable_positions": [
                    {
                        "position": 9,
                        "allowed": ["F"],
                        "reason": "helix_face expression",
                    }
                ],
            },
        ]
    }


def _write_scaffolds(directory: Path) -> Path:
    path = directory / "scaffolds.json"
    path.write_text(json.dumps(_scaffold_payload()), encoding="utf-8")
    return path


class GenerationTests(unittest.TestCase):
    def test_generate_candidates_returns_stats_tuple(self):
        with tempfile.TemporaryDirectory() as tmp:
            scaffolds = load_scaffolds(_write_scaffolds(Path(tmp)))

        candidates, stats = generate_candidates(scaffolds, max_mutations=1)

        self.assertEqual(stats.total_generated, len(candidates))
        self.assertGreater(stats.per_scaffold_generated["RhGC_A"], 0)

    def test_generate_candidates_includes_combinations_up_to_max_mutations(self):
        with tempfile.TemporaryDirectory() as tmp:
            scaffolds = load_scaffolds(_write_scaffolds(Path(tmp)))

        candidates, _ = generate_candidates(scaffolds, max_mutations=2)
        ids = {candidate.candidate_id for candidate in candidates}

        self.assertIn("RhGC_A_p8KtoF", ids)
        self.assertIn("RhGC_A_p8KtoF__p10RtoW", ids)
        self.assertIn("RhGC_A_p8KtoY__p10RtoW", ids)
        self.assertNotIn("RhGC_A_p5YtoF", ids)
        # positions in combined ids always ascending
        for candidate in candidates:
            positions = [m.position for m in candidate.mutations]
            self.assertEqual(positions, sorted(positions))

    def test_generation_cap_truncates_higher_order_first(self):
        with tempfile.TemporaryDirectory() as tmp:
            scaffolds = load_scaffolds(_write_scaffolds(Path(tmp)))

        # RhGC_A at Hamming 2 produces 3 singles + 2 doubles = 5. Cap to 3 -> keep singles only.
        candidates, stats = generate_candidates(
            scaffolds, max_mutations=2, max_combinations_per_scaffold=3
        )

        rhgc_a = [c for c in candidates if c.scaffold_name == "RhGC_A"]
        self.assertEqual(len(rhgc_a), 3)
        self.assertTrue(all(len(c.mutations) == 1 for c in rhgc_a))
        self.assertEqual(stats.per_scaffold_truncated["RhGC_A"], 2)

    def test_silent_mutations_are_skipped(self):
        payload = _scaffold_payload()
        # Position 8 is 'K'; adding K to allowed should produce no p8KtoK ghost
        payload["scaffolds"][0]["mutable_positions"][0]["allowed"] = ["K", "F"]

        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "s.json"
            path.write_text(json.dumps(payload), encoding="utf-8")
            scaffolds = load_scaffolds(path)

        candidates, _ = generate_candidates(scaffolds, max_mutations=2)

        self.assertFalse(any("p8KtoK" in c.candidate_id for c in candidates))
        self.assertTrue(any(c.candidate_id == "RhGC_A_p8KtoF" for c in candidates))

    def test_max_mutations_validated(self):
        with tempfile.TemporaryDirectory() as tmp:
            scaffolds = load_scaffolds(_write_scaffolds(Path(tmp)))
        with self.assertRaises(ValueError):
            generate_candidates(scaffolds, max_mutations=0)


class ScoringTests(unittest.TestCase):
    def test_rank_preserves_full_ranked_list(self):
        with tempfile.TemporaryDirectory() as tmp:
            scaffolds = load_scaffolds(_write_scaffolds(Path(tmp)))
        candidates, _ = generate_candidates(scaffolds, max_mutations=2)

        ranked = rank_candidates(
            candidates,
            scaffolds,
            target_family="RhGC",
            target_phenotype="spectral_tuning",
        )

        # rank_candidates no longer bakes diversity; full list survives
        self.assertEqual(len(ranked), len(candidates))
        self.assertGreaterEqual(ranked[0].scores["total"], ranked[-1].scores["total"])

    def test_protected_violation_flags_external_candidate(self):
        with tempfile.TemporaryDirectory() as tmp:
            scaffolds = load_scaffolds(_write_scaffolds(Path(tmp)))
        scaffold = scaffolds[0]
        external = Candidate(
            candidate_id="RhGC_A_p5YtoA",
            scaffold_name=scaffold.name,
            family=scaffold.family,
            mutations=[Mutation(position=5, from_aa="Y", to_aa="A")],
            target_phenotypes=scaffold.target_phenotypes,
            assay_architectures=scaffold.assay_architectures,
        )

        scored = score_candidate(external, scaffold)

        self.assertTrue(scored.has_protected_violation)
        self.assertEqual(scored.scores["protected_penalty"], -5.0)
        self.assertIn("protected_violation", scored.tags)


class DiversityTests(unittest.TestCase):
    def _ranked(self):
        with tempfile.TemporaryDirectory() as tmp:
            scaffolds = load_scaffolds(_write_scaffolds(Path(tmp)))
        candidates, _ = generate_candidates(scaffolds, max_mutations=2)
        return rank_candidates(
            candidates,
            scaffolds,
            target_family="RhGC",
            target_phenotype="spectral_tuning",
        )

    def test_per_scaffold_cap(self):
        diversified = diversify_ranked(self._ranked(), per_scaffold_cap=1, top_n=10)

        per_scaffold = {c.scaffold_name for c in diversified}
        self.assertEqual({"RhGC_A", "RhGC_B"}, per_scaffold)
        self.assertEqual(len(diversified), 2)

    def test_per_position_cap_counts_every_touched_position(self):
        diversified = diversify_ranked(self._ranked(), per_position_cap=1, top_n=20)

        touched: dict[tuple[str, int], int] = {}
        for candidate in diversified:
            for m in candidate.mutations:
                key = (candidate.scaffold_name, m.position)
                touched[key] = touched.get(key, 0) + 1
        self.assertTrue(all(count <= 1 for count in touched.values()))

    def test_top_n_respected(self):
        self.assertEqual(len(diversify_ranked(self._ranked(), top_n=2)), 2)


class CalibrationTests(unittest.TestCase):
    def _scaffolds(self):
        with tempfile.TemporaryDirectory() as tmp:
            return load_scaffolds(_write_scaffolds(Path(tmp)))

    def test_load_calibration_json_parses_all_categories(self):
        payload = {
            "calibration_sets": [
                {
                    "scaffold": "RhGC_A",
                    "source": "synthetic",
                    "known_useful": [
                        {
                            "label": "K8F",
                            "mutations": [{"position": 8, "to": "F"}],
                            "delta_lambda_nm": -12,
                        }
                    ],
                    "known_disruptive": [
                        {"label": "K8Y", "mutations": [{"position": 8, "to": "y"}]}
                    ],
                }
            ]
        }
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "cal.json"
            path.write_text(json.dumps(payload), encoding="utf-8")

            entries = load_calibration(path)

        categories = {e.category for e in entries}
        self.assertEqual(categories, {"useful", "disruptive"})
        useful = next(e for e in entries if e.category == "useful")
        disruptive = next(e for e in entries if e.category == "disruptive")
        self.assertEqual(useful.mutations, frozenset({(8, "F")}))
        self.assertEqual(useful.delta_lambda_nm, -12.0)
        # lowercase target amino acid normalized
        self.assertEqual(disruptive.mutations, frozenset({(8, "Y")}))

    def test_load_calibration_csv_parses_candidate_ids(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "cal.csv"
            with path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["candidate_id", "label", "evidence"])
                writer.writeheader()
                writer.writerows(
                    [
                        {"candidate_id": "RhGC_A_p8KtoF", "label": "positive", "evidence": ""},
                        {
                            "candidate_id": "RhGC_A_p8KtoF__p10RtoW",
                            "label": "positive",
                            "evidence": "toy double",
                        },
                        {"candidate_id": "RhGC_A_p10RtoW", "label": "neutral", "evidence": ""},
                        {"candidate_id": "RhGC_B_p9QtoF", "label": "negative", "evidence": ""},
                    ]
                )

            entries = load_calibration(path)

        by_category: dict[str, list] = {"useful": [], "neutral": [], "disruptive": []}
        for entry in entries:
            by_category[entry.category].append(entry)
        self.assertEqual(len(by_category["useful"]), 2)
        self.assertEqual(len(by_category["neutral"]), 1)
        self.assertEqual(len(by_category["disruptive"]), 1)
        # Multi-mutant id parsed back to the right set
        double = next(e for e in entries if len(e.mutations) == 2)
        self.assertEqual(double.mutations, frozenset({(8, "F"), (10, "W")}))
        self.assertEqual(double.scaffold_name, "RhGC_A")
        # Scaffold names containing underscores survive parsing
        self.assertTrue(all(e.scaffold_name in {"RhGC_A", "RhGC_B"} for e in entries))

    def test_load_calibration_csv_rejects_bad_label(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "cal.csv"
            path.write_text("candidate_id,label\nRhGC_A_p8KtoF,totally_made_up\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                load_calibration(path)

    def test_evaluate_ranking_separates_useful_from_disruptive(self):
        scaffolds = self._scaffolds()
        ranked = rank_candidates(
            generate_single_mutants(scaffolds),
            scaffolds,
            target_family="RhGC",
            target_phenotype="spectral_tuning",
        )

        # K8F is in a retinal_pocket position (bonus); p9F on RhGC_B scores lower (no bonus)
        calibration = [
            CalibrationEntry(
                scaffold_name="RhGC_A",
                label="K8F",
                mutations=frozenset({(8, "F")}),
                category="useful",
            ),
            CalibrationEntry(
                scaffold_name="RhGC_B",
                label="Q9F",
                mutations=frozenset({(9, "F")}),
                category="disruptive",
            ),
        ]

        report = evaluate_ranking(ranked, calibration, top_k=5)

        self.assertEqual(report.useful_matched, 1)
        self.assertEqual(report.disruptive_matched, 1)
        self.assertEqual(report.auroc_useful_vs_disruptive, 1.0)
        self.assertGreaterEqual(report.useful_in_top_k, 1)

    def test_evaluate_ranking_matches_multi_mutants(self):
        scaffolds = self._scaffolds()
        candidates, _ = generate_candidates(scaffolds, max_mutations=2)
        ranked = rank_candidates(candidates, scaffolds)

        calibration = [
            CalibrationEntry(
                scaffold_name="RhGC_A",
                label="K8F+R10W",
                mutations=frozenset({(8, "F"), (10, "W")}),
                category="useful",
            )
        ]

        report = evaluate_ranking(ranked, calibration)

        self.assertEqual(report.useful_matched, 1)

    def test_evaluate_ranking_reports_unmatched(self):
        scaffolds = self._scaffolds()
        ranked = rank_candidates(generate_single_mutants(scaffolds), scaffolds)
        calibration = [
            CalibrationEntry(
                scaffold_name="RhGC_A",
                label="nonexistent",
                mutations=frozenset({(99, "G")}),
                category="useful",
            )
        ]

        report = evaluate_ranking(ranked, calibration)

        self.assertEqual(report.matched, 0)
        self.assertEqual(report.unmatched, 1)
        self.assertEqual(report.unmatched_labels, ["nonexistent"])
        self.assertIsNone(report.auroc_useful_vs_disruptive)


class ReportIntegrationTests(unittest.TestCase):
    def test_decision_report_includes_scaffold_summary_generation_stats_and_calibration(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            scaffolds = load_scaffolds(_write_scaffolds(tmp_path))
            candidates, stats = generate_candidates(scaffolds, max_mutations=2)
            ranked = rank_candidates(
                candidates,
                scaffolds,
                target_family="RhGC",
                target_phenotype="spectral_tuning",
            )
            calibration = [
                CalibrationEntry(
                    scaffold_name="RhGC_A",
                    label="K8F",
                    mutations=frozenset({(8, "F")}),
                    category="useful",
                ),
                CalibrationEntry(
                    scaffold_name="RhGC_B",
                    label="Q9F",
                    mutations=frozenset({(9, "F")}),
                    category="disruptive",
                ),
            ]
            calibration_report = evaluate_ranking(ranked, calibration)

            csv_path = write_candidate_csv(ranked, scaffolds, tmp_path / "ranked.csv")
            report_path = write_decision_report(
                ranked,
                scaffolds,
                tmp_path / "report.md",
                target_family="RhGC",
                target_phenotype="spectral_tuning",
                top_n=5,
                per_scaffold_cap=2,
                generation_stats=stats,
                calibration_report=calibration_report,
            )
            report_text = report_path.read_text(encoding="utf-8")
            with csv_path.open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))

        self.assertIn("## Scaffolds", report_text)
        self.assertIn("## Generation", report_text)
        self.assertIn("## Calibration", report_text)
        self.assertIn("AUROC useful vs disruptive", report_text)
        # diversity cap: RhGC_A appears at most twice in the diversified top-N section
        top_section = report_text.split("## Top")[1]
        self.assertLessEqual(top_section.count("| RhGC_A |"), 2)
        # CSV has starting_lambda_nm populated from scaffold lookup
        rhgc_a_rows = [r for r in rows if r["scaffold_name"] == "RhGC_A"]
        self.assertTrue(all(r["starting_lambda_nm"] in {"527", "527.0"} for r in rhgc_a_rows))


if __name__ == "__main__":
    unittest.main()
