from __future__ import annotations

import argparse
from pathlib import Path

from .calibration import evaluate_ranking, load_calibration
from .generate import generate_candidates
from .ingest import load_scaffolds
from .position_map import apply_position_map_to_scaffolds, write_draft_position_map
from .report import write_candidate_csv, write_decision_report
from .score import rank_candidates


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the local opsin MVP pipeline")
    subparsers = parser.add_subparsers(dest="command")

    run_parser = subparsers.add_parser("run", help="Generate, score, and report candidates")
    run_parser.add_argument("--scaffolds", required=True, help="Path to scaffold JSON")
    run_parser.add_argument("--out", required=True, help="Output directory")
    run_parser.add_argument("--target-family", default=None)
    run_parser.add_argument("--target-phenotype", default=None)
    run_parser.add_argument("--top", type=int, default=10, help="Top-N size for the report")
    run_parser.add_argument(
        "--max-mutations",
        type=int,
        default=1,
        help="Maximum number of simultaneous mutations per candidate (1 = singles only).",
    )
    run_parser.add_argument(
        "--max-combinations-per-scaffold",
        type=int,
        default=None,
        help="Hard cap on candidates per scaffold; truncates higher-order mutants first.",
    )
    run_parser.add_argument(
        "--per-scaffold-cap",
        type=int,
        default=None,
        help="In the report's top-N, limit how many entries any one scaffold contributes.",
    )
    run_parser.add_argument(
        "--per-position-cap",
        type=int,
        default=None,
        help="In the report's top-N, limit how many entries touch any one (scaffold, position).",
    )
    run_parser.add_argument(
        "--calibration",
        default=None,
        help="Optional calibration file (.json or .csv) with known useful / disruptive mutations.",
    )
    run_parser.add_argument(
        "--calibration-top-k",
        type=int,
        default=20,
        help="Top-k used by the calibration evaluator.",
    )

    draft_parser = subparsers.add_parser("draft-position-map", help="Write a draft review CSV")
    draft_parser.add_argument("--scaffolds", required=True, help="Path to scaffold JSON")
    draft_parser.add_argument("--out", required=True, help="Output CSV path")

    apply_parser = subparsers.add_parser(
        "apply-position-map",
        help="Apply reviewed positions to a scaffold JSON file",
    )
    apply_parser.add_argument("--scaffolds", required=True, help="Source scaffold JSON")
    apply_parser.add_argument("--position-map", required=True, help="Reviewed position CSV")
    apply_parser.add_argument("--out", required=True, help="Output scaffold JSON")

    args = parser.parse_args()
    if args.command == "draft-position-map":
        scaffolds = load_scaffolds(args.scaffolds)
        path = write_draft_position_map(scaffolds, args.out)
        print(f"Wrote {path}")
        return

    if args.command == "apply-position-map":
        path = apply_position_map_to_scaffolds(args.scaffolds, args.position_map, args.out)
        print(f"Wrote {path}")
        return

    if args.command is None:
        parser.error("a command is required")

    scaffolds = load_scaffolds(args.scaffolds)
    candidates, gen_stats = generate_candidates(
        scaffolds,
        max_mutations=args.max_mutations,
        max_combinations_per_scaffold=args.max_combinations_per_scaffold,
    )
    ranked = rank_candidates(
        candidates,
        scaffolds,
        target_family=args.target_family,
        target_phenotype=args.target_phenotype,
    )

    calibration_report = None
    if args.calibration:
        entries = load_calibration(args.calibration)
        calibration_report = evaluate_ranking(
            ranked, entries, top_k=args.calibration_top_k
        )

    out_dir = Path(args.out)
    csv_path = write_candidate_csv(ranked, scaffolds, out_dir / "ranked_candidates.csv")
    report_path = write_decision_report(
        ranked,
        scaffolds,
        out_dir / "decision_report.md",
        target_family=args.target_family,
        target_phenotype=args.target_phenotype,
        top_n=args.top,
        per_scaffold_cap=args.per_scaffold_cap,
        per_position_cap=args.per_position_cap,
        generation_stats=gen_stats,
        calibration_report=calibration_report,
    )
    print(f"Wrote {csv_path}")
    print(f"Wrote {report_path}")
    if gen_stats.per_scaffold_truncated:
        for name, dropped in sorted(gen_stats.per_scaffold_truncated.items()):
            print(f"  truncated {dropped} candidates from {name} (cap)")
    if calibration_report is not None:
        auroc = calibration_report.auroc_useful_vs_disruptive
        print(
            "  calibration: useful {u}/{ut}, disruptive {d}/{dt}, AUROC {a}".format(
                u=calibration_report.useful_matched,
                ut=calibration_report.useful_total,
                d=calibration_report.disruptive_matched,
                dt=calibration_report.disruptive_total,
                a=auroc if auroc is not None else "n/a",
            )
        )


if __name__ == "__main__":
    main()
