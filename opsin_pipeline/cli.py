from __future__ import annotations

import argparse
from pathlib import Path

from .generate import generate_single_mutants
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
    run_parser.add_argument("--top", type=int, default=10)

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
    candidates = generate_single_mutants(scaffolds)
    ranked = rank_candidates(
        candidates,
        scaffolds,
        target_family=args.target_family,
        target_phenotype=args.target_phenotype,
    )

    out_dir = Path(args.out)
    csv_path = write_candidate_csv(ranked, out_dir / "ranked_candidates.csv")
    report_path = write_decision_report(
        ranked,
        scaffolds,
        out_dir / "decision_report.md",
        target_family=args.target_family,
        target_phenotype=args.target_phenotype,
        top_n=args.top,
    )
    print(f"Wrote {csv_path}")
    print(f"Wrote {report_path}")


if __name__ == "__main__":
    main()
