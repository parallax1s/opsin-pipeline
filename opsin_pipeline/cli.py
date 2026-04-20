from __future__ import annotations

import argparse
from pathlib import Path

from .generate import generate_single_mutants
from .ingest import load_scaffolds
from .report import write_candidate_csv, write_decision_report
from .score import rank_candidates


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the local opsin MVP pipeline")
    parser.add_argument("--scaffolds", required=True, help="Path to scaffold JSON")
    parser.add_argument("--out", required=True, help="Output directory")
    parser.add_argument("--target-family", default=None)
    parser.add_argument("--target-phenotype", default=None)
    parser.add_argument("--top", type=int, default=10)
    args = parser.parse_args()

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
