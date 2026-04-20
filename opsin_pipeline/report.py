from __future__ import annotations

import csv
from pathlib import Path

from .schemas import Candidate, Scaffold


def write_candidate_csv(candidates: list[Candidate], path: str | Path) -> Path:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "candidate_id",
        "scaffold_name",
        "family",
        "mutations",
        "total_score",
        "tags",
        "reason",
    ]
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for candidate in candidates:
            writer.writerow(
                {
                    "candidate_id": candidate.candidate_id,
                    "scaffold_name": candidate.scaffold_name,
                    "family": candidate.family,
                    "mutations": candidate.mutation_summary,
                    "total_score": candidate.scores.get("total", 0.0),
                    "tags": ";".join(candidate.tags),
                    "reason": candidate.reason_summary,
                }
            )
    return output_path


def write_decision_report(
    candidates: list[Candidate],
    scaffolds: list[Scaffold],
    path: str | Path,
    target_family: str | None = None,
    target_phenotype: str | None = None,
    top_n: int = 10,
) -> Path:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Opsin Pipeline Decision Report",
        "",
        f"- Target family: {target_family or 'any'}",
        f"- Target phenotype: {target_phenotype or 'any'}",
        f"- Scaffolds screened: {len(scaffolds)}",
        f"- Candidates generated: {len(candidates)}",
        "",
        "## Top Candidates",
        "",
        "| Rank | Candidate | Scaffold | Score | Tags |",
        "|---:|---|---|---:|---|",
    ]
    for index, candidate in enumerate(candidates[:top_n], start=1):
        lines.append(
            "| {rank} | {candidate} | {scaffold} | {score} | {tags} |".format(
                rank=index,
                candidate=candidate.candidate_id,
                scaffold=candidate.scaffold_name,
                score=candidate.scores.get("total", 0.0),
                tags=", ".join(candidate.tags),
            )
        )
    lines.extend(
        [
            "",
            "## Next Validation Layer",
            "",
            "- Run structure and retinal-pocket checks on the top candidates.",
            "- Reserve QM/MM or excited-state calculations for a small shortlist.",
            "- Treat this MVP score as triage, not a final biophysical prediction.",
            "",
        ]
    )
    output_path.write_text("\n".join(lines), encoding="utf-8")
    return output_path

