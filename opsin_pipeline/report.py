from __future__ import annotations

import csv
from pathlib import Path

from .calibration import CalibrationReport
from .diversify import diversify_ranked
from .generate import GenerationStats
from .schemas import Candidate, Scaffold


def write_candidate_csv(
    candidates: list[Candidate],
    scaffolds: list[Scaffold],
    path: str | Path,
) -> Path:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    scaffold_lambda = {s.name: s.starting_lambda_nm for s in scaffolds}
    fieldnames = [
        "candidate_id",
        "scaffold_name",
        "family",
        "mutations",
        "hamming",
        "total_score",
        "starting_lambda_nm",
        "min_distance_to_retinal_A",
        "min_plm_log_likelihood_delta",
        "has_protected_violation",
        "tags",
        "reason",
    ]
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for candidate in candidates:
            lam = scaffold_lambda.get(candidate.scaffold_name)
            min_dist = _min_distance(candidate)
            min_plm = _min_plm_delta(candidate)
            writer.writerow(
                {
                    "candidate_id": candidate.candidate_id,
                    "scaffold_name": candidate.scaffold_name,
                    "family": candidate.family,
                    "mutations": candidate.mutation_summary,
                    "hamming": len(candidate.mutations),
                    "total_score": candidate.scores.get("total", 0.0),
                    "starting_lambda_nm": "" if lam is None else lam,
                    "min_distance_to_retinal_A": "" if min_dist is None else f"{min_dist:.2f}",
                    "min_plm_log_likelihood_delta": "" if min_plm is None else f"{min_plm:.3f}",
                    "has_protected_violation": int(candidate.has_protected_violation),
                    "tags": ";".join(candidate.tags),
                    "reason": candidate.reason_summary,
                }
            )
    return output_path


def _min_distance(candidate: Candidate) -> float | None:
    distances = [
        m.distance_to_retinal for m in candidate.mutations if m.distance_to_retinal is not None
    ]
    return min(distances) if distances else None


def _min_plm_delta(candidate: Candidate) -> float | None:
    """Worst-case plausibility across the candidate's mutations (spec §12.3).

    Lower = less plausible. ``min`` is chosen so a candidate with any extremely
    implausible mutation is visibly flagged in the CSV; aligns with how
    ``min_distance_to_retinal_A`` already aggregates.
    """
    deltas = [
        m.plm_log_likelihood_delta
        for m in candidate.mutations
        if m.plm_log_likelihood_delta is not None
    ]
    return min(deltas) if deltas else None


def _pocket_signal_label(
    candidates: list[Candidate], *, use_graded_pocket: bool
) -> str:
    """One-line label summarizing which pocket-scoring mode produced these numbers.

    Makes it obvious at a glance whether a report is from the legacy reason-string
    scoring or the structure-grounded graded scoring (spec §7.5).
    """
    if not use_graded_pocket:
        return "legacy (reason-string binary, flag off)"
    any_distance = any(
        m.distance_to_retinal is not None for c in candidates for m in c.mutations
    )
    if not any_distance:
        return "graded requested but no pocket data — fell back to reason-string"
    fallback_count = sum(
        1
        for c in candidates
        if all(m.distance_to_retinal is None for m in c.mutations)
    )
    if fallback_count:
        return (
            f"graded (distance bands); {fallback_count} candidates used "
            "reason-string fallback for lack of pocket data"
        )
    return "graded (distance bands) on all candidates"


def write_decision_report(
    candidates: list[Candidate],
    scaffolds: list[Scaffold],
    path: str | Path,
    target_family: str | None = None,
    target_phenotype: str | None = None,
    top_n: int = 10,
    per_scaffold_cap: int | None = None,
    per_position_cap: int | None = None,
    generation_stats: GenerationStats | None = None,
    calibration_report: CalibrationReport | None = None,
    use_graded_pocket: bool = False,
) -> Path:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    diversified = diversify_ranked(
        candidates,
        per_scaffold_cap=per_scaffold_cap,
        per_position_cap=per_position_cap,
        top_n=top_n,
    )
    violation_count = sum(1 for c in candidates if c.has_protected_violation)
    pocket_signal = _pocket_signal_label(candidates, use_graded_pocket=use_graded_pocket)

    lines = [
        "# Opsin Pipeline Decision Report",
        "",
        f"- Target family: {target_family or 'any'}",
        f"- Target phenotype: {target_phenotype or 'any'}",
        f"- Scaffolds screened: {len(scaffolds)}",
        f"- Candidates generated: {len(candidates)}",
        f"- Protected-residue violations: {violation_count}",
        f"- Diversity caps: per_scaffold_cap={per_scaffold_cap}, per_position_cap={per_position_cap}",
        f"- Pocket signal: {pocket_signal}",
        "",
    ]

    if scaffolds:
        lines.extend(_render_scaffold_summary(scaffolds))

    if generation_stats is not None:
        lines.extend(_render_generation_stats(generation_stats))

    lines.extend(
        [
            f"## Top {len(diversified)} Candidates (diversified)",
            "",
            "| Rank | Candidate | Scaffold | Hamming | Score | λmax start (nm) | Tags |",
            "|---:|---|---|---:|---:|---:|---|",
        ]
    )
    scaffold_lambda = {s.name: s.starting_lambda_nm for s in scaffolds}
    for index, candidate in enumerate(diversified, start=1):
        lam = scaffold_lambda.get(candidate.scaffold_name)
        lines.append(
            "| {rank} | {candidate} | {scaffold} | {hamming} | {score} | {lam} | {tags} |".format(
                rank=index,
                candidate=candidate.candidate_id,
                scaffold=candidate.scaffold_name,
                hamming=len(candidate.mutations),
                score=candidate.scores.get("total", 0.0),
                lam="—" if lam is None else lam,
                tags=", ".join(candidate.tags),
            )
        )
    lines.append("")

    if calibration_report is not None:
        lines.extend(_render_calibration(calibration_report))

    lines.extend(
        [
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


def _render_scaffold_summary(scaffolds: list[Scaffold]) -> list[str]:
    lines = [
        "## Scaffolds",
        "",
        "| Scaffold | Family | λmax start (nm) | Protected | Mutable |",
        "|---|---|---:|---:|---:|",
    ]
    for scaffold in scaffolds:
        lam = scaffold.starting_lambda_nm
        lines.append(
            "| {name} | {family} | {lam} | {prot} | {mut} |".format(
                name=scaffold.name,
                family=scaffold.family,
                lam="—" if lam is None else lam,
                prot=len(scaffold.protected_positions),
                mut=len(scaffold.mutable_positions),
            )
        )
    lines.append("")
    return lines


def _render_generation_stats(stats: GenerationStats) -> list[str]:
    lines = ["## Generation", "", f"- Total candidates: {stats.total_generated}"]
    if stats.per_scaffold_truncated:
        lines.append("- Scaffolds truncated by combinatorial cap:")
        for name, dropped in sorted(stats.per_scaffold_truncated.items()):
            kept = stats.per_scaffold_generated.get(name, 0)
            lines.append(f"  - {name}: kept {kept}, dropped {dropped}")
    lines.append("")
    return lines


def _render_calibration(report: CalibrationReport) -> list[str]:
    auroc = report.auroc_useful_vs_disruptive
    mrr = report.mean_reciprocal_rank_useful
    lines = [
        "## Calibration",
        "",
        f"- Useful matched: {report.useful_matched} / {report.useful_total}",
        f"- Disruptive matched: {report.disruptive_matched} / {report.disruptive_total}",
        f"- Neutral matched: {report.neutral_matched} / {report.neutral_total}",
        f"- AUROC useful vs disruptive: {auroc if auroc is not None else 'n/a (need both sides matched)'}",
        f"- MRR of useful mutants: {mrr if mrr is not None else 'n/a'}",
        f"- Useful in top {report.top_k}: {report.useful_in_top_k} (random baseline ≈ {report.random_baseline_useful_in_top_k})",
        f"- Disruptive in top {report.top_k}: {report.disruptive_in_top_k}",
    ]
    if report.unmatched_labels:
        preview = ", ".join(report.unmatched_labels[:10])
        ellipsis = "…" if len(report.unmatched_labels) > 10 else ""
        lines.append(
            f"- Unmatched calibration entries ({len(report.unmatched_labels)}): {preview}{ellipsis}"
        )
    lines.append("")
    return lines
