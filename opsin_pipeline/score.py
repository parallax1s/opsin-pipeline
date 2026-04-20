from __future__ import annotations

from dataclasses import replace

from .schemas import Candidate, Scaffold


def rank_candidates(
    candidates: list[Candidate],
    scaffolds: list[Scaffold],
    target_family: str | None = None,
    target_phenotype: str | None = None,
) -> list[Candidate]:
    scaffold_by_name = {scaffold.name: scaffold for scaffold in scaffolds}
    scored = [
        score_candidate(
            candidate,
            scaffold_by_name[candidate.scaffold_name],
            target_family=target_family,
            target_phenotype=target_phenotype,
        )
        for candidate in candidates
    ]
    return sorted(
        scored,
        key=lambda candidate: (
            candidate.scores.get("total", 0.0),
            candidate.scaffold_name,
            candidate.candidate_id,
        ),
        reverse=True,
    )


def score_candidate(
    candidate: Candidate,
    scaffold: Scaffold,
    target_family: str | None = None,
    target_phenotype: str | None = None,
) -> Candidate:
    scores = {
        "family_match": 0.0,
        "phenotype_match": 0.0,
        "retinal_pocket": 0.0,
        "assay_readiness": 0.0,
        "mutation_penalty": -0.2 * len(candidate.mutations),
        "protected_penalty": -5.0 if candidate.has_protected_violation else 0.0,
    }

    if target_family and scaffold.family.lower() == target_family.lower():
        scores["family_match"] = 3.0
    if target_phenotype and target_phenotype in scaffold.target_phenotypes:
        scores["phenotype_match"] = 2.0

    reason_text = candidate.reason_summary.lower()
    if "retinal_pocket" in reason_text:
        scores["retinal_pocket"] = 2.0
    if {"growth_selection", "biochemical"} & set(scaffold.assay_architectures):
        scores["assay_readiness"] = 1.0

    scores["total"] = round(sum(scores.values()), 3)
    tags = sorted(set(candidate.tags + _score_tags(scores)))
    return replace(candidate, scores=scores, tags=tags)


def _score_tags(scores: dict[str, float]) -> list[str]:
    tags = []
    if scores.get("family_match", 0) > 0:
        tags.append("target_family")
    if scores.get("phenotype_match", 0) > 0:
        tags.append("target_phenotype")
    if scores.get("assay_readiness", 0) > 0:
        tags.append("assay_ready")
    return tags

