from __future__ import annotations

from dataclasses import replace

from .schemas import Candidate, Scaffold


# Spec §7.1 — module-level constants so they're trivially tunable.
# Default graded-pocket band points: strong <=4 A -> +2, medium 4-5.5 A -> +1, none -> 0.
DEFAULT_STRONG_POINTS = 2.0
DEFAULT_MEDIUM_POINTS = 1.0
DEFAULT_STRONG_MAX_A = 4.0
DEFAULT_MEDIUM_MAX_A = 5.5


def rank_candidates(
    candidates: list[Candidate],
    scaffolds: list[Scaffold],
    target_family: str | None = None,
    target_phenotype: str | None = None,
    use_graded_pocket: bool = False,
) -> list[Candidate]:
    scaffold_by_name = {scaffold.name: scaffold for scaffold in scaffolds}
    scored = [
        score_candidate(
            candidate,
            scaffold_by_name[candidate.scaffold_name],
            target_family=target_family,
            target_phenotype=target_phenotype,
            use_graded_pocket=use_graded_pocket,
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
    use_graded_pocket: bool = False,
    strong_max_A: float = DEFAULT_STRONG_MAX_A,
    medium_max_A: float = DEFAULT_MEDIUM_MAX_A,
    strong_points: float = DEFAULT_STRONG_POINTS,
    medium_points: float = DEFAULT_MEDIUM_POINTS,
) -> Candidate:
    """Score one candidate.

    When ``use_graded_pocket=False`` (default per spec §7.5), the retinal-pocket
    contribution is the legacy binary rule: +2.0 if any mutation's reason string
    contains ``retinal_pocket``.

    When ``use_graded_pocket=True``, the retinal_pocket key is replaced by two
    explicit band keys ``retinal_pocket_strong`` and ``retinal_pocket_medium``.
    Contribution = max over mutations of the band points for that mutation's
    ``distance_to_retinal`` (per spec §7.2: max aggregation, not sum). Mutations
    without distance data fall back to the reason-string rule only when the
    candidate has *no* distance data at all on any of its mutations.
    """
    violation = any(m.position in scaffold.protected_positions for m in candidate.mutations)

    scores: dict[str, float] = {
        "family_match": 0.0,
        "phenotype_match": 0.0,
        "assay_readiness": 0.0,
        "mutation_penalty": -0.2 * len(candidate.mutations),
        "protected_penalty": -5.0 if violation else 0.0,
    }

    if target_family and scaffold.family.lower() == target_family.lower():
        scores["family_match"] = 3.0
    if target_phenotype and target_phenotype in scaffold.target_phenotypes:
        scores["phenotype_match"] = 2.0

    pocket_mode = _apply_pocket_signal(
        scores,
        candidate,
        use_graded_pocket=use_graded_pocket,
        strong_max_A=strong_max_A,
        medium_max_A=medium_max_A,
        strong_points=strong_points,
        medium_points=medium_points,
    )

    if {"growth_selection", "biochemical"} & set(scaffold.assay_architectures):
        scores["assay_readiness"] = 1.0

    scores["total"] = round(sum(scores.values()), 3)
    tags = sorted(set(candidate.tags + _score_tags(scores, violation, pocket_mode)))
    return replace(candidate, scores=scores, tags=tags, has_protected_violation=violation)


def _apply_pocket_signal(
    scores: dict[str, float],
    candidate: Candidate,
    *,
    use_graded_pocket: bool,
    strong_max_A: float,
    medium_max_A: float,
    strong_points: float,
    medium_points: float,
) -> str:
    """Fill in pocket score keys. Returns the mode label used for tag reporting.

    Modes:
    - ``"legacy_reason"`` — the original reason-string rule (binary +2 if tag matches).
    - ``"graded_distance"`` — graded bands from distance_to_retinal on mutations.
    - ``"graded_fallback_reason"`` — user asked for graded but no mutation has distance
      data, so we fall back to the reason-string rule (still a single ``retinal_pocket``
      key so output stays readable).
    """
    distances = [
        m.distance_to_retinal for m in candidate.mutations if m.distance_to_retinal is not None
    ]

    if not use_graded_pocket:
        scores["retinal_pocket"] = 2.0 if _reason_has_pocket_tag(candidate) else 0.0
        return "legacy_reason"

    if not distances:
        # graded requested but no distance evidence on this candidate — stay in
        # backward-compat mode so untouched scaffolds produce identical scores
        scores["retinal_pocket"] = 2.0 if _reason_has_pocket_tag(candidate) else 0.0
        return "graded_fallback_reason"

    best = max(
        _band_points(
            d,
            strong_max_A=strong_max_A,
            medium_max_A=medium_max_A,
            strong_points=strong_points,
            medium_points=medium_points,
        )
        for d in distances
    )
    scores["retinal_pocket_strong"] = strong_points if best == strong_points else 0.0
    scores["retinal_pocket_medium"] = medium_points if best == medium_points else 0.0
    return "graded_distance"


def _reason_has_pocket_tag(candidate: Candidate) -> bool:
    return "retinal_pocket" in candidate.reason_summary.lower()


def _band_points(
    distance_A: float,
    *,
    strong_max_A: float,
    medium_max_A: float,
    strong_points: float,
    medium_points: float,
) -> float:
    if distance_A <= strong_max_A:
        return strong_points
    if distance_A <= medium_max_A:
        return medium_points
    return 0.0


def _score_tags(scores: dict[str, float], violation: bool, pocket_mode: str) -> list[str]:
    tags = []
    if scores.get("family_match", 0) > 0:
        tags.append("target_family")
    if scores.get("phenotype_match", 0) > 0:
        tags.append("target_phenotype")
    if scores.get("assay_readiness", 0) > 0:
        tags.append("assay_ready")
    if scores.get("retinal_pocket_strong", 0) > 0:
        tags.append("retinal_pocket_strong")
    if scores.get("retinal_pocket_medium", 0) > 0:
        tags.append("retinal_pocket_medium")
    if pocket_mode == "graded_fallback_reason":
        tags.append("pocket_graded_fallback")
    if violation:
        tags.append("protected_violation")
    return tags
