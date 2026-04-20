from __future__ import annotations

from .schemas import Candidate


def diversify_ranked(
    candidates: list[Candidate],
    per_scaffold_cap: int | None = None,
    per_position_cap: int | None = None,
    top_n: int | None = None,
) -> list[Candidate]:
    """Greedy diversity filter over a pre-sorted ranked list.

    Walks candidates in the given order and admits each one unless it would exceed a cap.
    `per_position_cap` counts against every position a candidate touches, so a double
    mutant at (8, 19) consumes one slot at each of those positions.
    """
    selected: list[Candidate] = []
    scaffold_counts: dict[str, int] = {}
    position_counts: dict[tuple[str, int], int] = {}

    for candidate in candidates:
        if top_n is not None and len(selected) >= top_n:
            break

        if (
            per_scaffold_cap is not None
            and scaffold_counts.get(candidate.scaffold_name, 0) >= per_scaffold_cap
        ):
            continue

        if per_position_cap is not None and any(
            position_counts.get((candidate.scaffold_name, m.position), 0) >= per_position_cap
            for m in candidate.mutations
        ):
            continue

        selected.append(candidate)
        scaffold_counts[candidate.scaffold_name] = (
            scaffold_counts.get(candidate.scaffold_name, 0) + 1
        )
        for mutation in candidate.mutations:
            key = (candidate.scaffold_name, mutation.position)
            position_counts[key] = position_counts.get(key, 0) + 1

    return selected
