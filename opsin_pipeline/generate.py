from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations, product

from .schemas import Candidate, MutablePosition, Mutation, Scaffold


@dataclass(frozen=True)
class GenerationStats:
    total_generated: int
    per_scaffold_generated: dict[str, int] = field(default_factory=dict)
    per_scaffold_truncated: dict[str, int] = field(default_factory=dict)


def generate_candidates(
    scaffolds: list[Scaffold],
    max_mutations: int = 1,
    max_combinations_per_scaffold: int | None = None,
) -> tuple[list[Candidate], GenerationStats]:
    if max_mutations < 1:
        raise ValueError("max_mutations must be >= 1")

    candidates: list[Candidate] = []
    per_generated: dict[str, int] = {}
    per_truncated: dict[str, int] = {}

    for scaffold in scaffolds:
        usable = sorted(
            (
                mutable
                for mutable in scaffold.mutable_positions
                if mutable.position not in scaffold.protected_positions
            ),
            key=lambda m: m.position,
        )

        scaffold_candidates: list[Candidate] = []
        limit = min(max_mutations, len(usable))
        # walk Hamming levels in order so the cap keeps singles before doubles, etc.
        for mutation_count in range(1, limit + 1):
            for option_group in combinations(usable, mutation_count):
                aa_choices = [_normalized_allowed(m) for m in option_group]
                for aa_combo in product(*aa_choices):
                    mutations = _build_mutations(scaffold, option_group, aa_combo)
                    if len(mutations) != mutation_count:
                        continue  # one or more positions collapsed to silent (to_aa == from_aa)
                    tags = _tags_from_reasons(m.reason for m in option_group)
                    label = "__".join(mutation.label for mutation in mutations)
                    scaffold_candidates.append(
                        Candidate(
                            candidate_id=f"{scaffold.name}_{label}",
                            scaffold_name=scaffold.name,
                            family=scaffold.family,
                            mutations=mutations,
                            target_phenotypes=scaffold.target_phenotypes,
                            assay_architectures=scaffold.assay_architectures,
                            tags=tags,
                        )
                    )

        total = len(scaffold_candidates)
        if (
            max_combinations_per_scaffold is not None
            and total > max_combinations_per_scaffold
        ):
            per_truncated[scaffold.name] = total - max_combinations_per_scaffold
            scaffold_candidates = scaffold_candidates[:max_combinations_per_scaffold]

        per_generated[scaffold.name] = len(scaffold_candidates)
        candidates.extend(scaffold_candidates)

    stats = GenerationStats(
        total_generated=len(candidates),
        per_scaffold_generated=per_generated,
        per_scaffold_truncated=per_truncated,
    )
    return candidates, stats


def generate_single_mutants(scaffolds: list[Scaffold]) -> list[Candidate]:
    candidates, _ = generate_candidates(scaffolds, max_mutations=1)
    return candidates


def _normalized_allowed(mutable: MutablePosition) -> list[str]:
    return sorted({aa.upper() for aa in mutable.allowed})


def _build_mutations(
    scaffold: Scaffold,
    option_group: tuple[MutablePosition, ...],
    aa_combo: tuple[str, ...],
) -> list[Mutation]:
    mutations: list[Mutation] = []
    for mutable, to_aa in zip(option_group, aa_combo):
        from_aa = scaffold.residue_at(mutable.position)
        if to_aa == from_aa:
            continue
        mutations.append(
            Mutation(
                position=mutable.position,
                from_aa=from_aa,
                to_aa=to_aa,
                reason=mutable.reason,
                distance_to_retinal=mutable.distance_to_retinal,
                role=mutable.role,
            )
        )
    return mutations


def _tags_from_reasons(reasons) -> list[str]:
    tags: set[str] = set()
    for reason in reasons:
        text = reason.lower().replace("-", "_")
        for tag in ("retinal_pocket", "spectral_tuning", "cGMP_activity", "photocurrent"):
            if tag.lower() in text:
                tags.add(tag)
    return sorted(tags)
