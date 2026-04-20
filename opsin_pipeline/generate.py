from __future__ import annotations

from itertools import combinations, product

from .schemas import Candidate, Mutation, Scaffold


def generate_single_mutants(scaffolds: list[Scaffold]) -> list[Candidate]:
    return generate_candidates(scaffolds, max_mutations=1)


def generate_candidates(scaffolds: list[Scaffold], max_mutations: int = 1) -> list[Candidate]:
    candidates: list[Candidate] = []
    for scaffold in scaffolds:
        mutation_options = _mutation_options(scaffold)
        limit = min(max_mutations, len(mutation_options))
        for mutation_count in range(1, limit + 1):
            for option_group in combinations(mutation_options, mutation_count):
                for mutation_group in product(*option_group):
                    mutations = list(mutation_group)
                    tags = sorted(
                        {
                            tag
                            for mutation in mutations
                            for tag in _tags_from_reason(mutation.reason)
                        }
                    )
                    candidate_id = f"{scaffold.name}_{'__'.join(m.label for m in mutations)}"
                    candidates.append(
                        Candidate(
                            candidate_id=candidate_id,
                            scaffold_name=scaffold.name,
                            family=scaffold.family,
                            mutations=mutations,
                            target_phenotypes=scaffold.target_phenotypes,
                            assay_architectures=scaffold.assay_architectures,
                            tags=tags,
                        )
                    )
    return candidates


def _mutation_options(scaffold: Scaffold) -> list[list[Mutation]]:
    options: list[list[Mutation]] = []
    for mutable in scaffold.mutable_positions:
        if mutable.position in scaffold.protected_positions:
            continue
        from_aa = scaffold.residue_at(mutable.position)
        mutations = [
            Mutation(
                position=mutable.position,
                from_aa=from_aa,
                to_aa=to_aa.upper(),
                reason=mutable.reason,
            )
            for to_aa in mutable.allowed
        ]
        if mutations:
            options.append(mutations)
    return options


def _tags_from_reason(reason: str) -> list[str]:
    text = reason.lower().replace("-", "_")
    tags: list[str] = []
    for tag in ("retinal_pocket", "spectral_tuning", "cGMP_activity", "photocurrent"):
        if tag.lower() in text:
            tags.append(tag)
    return tags
