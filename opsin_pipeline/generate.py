from __future__ import annotations

from .schemas import Candidate, Mutation, Scaffold


def generate_single_mutants(scaffolds: list[Scaffold]) -> list[Candidate]:
    candidates: list[Candidate] = []
    for scaffold in scaffolds:
        for mutable in scaffold.mutable_positions:
            if mutable.position in scaffold.protected_positions:
                continue
            from_aa = scaffold.residue_at(mutable.position)
            for to_aa in mutable.allowed:
                mutation = Mutation(
                    position=mutable.position,
                    from_aa=from_aa,
                    to_aa=to_aa.upper(),
                    reason=mutable.reason,
                )
                tags = _tags_from_reason(mutable.reason)
                candidates.append(
                    Candidate(
                        candidate_id=f"{scaffold.name}_{mutation.label}",
                        scaffold_name=scaffold.name,
                        family=scaffold.family,
                        mutations=[mutation],
                        target_phenotypes=scaffold.target_phenotypes,
                        assay_architectures=scaffold.assay_architectures,
                        tags=tags,
                    )
                )
    return candidates


def _tags_from_reason(reason: str) -> list[str]:
    text = reason.lower().replace("-", "_")
    tags: list[str] = []
    for tag in ("retinal_pocket", "spectral_tuning", "cGMP_activity", "photocurrent"):
        if tag.lower() in text:
            tags.append(tag)
    return tags

