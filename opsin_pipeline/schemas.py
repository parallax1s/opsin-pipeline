from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass(frozen=True)
class MutablePosition:
    position: int
    allowed: list[str]
    reason: str = ""


@dataclass(frozen=True)
class Mutation:
    position: int
    from_aa: str
    to_aa: str
    reason: str = ""

    @property
    def label(self) -> str:
        return f"p{self.position}{self.from_aa}to{self.to_aa}"


@dataclass(frozen=True)
class Scaffold:
    name: str
    family: str
    sequence: str
    target_phenotypes: list[str] = field(default_factory=list)
    assay_architectures: list[str] = field(default_factory=list)
    protected_positions: set[int] = field(default_factory=set)
    mutable_positions: list[MutablePosition] = field(default_factory=list)
    starting_lambda_nm: float | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def residue_at(self, position: int) -> str:
        if position < 1 or position > len(self.sequence):
            raise ValueError(
                f"Position {position} is outside {self.name} sequence length {len(self.sequence)}"
            )
        return self.sequence[position - 1]


@dataclass(frozen=True)
class Candidate:
    candidate_id: str
    scaffold_name: str
    family: str
    mutations: list[Mutation]
    target_phenotypes: list[str]
    assay_architectures: list[str]
    has_protected_violation: bool = False
    tags: list[str] = field(default_factory=list)
    scores: dict[str, float] = field(default_factory=dict)

    @property
    def mutation_summary(self) -> str:
        return ";".join(m.label for m in self.mutations)

    @property
    def reason_summary(self) -> str:
        reasons = [m.reason for m in self.mutations if m.reason]
        return "; ".join(reasons)

