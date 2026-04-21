"""PLM scorer interface + MockPLMScorer per spec §5.

The abstract ``PLMScorer`` defines the minimal contract: given a scaffold
sequence and a set of (position, from_aa, to_aa) mutations, return a
``PLMPredictionSet``. Concrete backends:

- ``MockPLMScorer`` (this module) — deterministic pseudo-random deltas keyed
  on ``(position, from_aa, to_aa, seed)``. Zero dependencies. Powers all
  unit tests and the ``--mock`` CLI path.
- ``ESM2Scorer`` (opsin_pipeline/plm/esm.py — to land in a later chunk)
  wraps the HuggingFace ``esm2_t12_35M_UR50D`` model behind a lazy import.

Any future backend (SaProt, ESM-C, fine-tuned) subclasses ``PLMScorer`` and
wires a new ``--model`` value in the CLI.
"""
from __future__ import annotations

import hashlib
from abc import ABC, abstractmethod
from dataclasses import dataclass, field

from .predictions import PLMPrediction, PLMPredictionSet, sequence_sha256


STANDARD_AAS: frozenset[str] = frozenset("ACDEFGHIKLMNPQRSTVWY")


@dataclass(frozen=True)
class MutationRequest:
    position: int                     # 1-indexed
    from_aa: str                      # one-letter, uppercase, expected to match sequence
    to_aa: str                        # one-letter, uppercase


class PLMScorer(ABC):
    """Abstract base for PLM log-likelihood scoring of point mutations."""

    model_id: str = "unknown"

    @abstractmethod
    def score(
        self,
        *,
        scaffold_name: str,
        sequence: str,
        mutations: list[MutationRequest],
    ) -> PLMPredictionSet:
        """Return a PLMPredictionSet containing one PLMPrediction per input mutation.

        Implementations must:

        - verify that sequence[mutation.position - 1] == mutation.from_aa (raise if not).
        - verify that mutation.to_aa is in STANDARD_AAS (raise if not).
        - leave silent mutations (from_aa == to_aa) as 0.0 deltas.
        - raise a suitable error for out-of-range positions.
        """


class UnsupportedResidueError(ValueError):
    """Raised when a scaffold has a non-standard residue at a requested position."""


class PositionOutOfRangeError(ValueError):
    """Raised when the mutation position is outside [1, len(sequence)]."""


class WildTypeMismatchError(ValueError):
    """Raised when the request's from_aa doesn't match the scaffold sequence."""


class MockPLMScorer(PLMScorer):
    """Deterministic pseudo-random scorer for tests.

    The delta for ``(scaffold, position, from_aa, to_aa, seed)`` is stable
    across runs and across platforms. Range is roughly [-5, +2] centred near
    -1.5, which is the rough shape of stock ESM2 deltas on opsin positions.
    """

    model_id = "mock"

    def __init__(self, seed: int = 0):
        self.seed = seed

    def score(
        self,
        *,
        scaffold_name: str,
        sequence: str,
        mutations: list[MutationRequest],
    ) -> PLMPredictionSet:
        predictions: list[PLMPrediction] = []
        for req in mutations:
            self._validate(sequence, req)
            delta = 0.0 if req.from_aa == req.to_aa else self._delta_for(scaffold_name, req)
            predictions.append(
                PLMPrediction(
                    position=req.position,
                    from_aa=req.from_aa,
                    to_aa=req.to_aa,
                    log_likelihood_delta=delta,
                    model_id=self.model_id,
                )
            )
        return PLMPredictionSet(
            scaffold_name=scaffold_name,
            sequence_sha256=sequence_sha256(sequence),
            model_id=self.model_id,
            predictions=predictions,
        )

    def _validate(self, sequence: str, req: MutationRequest) -> None:
        if req.position < 1 or req.position > len(sequence):
            raise PositionOutOfRangeError(
                f"position {req.position} is outside sequence length {len(sequence)}"
            )
        actual_wt = sequence[req.position - 1].upper()
        if actual_wt not in STANDARD_AAS:
            raise UnsupportedResidueError(
                f"scaffold residue at position {req.position} is {actual_wt!r}; "
                "only the standard 20 amino acids are supported"
            )
        if req.to_aa not in STANDARD_AAS:
            raise UnsupportedResidueError(
                f"target residue {req.to_aa!r} at position {req.position} is not one of the standard 20"
            )
        if req.from_aa != actual_wt:
            raise WildTypeMismatchError(
                f"mutation request says from_aa={req.from_aa} at position {req.position}, "
                f"but scaffold sequence has {actual_wt}"
            )

    def _delta_for(self, scaffold_name: str, req: MutationRequest) -> float:
        # Deterministic hash -> bounded float. Zero dependencies, works offline.
        key = f"{scaffold_name}:{req.position}:{req.from_aa}:{req.to_aa}:{self.seed}".encode()
        digest = hashlib.sha256(key).digest()
        # Take the first 4 bytes as an unsigned int and map to [-5.0, +2.0].
        raw = int.from_bytes(digest[:4], "big")
        unit = raw / 0xFFFFFFFF  # in [0, 1]
        return round(-5.0 + unit * 7.0, 4)
