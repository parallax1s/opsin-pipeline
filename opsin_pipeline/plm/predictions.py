"""PLM prediction data model and JSON I/O per spec §4.

``PLMPrediction`` records one mutation's log-likelihood delta under a named
model. ``PLMPredictionSet`` bundles predictions per scaffold with a
sequence_sha256 for reproducibility.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class PLMPrediction:
    position: int                      # 1-indexed scaffold position
    from_aa: str                       # one-letter
    to_aa: str                         # one-letter
    log_likelihood_delta: float        # log P(to_aa|ctx) - log P(from_aa|ctx)
    model_id: str


@dataclass(frozen=True)
class PLMPredictionSet:
    scaffold_name: str
    sequence_sha256: str
    model_id: str
    predictions: list[PLMPrediction] = field(default_factory=list)


def sequence_sha256(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("utf-8")).hexdigest()


def predictions_to_dict(pset: PLMPredictionSet) -> dict:
    return {
        "scaffold_name": pset.scaffold_name,
        "sequence_sha256": pset.sequence_sha256,
        "model_id": pset.model_id,
        "predictions": [
            {
                "position": p.position,
                "from_aa": p.from_aa,
                "to_aa": p.to_aa,
                "log_likelihood_delta": p.log_likelihood_delta,
                "model_id": p.model_id,
            }
            for p in pset.predictions
        ],
    }


def predictions_from_dict(data: dict) -> PLMPredictionSet:
    model_id = str(data.get("model_id", ""))
    return PLMPredictionSet(
        scaffold_name=str(data["scaffold_name"]),
        sequence_sha256=str(data.get("sequence_sha256", "")),
        model_id=model_id,
        predictions=[
            PLMPrediction(
                position=int(p["position"]),
                from_aa=str(p["from_aa"]).upper(),
                to_aa=str(p["to_aa"]).upper(),
                log_likelihood_delta=float(p["log_likelihood_delta"]),
                model_id=str(p.get("model_id") or model_id),
            )
            for p in data.get("predictions", [])
        ],
    )


def write_predictions(pset: PLMPredictionSet, path: str | Path) -> Path:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(predictions_to_dict(pset), indent=2) + "\n", encoding="utf-8")
    return p


def read_predictions(path: str | Path) -> PLMPredictionSet:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    return predictions_from_dict(data)


def deltas_by_position(pset: PLMPredictionSet) -> dict[int, dict[str, float]]:
    """Group predictions into `{position: {to_aa: delta}}` for fast lookup during
    ingest merge. Later predictions for the same (position, to_aa) overwrite
    earlier ones (last-write-wins; usually redundant entries would be an upstream
    bug and tests should catch them)."""
    grouped: dict[int, dict[str, float]] = {}
    for p in pset.predictions:
        grouped.setdefault(p.position, {})[p.to_aa] = p.log_likelihood_delta
    return grouped
