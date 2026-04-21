from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path
from typing import Any

from .plm.predictions import PLMPredictionSet, deltas_by_position, read_predictions
from .schemas import MutablePosition, Scaffold
from .structure.pocket import PocketMap, read_pocket_map


def load_scaffolds(path: str | Path) -> list[Scaffold]:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    records = data.get("scaffolds")
    if not isinstance(records, list):
        raise ValueError("Scaffold JSON must contain a top-level 'scaffolds' list")
    scaffolds_path = Path(path)
    return [_parse_scaffold(record, base=scaffolds_path.parent) for record in records]


def _parse_scaffold(record: dict[str, Any], *, base: Path) -> Scaffold:
    required = ("name", "family", "sequence")
    missing = [key for key in required if not record.get(key)]
    if missing:
        raise ValueError(f"Scaffold record missing required fields: {', '.join(missing)}")

    mutable_positions = [
        MutablePosition(
            position=int(item["position"]),
            allowed=[str(aa) for aa in item.get("allowed", [])],
            reason=str(item.get("reason", "")),
        )
        for item in record.get("mutable_positions", [])
    ]

    pocket_map_path = record.get("pocket_map_path")
    if pocket_map_path:
        resolved = (base / pocket_map_path) if not Path(pocket_map_path).is_absolute() else Path(pocket_map_path)
        pocket_map = read_pocket_map(resolved)
        mutable_positions = _merge_pocket_map(mutable_positions, pocket_map)

    plm_predictions_path = record.get("plm_predictions_path")
    if plm_predictions_path:
        resolved = (
            (base / plm_predictions_path)
            if not Path(plm_predictions_path).is_absolute()
            else Path(plm_predictions_path)
        )
        pset = read_predictions(resolved)
        mutable_positions = _merge_plm_predictions(mutable_positions, pset)

    metadata = {
        key: value
        for key, value in record.items()
        if key
        not in {
            "name",
            "family",
            "sequence",
            "target_phenotypes",
            "assay_architectures",
            "protected_positions",
            "mutable_positions",
            "starting_lambda_nm",
            "pocket_map_path",
            "plm_predictions_path",
        }
    }

    return Scaffold(
        name=str(record["name"]),
        family=str(record["family"]),
        sequence=str(record["sequence"]).strip().upper(),
        target_phenotypes=[str(item) for item in record.get("target_phenotypes", [])],
        assay_architectures=[str(item) for item in record.get("assay_architectures", [])],
        protected_positions={int(item) for item in record.get("protected_positions", [])},
        mutable_positions=mutable_positions,
        starting_lambda_nm=(
            float(record["starting_lambda_nm"])
            if record.get("starting_lambda_nm") is not None
            else None
        ),
        pocket_map_path=str(pocket_map_path) if pocket_map_path else None,
        plm_predictions_path=str(plm_predictions_path) if plm_predictions_path else None,
        metadata=metadata,
    )


def _merge_pocket_map(
    positions: list[MutablePosition], pocket_map: PocketMap
) -> list[MutablePosition]:
    """Attach distance_to_retinal / role from the pocket map to each MutablePosition.

    Only pocket residues with a ``seq_index`` (i.e. mapping has been applied) are used.
    Positions without a matching pocket entry keep ``distance_to_retinal=None``.
    """
    by_seq_index = {
        r.seq_index: r for r in pocket_map.pocket_residues if r.seq_index is not None
    }
    merged: list[MutablePosition] = []
    for pos in positions:
        pocket_residue = by_seq_index.get(pos.position)
        if pocket_residue is None:
            merged.append(pos)
            continue
        merged.append(
            replace(
                pos,
                distance_to_retinal=pocket_residue.min_distance_A,
                role=pocket_residue.role,
            )
        )
    return merged


def _merge_plm_predictions(
    positions: list[MutablePosition], pset: PLMPredictionSet
) -> list[MutablePosition]:
    """Attach a per-position ``{to_aa: delta}`` dict from the PLM prediction set.

    Each position gets a filtered dict containing only the target AAs listed in
    that position's ``allowed`` list; deltas for unlisted AAs are dropped so the
    generator doesn't see stale predictions. Positions with no matching entries
    keep ``plm_log_likelihood_deltas=None`` rather than gain an empty dict — this
    lets the scorer distinguish "no PLM data" (fallback) from "PLM ran, nothing
    interesting" (would be weird but legal).
    """
    grouped = deltas_by_position(pset)
    merged: list[MutablePosition] = []
    for pos in positions:
        position_deltas = grouped.get(pos.position)
        if not position_deltas:
            merged.append(pos)
            continue
        filtered = {
            to_aa.upper(): position_deltas[to_aa.upper()]
            for to_aa in pos.allowed
            if to_aa.upper() in position_deltas
        }
        if not filtered:
            merged.append(pos)
            continue
        merged.append(replace(pos, plm_log_likelihood_deltas=filtered))
    return merged
