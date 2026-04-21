from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path
from typing import Any

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
