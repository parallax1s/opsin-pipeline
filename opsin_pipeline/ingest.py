from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .schemas import MutablePosition, Scaffold


def load_scaffolds(path: str | Path) -> list[Scaffold]:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    records = data.get("scaffolds")
    if not isinstance(records, list):
        raise ValueError("Scaffold JSON must contain a top-level 'scaffolds' list")
    return [_parse_scaffold(record) for record in records]


def _parse_scaffold(record: dict[str, Any]) -> Scaffold:
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
        metadata=metadata,
    )

