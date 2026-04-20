from __future__ import annotations

import csv
import json
from pathlib import Path

from .ingest import load_scaffolds
from .schemas import Scaffold


POSITION_MAP_FIELDS = [
    "scaffold",
    "family",
    "seq_index",
    "aa",
    "alignment_col",
    "pdb_chain",
    "pdb_resnum",
    "region",
    "role",
    "protected",
    "mutable",
    "allowed_mutations",
    "review_status",
    "notes",
]


def write_draft_position_map(scaffolds: list[Scaffold], path: str | Path) -> Path:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for scaffold in scaffolds:
        for index, aa in enumerate(scaffold.sequence, start=1):
            rows.append(
                {
                    "scaffold": scaffold.name,
                    "family": scaffold.family,
                    "seq_index": str(index),
                    "aa": aa,
                    "alignment_col": str(index),
                    "pdb_chain": "",
                    "pdb_resnum": "",
                    "region": "",
                    "role": "",
                    "protected": "false",
                    "mutable": "false",
                    "allowed_mutations": "",
                    "review_status": "draft",
                    "notes": "",
                }
            )
    _write_rows(output_path, rows)
    return output_path


def apply_position_map_to_scaffolds(
    scaffolds_path: str | Path,
    position_map_path: str | Path,
    output_path: str | Path,
) -> Path:
    source_data = json.loads(Path(scaffolds_path).read_text(encoding="utf-8"))
    rows = _read_rows(position_map_path)
    grouped = _group_reviewed_rows(rows)

    updated = []
    for scaffold in source_data.get("scaffolds", []):
        name = scaffold["name"]
        reviewed_rows = grouped.get(name, [])
        protected = sorted(
            int(row["seq_index"]) for row in reviewed_rows if _is_true(row["protected"])
        )
        mutable = []
        for row in reviewed_rows:
            if not _is_true(row["mutable"]):
                continue
            allowed = _split_allowed(row.get("allowed_mutations", ""))
            if not allowed:
                continue
            reason_parts = [
                row.get("region", ""),
                row.get("role", ""),
                row.get("notes", ""),
            ]
            reason = " ".join(part for part in reason_parts if part).strip()
            mutable.append(
                {
                    "position": int(row["seq_index"]),
                    "allowed": allowed,
                    "reason": reason,
                }
            )

        next_scaffold = dict(scaffold)
        next_scaffold["protected_positions"] = protected
        next_scaffold["mutable_positions"] = mutable
        updated.append(next_scaffold)

    output = {"scaffolds": updated}
    output_file = Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    output_file.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    return output_file


def _group_reviewed_rows(rows: list[dict[str, str]]) -> dict[str, list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        status = row.get("review_status", "").strip().lower()
        if status not in {"reviewed", "approved"}:
            continue
        grouped.setdefault(row["scaffold"], []).append(row)
    return grouped


def _split_allowed(value: str) -> list[str]:
    return [item.strip().upper() for item in value.replace(",", ";").split(";") if item.strip()]


def _is_true(value: str) -> bool:
    return value.strip().lower() in {"1", "true", "yes", "y"}


def _read_rows(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _write_rows(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=POSITION_MAP_FIELDS)
        writer.writeheader()
        writer.writerows(rows)

