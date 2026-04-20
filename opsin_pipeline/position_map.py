from __future__ import annotations

import csv
import json
from pathlib import Path

from .ingest import load_scaffolds
from .schemas import Scaffold
from .structure.pocket import PocketMap, read_pocket_map


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
    "distance_to_retinal_A",
    "pocket_band",
    "review_needed",
    "review_status",
    "notes",
]


class UnmappedPocketMapError(ValueError):
    """Raised by pocket-annotate when the PocketMap has no seq_index on any residue."""


def write_draft_position_map(scaffolds: list[Scaffold], path: str | Path) -> Path:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for scaffold in scaffolds:
        for index, aa in enumerate(scaffold.sequence, start=1):
            rows.append(_blank_row(scaffold.name, scaffold.family, index, aa))
    _write_rows(output_path, rows)
    return output_path


def _blank_row(scaffold_name: str, family: str, seq_index: int, aa: str) -> dict[str, str]:
    return {
        "scaffold": scaffold_name,
        "family": family,
        "seq_index": str(seq_index),
        "aa": aa,
        "alignment_col": str(seq_index),
        "pdb_chain": "",
        "pdb_resnum": "",
        "region": "",
        "role": "",
        "protected": "false",
        "mutable": "false",
        "allowed_mutations": "",
        "distance_to_retinal_A": "",
        "pocket_band": "",
        "review_needed": "false",
        "review_status": "draft",
        "notes": "",
    }


def annotate_with_pocket(
    position_map_path: str | Path,
    pocket_map_path: str | Path,
    output_path: str | Path,
) -> Path:
    """Merge pocket evidence into a draft position-map CSV (Case A in spec §8.2).

    The PocketMap must carry ``seq_index`` on its residues (i.e. an offset or
    mapping was applied when it was created). If not, UnmappedPocketMapError is
    raised; the caller must re-run ``cli pocket`` with --pdb-offset or
    --pdb-mapping first.

    ``mutable`` and ``protected`` are never written by this step — they remain
    human-controlled.
    """
    pocket_map = read_pocket_map(pocket_map_path)
    pocket_by_index = {
        r.seq_index: r for r in pocket_map.pocket_residues if r.seq_index is not None
    }
    if not pocket_by_index:
        raise UnmappedPocketMapError(
            f"Pocket map {pocket_map_path} has no seq_index on any residue. "
            "Re-run `cli pocket` with --pdb-offset or --pdb-mapping so the "
            "evidence can be joined to scaffold rows; the unmapped evidence "
            "file stays in PDB numbering and is not merged into scaffold CSVs."
        )

    rows = _read_rows(position_map_path)
    for row in rows:
        if row.get("scaffold") != pocket_map.scaffold_name:
            continue
        try:
            seq_index = int(row["seq_index"])
        except (KeyError, ValueError):
            continue
        residue = pocket_by_index.get(seq_index)
        if residue is None:
            continue

        row["pdb_chain"] = residue.pdb_chain
        row["pdb_resnum"] = str(residue.pdb_resnum)
        row["distance_to_retinal_A"] = f"{residue.min_distance_A:.2f}"
        row["pocket_band"] = residue.band
        if residue.band in {"strong", "medium"}:
            row["region"] = _merge_region(row.get("region", ""), "retinal_pocket")
        if residue.role:
            row["role"] = _merge_region(row.get("role", ""), residue.role)
        if residue.mapping_note:
            row["review_needed"] = "true"
            row["notes"] = _append_note(row.get("notes", ""), residue.mapping_note)

    output_file = Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    _write_rows(output_file, rows)
    return output_file


def _merge_region(existing: str, addition: str) -> str:
    parts = [p.strip() for p in existing.replace(",", ";").split(";") if p.strip()]
    if addition not in parts:
        parts.append(addition)
    return "; ".join(parts)


def _append_note(existing: str, addition: str) -> str:
    if not existing:
        return addition
    if addition in existing:
        return existing
    return f"{existing}; {addition}"


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

