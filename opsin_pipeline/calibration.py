from __future__ import annotations

import csv
import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .schemas import Candidate

Category = str  # "useful" | "neutral" | "disruptive"
MutationKey = frozenset  # frozenset[tuple[int, str]] — (position, to_aa)


@dataclass(frozen=True)
class CalibrationEntry:
    scaffold_name: str
    label: str
    mutations: MutationKey
    category: Category
    note: str = ""
    delta_lambda_nm: float | None = None
    source: str = ""


@dataclass(frozen=True)
class CalibrationReport:
    matched: int
    unmatched: int
    useful_total: int
    useful_matched: int
    disruptive_total: int
    disruptive_matched: int
    neutral_total: int
    neutral_matched: int
    auroc_useful_vs_disruptive: float | None
    mean_reciprocal_rank_useful: float | None
    top_k: int
    useful_in_top_k: int
    disruptive_in_top_k: int
    random_baseline_useful_in_top_k: float
    unmatched_labels: list[str] = field(default_factory=list)


_LABEL_TO_CATEGORY = {
    "useful": "useful",
    "positive": "useful",
    "good": "useful",
    "red_shift": "useful",
    "blue_shift": "useful",
    "neutral": "neutral",
    "none": "neutral",
    "null": "neutral",
    "disruptive": "disruptive",
    "negative": "disruptive",
    "bad": "disruptive",
    "loss_of_function": "disruptive",
}

_ID_SPLIT = re.compile(r"^(?P<scaffold>.+?)_(?P<muts>p\d+[A-Z]to[A-Z](?:__p\d+[A-Z]to[A-Z])*)$")
_MUTATION = re.compile(r"p(\d+)([A-Z])to([A-Z])")


def load_calibration(path: str | Path) -> list[CalibrationEntry]:
    """Load calibration entries from JSON or CSV.

    JSON: ``{"calibration_sets": [{"scaffold": str, "known_useful": [...],
    "known_neutral": [...], "known_disruptive": [...]}]}``. Each entry:
    ``{"label": str, "mutations": [{"position": int, "to": str}], ...}``.

    CSV: ``candidate_id, label[, evidence]``. Scaffold and mutation set are parsed
    out of the candidate_id. Labels accept aliases (``positive`` -> useful,
    ``negative``/``disruptive`` -> disruptive, ``neutral`` -> neutral).
    """
    suffix = Path(path).suffix.lower()
    if suffix == ".json":
        return _load_json(path)
    if suffix == ".csv":
        return _load_csv(path)
    raise ValueError(f"Calibration path must end in .json or .csv (got {suffix!r})")


def _load_json(path: str | Path) -> list[CalibrationEntry]:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    records = data.get("calibration_sets")
    if not isinstance(records, list):
        raise ValueError("Calibration JSON must contain a top-level 'calibration_sets' list")

    entries: list[CalibrationEntry] = []
    for record in records:
        scaffold = record.get("scaffold")
        if not scaffold:
            raise ValueError("Calibration set missing 'scaffold'")
        source = str(record.get("source", ""))
        for raw_category, category in (
            ("known_useful", "useful"),
            ("known_neutral", "neutral"),
            ("known_disruptive", "disruptive"),
        ):
            for item in record.get(raw_category, []):
                entries.append(
                    _parse_json_entry(
                        item, scaffold=str(scaffold), category=category, source=source
                    )
                )
    return entries


def _parse_json_entry(
    item: dict[str, Any], *, scaffold: str, category: Category, source: str
) -> CalibrationEntry:
    raw_mutations = item.get("mutations")
    if not raw_mutations:
        raise ValueError(f"Calibration entry in {scaffold} missing 'mutations'")
    mutations: MutationKey = frozenset(
        (int(m["position"]), str(m["to"]).upper()) for m in raw_mutations
    )
    delta = item.get("delta_lambda_nm")
    return CalibrationEntry(
        scaffold_name=scaffold,
        label=str(item.get("label", "")),
        mutations=mutations,
        category=category,
        note=str(item.get("note", "")),
        delta_lambda_nm=float(delta) if delta is not None else None,
        source=source,
    )


def _load_csv(path: str | Path) -> list[CalibrationEntry]:
    entries: list[CalibrationEntry] = []
    with Path(path).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or "candidate_id" not in reader.fieldnames:
            raise ValueError("Calibration CSV must have a 'candidate_id' column")
        if "label" not in reader.fieldnames:
            raise ValueError("Calibration CSV must have a 'label' column")
        for row in reader:
            candidate_id = (row.get("candidate_id") or "").strip()
            raw_label = (row.get("label") or "").strip().lower().replace("-", "_")
            if not candidate_id or not raw_label:
                continue
            category = _LABEL_TO_CATEGORY.get(raw_label)
            if category is None:
                raise ValueError(
                    f"Unrecognized calibration label {raw_label!r} for {candidate_id}"
                )
            scaffold, mutations = _parse_candidate_id(candidate_id)
            entries.append(
                CalibrationEntry(
                    scaffold_name=scaffold,
                    label=candidate_id,
                    mutations=mutations,
                    category=category,
                    note=(row.get("evidence") or row.get("note") or "").strip(),
                )
            )
    return entries


def _parse_candidate_id(candidate_id: str) -> tuple[str, MutationKey]:
    match = _ID_SPLIT.match(candidate_id)
    if match is None:
        raise ValueError(
            f"Candidate id {candidate_id!r} does not match expected format "
            "'<scaffold>_p<pos><from>to<to>[__p<pos><from>to<to>...]'"
        )
    scaffold = match.group("scaffold")
    mutations = frozenset(
        (int(pos), to) for pos, _from, to in _MUTATION.findall(match.group("muts"))
    )
    return scaffold, mutations


def candidate_key(candidate: Candidate) -> tuple[str, MutationKey]:
    return (
        candidate.scaffold_name,
        frozenset((m.position, m.to_aa) for m in candidate.mutations),
    )


def evaluate_ranking(
    ranked: list[Candidate],
    calibration: list[CalibrationEntry],
    top_k: int = 20,
) -> CalibrationReport:
    lookup: dict[tuple[str, MutationKey], CalibrationEntry] = {
        (entry.scaffold_name, entry.mutations): entry for entry in calibration
    }

    totals = {"useful": 0, "neutral": 0, "disruptive": 0}
    for entry in calibration:
        totals[entry.category] += 1

    matched_by_category: dict[Category, list[tuple[int, float]]] = {
        "useful": [],
        "neutral": [],
        "disruptive": [],
    }
    matched_keys: set[tuple[str, MutationKey]] = set()

    for rank_index, candidate in enumerate(ranked, start=1):
        key = candidate_key(candidate)
        entry = lookup.get(key)
        if entry is None:
            continue
        matched_keys.add(key)
        matched_by_category[entry.category].append(
            (rank_index, candidate.scores.get("total", 0.0))
        )

    unmatched_labels = sorted(
        entry.label or f"{entry.scaffold_name}:{sorted(entry.mutations)}"
        for entry in calibration
        if (entry.scaffold_name, entry.mutations) not in matched_keys
    )

    useful_scores = [s for _, s in matched_by_category["useful"]]
    disruptive_scores = [s for _, s in matched_by_category["disruptive"]]
    auroc = _auroc(useful_scores, disruptive_scores)

    useful_ranks = [r for r, _ in matched_by_category["useful"]]
    mrr = (
        round(sum(1.0 / r for r in useful_ranks) / len(useful_ranks), 4)
        if useful_ranks
        else None
    )

    useful_in_top_k = sum(1 for r in useful_ranks if r <= top_k)
    disruptive_in_top_k = sum(1 for r, _ in matched_by_category["disruptive"] if r <= top_k)
    total_ranked = len(ranked)
    baseline = (
        (totals["useful"] / total_ranked) * min(top_k, total_ranked)
        if total_ranked > 0
        else 0.0
    )

    matched = sum(len(v) for v in matched_by_category.values())
    return CalibrationReport(
        matched=matched,
        unmatched=len(calibration) - matched,
        useful_total=totals["useful"],
        useful_matched=len(matched_by_category["useful"]),
        disruptive_total=totals["disruptive"],
        disruptive_matched=len(matched_by_category["disruptive"]),
        neutral_total=totals["neutral"],
        neutral_matched=len(matched_by_category["neutral"]),
        auroc_useful_vs_disruptive=auroc,
        mean_reciprocal_rank_useful=mrr,
        top_k=top_k,
        useful_in_top_k=useful_in_top_k,
        disruptive_in_top_k=disruptive_in_top_k,
        random_baseline_useful_in_top_k=round(baseline, 3),
        unmatched_labels=unmatched_labels,
    )


def _auroc(positives: list[float], negatives: list[float]) -> float | None:
    if not positives or not negatives:
        return None
    total = len(positives) * len(negatives)
    wins = 0.0
    for pos in positives:
        for neg in negatives:
            if pos > neg:
                wins += 1.0
            elif pos == neg:
                wins += 0.5
    return round(wins / total, 4)
