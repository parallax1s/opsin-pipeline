from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from .schemas import Candidate


@dataclass(frozen=True)
class CalibrationResult:
    known_positive_count: int
    top_n: int
    positive_hits_at_n: int
    positive_recall_at_n: float


def evaluate_calibration(
    ranked_candidates: list[Candidate],
    calibration_path: str | Path,
    top_n: int = 10,
) -> CalibrationResult:
    positives = _load_positive_ids(calibration_path)
    top_ids = {candidate.candidate_id for candidate in ranked_candidates[:top_n]}
    hits = len(positives & top_ids)
    recall = hits / len(positives) if positives else 0.0
    return CalibrationResult(
        known_positive_count=len(positives),
        top_n=top_n,
        positive_hits_at_n=hits,
        positive_recall_at_n=round(recall, 3),
    )


def _load_positive_ids(path: str | Path) -> set[str]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        rows = csv.DictReader(handle)
        return {
            row["candidate_id"]
            for row in rows
            if row.get("label", "").strip().lower() == "positive"
        }
