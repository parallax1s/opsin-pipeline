from __future__ import annotations

import argparse
from pathlib import Path

from .calibration import evaluate_ranking, load_calibration
from .generate import generate_candidates
from .ingest import load_scaffolds
from .plm.predictions import write_predictions
from .plm.scorer import MockPLMScorer, MutationRequest, STANDARD_AAS
from .position_map import (
    UnmappedPocketMapError,
    annotate_with_pocket,
    apply_position_map_to_scaffolds,
    write_draft_position_map,
)
from .report import write_candidate_csv, write_decision_report
from .score import rank_candidates
from .structure.ligands import identify_ligands
from .structure.pdb import parse_pdb
from .structure.pocket import (
    DEFAULT_CUTOFF_A,
    DEFAULT_MEDIUM_MAX_A,
    DEFAULT_STRONG_MAX_A,
    apply_offset_mapping,
    compute_pocket,
    write_pocket_map,
)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the local opsin MVP pipeline")
    subparsers = parser.add_subparsers(dest="command")

    _add_run_parser(subparsers)
    _add_pocket_parser(subparsers)
    _add_pocket_annotate_parser(subparsers)
    _add_plm_parser(subparsers)
    _add_draft_map_parser(subparsers)
    _add_apply_map_parser(subparsers)

    args = parser.parse_args()

    if args.command == "draft-position-map":
        scaffolds = load_scaffolds(args.scaffolds)
        path = write_draft_position_map(scaffolds, args.out)
        print(f"Wrote {path}")
        return

    if args.command == "apply-position-map":
        path = apply_position_map_to_scaffolds(args.scaffolds, args.position_map, args.out)
        print(f"Wrote {path}")
        return

    if args.command == "pocket":
        _run_pocket(args)
        return

    if args.command == "pocket-annotate":
        _run_pocket_annotate(args)
        return

    if args.command == "plm":
        _run_plm(args)
        return

    if args.command is None:
        parser.error("a command is required")

    _run_pipeline(args)


def _add_run_parser(subparsers: argparse._SubParsersAction) -> None:
    run_parser = subparsers.add_parser("run", help="Generate, score, and report candidates")
    run_parser.add_argument("--scaffolds", required=True, help="Path to scaffold JSON")
    run_parser.add_argument("--out", required=True, help="Output directory")
    run_parser.add_argument("--target-family", default=None)
    run_parser.add_argument("--target-phenotype", default=None)
    run_parser.add_argument("--top", type=int, default=10, help="Top-N size for the report")
    run_parser.add_argument(
        "--max-mutations",
        type=int,
        default=1,
        help="Maximum number of simultaneous mutations per candidate (1 = singles only).",
    )
    run_parser.add_argument(
        "--max-combinations-per-scaffold",
        type=int,
        default=None,
        help="Hard cap on candidates per scaffold; truncates higher-order mutants first.",
    )
    run_parser.add_argument(
        "--per-scaffold-cap",
        type=int,
        default=None,
        help="In the report's top-N, limit how many entries any one scaffold contributes.",
    )
    run_parser.add_argument(
        "--per-position-cap",
        type=int,
        default=None,
        help="In the report's top-N, limit how many entries touch any one (scaffold, position).",
    )
    run_parser.add_argument(
        "--calibration",
        default=None,
        help="Optional calibration file (.json or .csv) with known useful / disruptive mutations.",
    )
    run_parser.add_argument(
        "--calibration-top-k",
        type=int,
        default=20,
        help="Top-k used by the calibration evaluator.",
    )
    run_parser.add_argument(
        "--graded-pocket",
        action="store_true",
        help=(
            "Use structure-grounded graded pocket scoring (bands strong <=4 A / "
            "medium <=5.5 A / none) instead of the legacy reason-string rule. "
            "Requires scaffolds with a pocket_map_path. Default OFF until the "
            "literature calibration set demonstrates it improves AUROC."
        ),
    )


def _add_pocket_parser(subparsers: argparse._SubParsersAction) -> None:
    pocket_parser = subparsers.add_parser(
        "pocket",
        help="Parse a PDB and write a PocketMap JSON for review / scaffold merge",
    )
    pocket_parser.add_argument("--scaffold", required=True, help="Scaffold name (goes into PocketMap)")
    pocket_parser.add_argument("--pdb", required=True, help="Path to a PDB file")
    pocket_parser.add_argument("--out", required=True, help="Output JSON path")
    pocket_parser.add_argument(
        "--pdb-chain",
        default=None,
        help="Chain to analyze (ligand + protein residues). Required if the PDB has multiple chains with retinal.",
    )
    pocket_parser.add_argument(
        "--pdb-offset",
        type=int,
        default=None,
        help="Apply seq_index = pdb_resnum + offset to every pocket residue.",
    )
    pocket_parser.add_argument(
        "--scaffolds",
        default=None,
        help="Optional scaffold JSON for AA-identity verification when --pdb-offset is given.",
    )
    pocket_parser.add_argument(
        "--strict-mapping",
        action="store_true",
        help="With --pdb-offset + --scaffolds, promote AA mismatches from warnings to hard errors.",
    )
    pocket_parser.add_argument(
        "--cutoff",
        type=float,
        default=DEFAULT_CUTOFF_A,
        help=f"Outer radius for pocket residue inclusion (A). Default {DEFAULT_CUTOFF_A}.",
    )
    pocket_parser.add_argument(
        "--strong-max",
        type=float,
        default=DEFAULT_STRONG_MAX_A,
        help=f"Upper bound of the 'strong' band (A). Default {DEFAULT_STRONG_MAX_A}.",
    )
    pocket_parser.add_argument(
        "--medium-max",
        type=float,
        default=DEFAULT_MEDIUM_MAX_A,
        help=f"Upper bound of the 'medium' band (A). Default {DEFAULT_MEDIUM_MAX_A}.",
    )
    pocket_parser.add_argument(
        "--allow-weird-ligand",
        action="store_true",
        help="Skip the 15-25 heavy-atom sanity check for the retinal ligand.",
    )


def _add_plm_parser(subparsers: argparse._SubParsersAction) -> None:
    plm_parser = subparsers.add_parser(
        "plm",
        help="Produce a PLMPredictionSet JSON for a scaffold's mutable positions",
    )
    plm_parser.add_argument("--scaffold", required=True, help="Scaffold name to score")
    plm_parser.add_argument("--scaffolds", required=True, help="Scaffold JSON file")
    plm_parser.add_argument("--out", required=True, help="Output JSON path")
    plm_parser.add_argument(
        "--mock",
        action="store_true",
        help="Use the deterministic MockPLMScorer (no torch, no network).",
    )
    plm_parser.add_argument(
        "--mock-seed",
        type=int,
        default=0,
        help="Seed for MockPLMScorer (default 0).",
    )
    plm_parser.add_argument(
        "--model",
        default="esm2_t12_35M_UR50D",
        help=(
            "Model identifier for the real scorer (default esm2_t12_35M_UR50D). "
            "Ignored when --mock is passed."
        ),
    )


def _add_pocket_annotate_parser(subparsers: argparse._SubParsersAction) -> None:
    annotate = subparsers.add_parser(
        "pocket-annotate",
        help="Merge PocketMap evidence into an existing draft position-map CSV",
    )
    annotate.add_argument("--position-map", required=True, help="Input draft position-map CSV")
    annotate.add_argument("--pocket-map", required=True, help="PocketMap JSON (must have seq_index)")
    annotate.add_argument("--out", required=True, help="Output CSV path")


def _add_draft_map_parser(subparsers: argparse._SubParsersAction) -> None:
    draft_parser = subparsers.add_parser("draft-position-map", help="Write a draft review CSV")
    draft_parser.add_argument("--scaffolds", required=True, help="Path to scaffold JSON")
    draft_parser.add_argument("--out", required=True, help="Output CSV path")


def _add_apply_map_parser(subparsers: argparse._SubParsersAction) -> None:
    apply_parser = subparsers.add_parser(
        "apply-position-map",
        help="Apply reviewed positions to a scaffold JSON file",
    )
    apply_parser.add_argument("--scaffolds", required=True, help="Source scaffold JSON")
    apply_parser.add_argument("--position-map", required=True, help="Reviewed position CSV")
    apply_parser.add_argument("--out", required=True, help="Output scaffold JSON")


def _run_pocket(args: argparse.Namespace) -> None:
    pdb_path = Path(args.pdb)
    atoms = parse_pdb(pdb_path)
    matches = identify_ligands(
        atoms,
        preferred_chain=args.pdb_chain,
        allow_weird_ligand=args.allow_weird_ligand,
        path_hint=str(pdb_path),
    )
    pocket_map = compute_pocket(
        atoms,
        matches,
        scaffold_name=args.scaffold,
        pdb_path=str(pdb_path),
        pdb_chain=args.pdb_chain,
        cutoff_A=args.cutoff,
        strong_max_A=args.strong_max,
        medium_max_A=args.medium_max,
    )

    if args.pdb_offset is not None:
        scaffold_sequence: str | None = None
        if args.scaffolds:
            scaffold_by_name = {s.name: s for s in load_scaffolds(args.scaffolds)}
            if args.scaffold not in scaffold_by_name:
                raise SystemExit(
                    f"Scaffold {args.scaffold!r} not found in {args.scaffolds}; cannot verify AA identity."
                )
            scaffold_sequence = scaffold_by_name[args.scaffold].sequence
        pocket_map = apply_offset_mapping(
            pocket_map,
            offset=args.pdb_offset,
            scaffold_sequence=scaffold_sequence,
            strict=args.strict_mapping,
        )

    out_path = write_pocket_map(pocket_map, args.out)
    print(f"Wrote {out_path}")
    print(
        f"  {len(pocket_map.pocket_residues)} pocket residues, "
        f"{sum(1 for r in pocket_map.pocket_residues if r.band == 'strong')} strong, "
        f"{sum(1 for r in pocket_map.pocket_residues if r.band == 'medium')} medium."
    )
    mismatches = [r for r in pocket_map.pocket_residues if r.mapping_note]
    if mismatches:
        print(f"  {len(mismatches)} mapping mismatches — review_needed")
        for r in mismatches[:5]:
            print(f"    {r.mapping_note}")


def _run_plm(args: argparse.Namespace) -> None:
    scaffolds = load_scaffolds(args.scaffolds)
    scaffold_by_name = {s.name: s for s in scaffolds}
    if args.scaffold not in scaffold_by_name:
        raise SystemExit(
            f"plm: scaffold {args.scaffold!r} not found in {args.scaffolds}. "
            f"Available: {sorted(scaffold_by_name)}"
        )
    scaffold = scaffold_by_name[args.scaffold]

    if len(scaffold.sequence) > 1024:
        raise SystemExit(
            f"plm: scaffold sequence length {len(scaffold.sequence)} exceeds the "
            "1024-residue ESM2 context window. Long-sequence handling is deferred "
            "per spec §12b."
        )

    # Build mutation requests from mutable_positions. Skip silent ones so downstream
    # JSON doesn't carry pointless zero deltas.
    requests: list[MutationRequest] = []
    for mutable in scaffold.mutable_positions:
        from_aa = scaffold.residue_at(mutable.position)
        if from_aa not in STANDARD_AAS:
            raise SystemExit(
                f"plm: scaffold {scaffold.name!r} has non-standard residue {from_aa!r} "
                f"at position {mutable.position}; only the standard 20 amino acids "
                "are supported (spec §9 hard error)."
            )
        for to_aa in mutable.allowed:
            to_upper = to_aa.upper()
            if to_upper == from_aa:
                continue
            if to_upper not in STANDARD_AAS:
                raise SystemExit(
                    f"plm: target residue {to_upper!r} at position {mutable.position} "
                    "is not one of the standard 20 amino acids."
                )
            requests.append(MutationRequest(position=mutable.position, from_aa=from_aa, to_aa=to_upper))

    if not requests:
        print(
            f"plm: scaffold {scaffold.name} has no non-silent mutations in its "
            "mutable_positions; writing empty PLMPredictionSet."
        )

    if args.mock:
        scorer = MockPLMScorer(seed=args.mock_seed)
    else:
        raise SystemExit(
            "plm: only --mock is wired in this chunk. The ESM2 backend "
            "(opsin_pipeline/plm/esm.py) is a follow-up commit gated on "
            "an optional torch install."
        )

    pset = scorer.score(
        scaffold_name=scaffold.name,
        sequence=scaffold.sequence,
        mutations=requests,
    )
    out_path = write_predictions(pset, args.out)
    print(f"Wrote {out_path}")
    print(
        f"  {len(pset.predictions)} predictions via {pset.model_id}; "
        f"median delta {_median([p.log_likelihood_delta for p in pset.predictions])}"
    )


def _median(values: list[float]) -> float:
    if not values:
        return 0.0
    sorted_vals = sorted(values)
    mid = len(sorted_vals) // 2
    if len(sorted_vals) % 2:
        return round(sorted_vals[mid], 4)
    return round((sorted_vals[mid - 1] + sorted_vals[mid]) / 2, 4)


def _run_pocket_annotate(args: argparse.Namespace) -> None:
    try:
        out_path = annotate_with_pocket(
            position_map_path=args.position_map,
            pocket_map_path=args.pocket_map,
            output_path=args.out,
        )
    except UnmappedPocketMapError as exc:
        raise SystemExit(f"pocket-annotate: {exc}")
    print(f"Wrote {out_path}")


def _run_pipeline(args: argparse.Namespace) -> None:
    scaffolds = load_scaffolds(args.scaffolds)
    candidates, gen_stats = generate_candidates(
        scaffolds,
        max_mutations=args.max_mutations,
        max_combinations_per_scaffold=args.max_combinations_per_scaffold,
    )
    ranked = rank_candidates(
        candidates,
        scaffolds,
        target_family=args.target_family,
        target_phenotype=args.target_phenotype,
        use_graded_pocket=args.graded_pocket,
    )

    calibration_report = None
    if args.calibration:
        entries = load_calibration(args.calibration)
        calibration_report = evaluate_ranking(
            ranked, entries, top_k=args.calibration_top_k
        )

    out_dir = Path(args.out)
    csv_path = write_candidate_csv(ranked, scaffolds, out_dir / "ranked_candidates.csv")
    report_path = write_decision_report(
        ranked,
        scaffolds,
        out_dir / "decision_report.md",
        target_family=args.target_family,
        target_phenotype=args.target_phenotype,
        top_n=args.top,
        per_scaffold_cap=args.per_scaffold_cap,
        per_position_cap=args.per_position_cap,
        generation_stats=gen_stats,
        calibration_report=calibration_report,
        use_graded_pocket=args.graded_pocket,
    )
    print(f"Wrote {csv_path}")
    print(f"Wrote {report_path}")
    if gen_stats.per_scaffold_truncated:
        for name, dropped in sorted(gen_stats.per_scaffold_truncated.items()):
            print(f"  truncated {dropped} candidates from {name} (cap)")
    if calibration_report is not None:
        auroc = calibration_report.auroc_useful_vs_disruptive
        print(
            "  calibration: useful {u}/{ut}, disruptive {d}/{dt}, AUROC {a}".format(
                u=calibration_report.useful_matched,
                ut=calibration_report.useful_total,
                d=calibration_report.disruptive_matched,
                dt=calibration_report.disruptive_total,
                a=auroc if auroc is not None else "n/a",
            )
        )


if __name__ == "__main__":
    main()
