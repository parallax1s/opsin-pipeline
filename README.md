# Opsin Pipeline MVP

Local, dependency-free MVP for broad in silico opsin candidate triage.

The package does not run LigandMPNN, Caliby, RhoMax, METL, or QM/MM directly yet. It defines the data formats, candidate generation, transparent scoring, ranking, diversity filtering, and calibration hooks those adapters can plug into later.

## Run Tests

```bash
PYTHONPATH=. python3 -m unittest discover -s tests -v
```

## Run Example

Singles only (original MVP behavior):

```bash
PYTHONPATH=. python3 -m opsin_pipeline.cli run \
  --scaffolds configs/example_scaffolds.json \
  --out out/example \
  --target-family RhGC \
  --target-phenotype spectral_tuning
```

Singles + doubles with diversity caps and a calibration check (JSON):

```bash
PYTHONPATH=. python3 -m opsin_pipeline.cli run \
  --scaffolds configs/example_scaffolds.json \
  --calibration configs/example_calibration.json \
  --out out/example_multi \
  --target-family RhGC \
  --target-phenotype spectral_tuning \
  --max-mutations 2 \
  --max-combinations-per-scaffold 500 \
  --per-scaffold-cap 5 \
  --per-position-cap 2 \
  --top 20
```

The same run with a CSV calibration file (the format Codex added):

```bash
PYTHONPATH=. python3 -m opsin_pipeline.cli run \
  --scaffolds configs/example_scaffolds.json \
  --calibration configs/example_calibration.csv \
  --out out/example_csv \
  --target-family RhGC \
  --target-phenotype spectral_tuning \
  --max-mutations 2 \
  --per-scaffold-cap 5 \
  --top 20
```

Outputs:

- `ranked_candidates.csv` — full ranked list (never truncated by diversity caps)
- `decision_report.md` — diversified top-N, scaffold summary, generation stats, calibration section

## Real Scaffold Seed Set

`configs/real_scaffolds_seed.json` and `configs/real_scaffolds_seed.fasta` contain initial scaffold sequences pulled from UniProt and RCSB:

- ChR1 and ChR2 full UniProt sequences
- CaRhGC and BeRhGC1 UniProt sequences
- BeCNG1 as a cGMP-gated channel partner candidate
- PDB-backed ChR2, C1C2, Chrimson, GtACR1, and ChRmine structure sequences

The real seed config intentionally leaves `mutable_positions` empty. Add protected and mutable positions only after sequence/structure numbering is reviewed.

## Review Alignment / Numbering

Start by creating a draft position map:

```bash
PYTHONPATH=. python3 -m opsin_pipeline.cli draft-position-map \
  --scaffolds configs/real_scaffolds_seed.json \
  --out out/review/position_map_draft.csv
```

Review the CSV manually or with an external alignment workflow. Mark rows as `review_status=reviewed`, then set:

- `protected=true` for residues that should not be mutated
- `mutable=true` plus `allowed_mutations=F;Y;W` for reviewed mutation sites
- `region`, `role`, and `notes` for traceability

Apply the reviewed map back into scaffold JSON:

```bash
PYTHONPATH=. python3 -m opsin_pipeline.cli apply-position-map \
  --scaffolds configs/real_scaffolds_seed.json \
  --position-map out/review/position_map_draft.csv \
  --out configs/real_scaffolds_reviewed.json
```

Candidate generation should use the reviewed scaffold JSON, not the seed file.

## Candidate generation

`generate_candidates(scaffolds, max_mutations=N, max_combinations_per_scaffold=K)`:

- Enumerates combinations of mutable positions up to `max_mutations` (default 1 = singles).
- Skips protected positions and silent (wild-type → wild-type) mutations.
- Orders candidates by Hamming level then by position — a per-scaffold cap truncates higher-order mutants first and reports the drop count via `GenerationStats`.

## Diversity-aware reporting

`diversify_ranked(candidates, per_scaffold_cap, per_position_cap, top_n)` walks the ranked list greedily and admits each candidate unless it would exceed a cap. `per_position_cap` counts every position a candidate touches, so a double at (8, 19) uses one slot at each.

The full CSV still contains every ranked candidate; the diversity cap only shapes the report's shortlist.

## Calibration

Two accepted input formats, both routed through the same evaluator.

### JSON (structured)

```json
{
  "calibration_sets": [
    {
      "scaffold": "CaRhGC",
      "source": "DOI or internal notebook",
      "known_useful":     [{"label": "K8F", "mutations": [{"position": 8, "to": "F"}], "delta_lambda_nm": -12}],
      "known_neutral":    [{"label": "...", "mutations": [...]}],
      "known_disruptive": [{"label": "...", "mutations": [...]}]
    }
  ]
}
```

### CSV (spreadsheet-friendly)

```csv
candidate_id,label,evidence
CaRhGC_p8KtoF,positive,DOI 10.xxxx
CaRhGC_p8KtoW,negative,steric clash
```

Labels accept aliases: `positive`/`useful`/`good` → useful, `negative`/`disruptive`/`bad` → disruptive, `neutral`/`none` → neutral. Scaffold name and mutation set are parsed out of `candidate_id`, so multi-mutant rows (`CaRhGC_p8KtoF__p19FtoA,positive,...`) match correctly.

`evaluate_ranking(ranked, entries, top_k=20)` matches entries against ranked candidates by `(scaffold, frozenset((position, to_aa)))` — naming-independent. Metrics:

- `useful_matched / useful_total` and the same for disruptive/neutral.
- `auroc_useful_vs_disruptive` — Mann-Whitney U over total score; ties counted as 0.5.
- `mean_reciprocal_rank_useful` — 1/rank averaged over useful matches.
- `useful_in_top_k`, `disruptive_in_top_k`, plus `random_baseline_useful_in_top_k` so you can tell whether the ranking beats chance.
- `unmatched_labels` — calibration entries that did not appear in the ranked list (often because they hit protected positions and were never generated).

`configs/example_calibration.json` and `configs/example_calibration.csv` are synthetic — replace with a curated set (BR D85/T89, ChR C128/H134, RhGC literature shifts, etc.) before trusting the numbers.

## Scoring signals

- target family match
- target phenotype match
- retinal-pocket reason tag
- assay-readiness tag (`growth_selection` or `biochemical`)
- mutation-count penalty: −0.2 per mutation
- protected-residue penalty: −5 if any mutation hits a protected position. The flag is re-derived from the scaffold at scoring time, so externally supplied candidates (calibration entries, adapter outputs) are flagged too.
- optional multi-mutant generation (`--max-mutations`)
- optional diversity caps (`--per-scaffold-cap`, `--per-position-cap`)
- optional calibration check (`--calibration`)

Treat the total as a triage heuristic — a reproducible queue for deeper structure, physics, and wet-lab validation. The calibration harness is the honest way to tell whether changes to the score function actually improve the ranking.

## Roadmap

- Structure-grounded retinal-pocket mask (parse PDB/mmCIF, identify retinal, flag residues within 4–5 Å) — replaces the reason-string regex with geometry.
- Optional PLM adapter (ESM / ESM-C / SaProt) as a sequence-plausibility prior.
- Model adapter layer: LigandMPNN, Caliby, RhoMax / OPTICS, QM/MM shortlist hooks.
