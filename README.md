# Opsin Pipeline MVP

Local dependency-free MVP for broad in silico opsin candidate triage.

The package does not run LigandMPNN, Caliby, RhoMax, METL, or QM/MM directly yet. It defines the data formats, candidate generation, transparent scoring, ranking, and report outputs that those adapters can plug into later.

## Run Tests

```bash
python3 -m unittest discover -s opsin_pipeline/tests -v
```

## Run Example

```bash
PYTHONPATH=opsin_pipeline python3 -m opsin_pipeline.cli \
  --scaffolds opsin_pipeline/configs/example_scaffolds.json \
  --out opsin_pipeline/out/example \
  --target-family RhGC \
  --target-phenotype spectral_tuning
```

Outputs:

- `ranked_candidates.csv`
- `decision_report.md`

## Real Scaffold Seed Set

`configs/real_scaffolds_seed.json` and `configs/real_scaffolds_seed.fasta` contain initial scaffold sequences pulled from UniProt and RCSB:

- ChR1 and ChR2 full UniProt sequences
- CaRhGC and BeRhGC1 UniProt sequences
- BeCNG1 as a cGMP-gated channel partner candidate
- PDB-backed ChR2, C1C2, Chrimson, GtACR1, and ChRmine structure sequences

The real seed config intentionally leaves `mutable_positions` empty. Add protected and mutable positions only after sequence/structure numbering is reviewed.

## Current Scoring Signals

- target family match
- target phenotype match
- retinal-pocket mutation tag
- assay-readiness tag (`growth_selection` or `biochemical`)
- small mutation-count penalty
- protected-residue violation penalty

This score is a triage heuristic. Treat it as a reproducible queue for deeper structure, physics, and wet-lab validation.
