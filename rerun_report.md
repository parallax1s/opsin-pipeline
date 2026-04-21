# Calibration Rerun Report

Branch: `literature-calibration-set`  
Date: 2026-04-21

## Fixes Applied

- Shifted bacteriorhodopsin mature-protein calibration positions onto full UniProt `P02945` sequence coordinates:
  - `R82 -> full p95`
  - `D85 -> full p98`
  - `D96 -> full p109`
  - `K216 -> full p229`
- Verified scaffold wild-type identities for all referenced BR and ChR2 calibration positions.
- Removed the unverified `Chrimson K176R` seed from the calibration set pending full-reference numbering resolution.
- Added current-score-readable reason tags to mutable positions, including `retinal_pocket`, `spectral_tuning`, `photocurrent`, and related role terms.

## Command

```bash
PYTHONPATH=. python3 -m opsin_pipeline.cli run \
  --scaffolds configs/calibration_scaffolds.json \
  --calibration configs/calibration_literature.json \
  --target-family bacteriorhodopsin \
  --target-phenotype spectral_tuning \
  --max-mutations 2 \
  --max-combinations-per-scaffold 500 \
  --per-scaffold-cap 50 \
  --per-position-cap 20 \
  --calibration-top-k 20 \
  --out out/calibration_rerun
```

## Results

- Scaffolds screened: 5
- Candidates generated: 156
- Useful matched: 17 / 17
- Disruptive matched: 5 / 5
- Neutral matched: 4 / 4
- AUROC useful vs disruptive: 0.4588
- MRR of useful mutants: 0.1376
- Useful in top 20: 4
- Random baseline useful in top 20: approximately 2.179
- Disruptive in top 20: 2

## Interpretation

This is now a valid Stage 1 baseline rather than a numbering artifact. The below-random AUROC means the current reason-tag heuristic still cannot separate useful retinal-pocket mutations from disruptive retinal-pocket mutations. That is expected. Both classes often share the same coarse tags, so Stage 2 needs the structure-grounded pocket distance signal and later PLM or biophysical priors.

The calibration set is still a starter set. Chrimson and RhGC entries need expansion after numbering/domain review.
