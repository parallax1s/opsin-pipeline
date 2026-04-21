# Graded Pocket Scoring — Calibration Comparison

Branch: `pocket-maps-for-calibration`
Date: 2026-04-21

## Setup

Produced scaffold-numbered `PocketMap` JSONs for 4 of the 5 calibration
scaffolds by running `cli pocket` against real RCSB PDB structures with
`--pdb-offset` chosen so PDB author numbering aligns onto scaffold sequence
coordinates (AA verification clean, no `mapping_note`s):

| Scaffold | PDB | Offset | Ligand | Pocket residues | Schiff-base K |
|---|---|---:|---|---:|---|
| `BR_Hsal_BACR` | 1C3W | +13 | `RET` | 29 (17 strong, 7 medium, 5 none) | K229 @ 1.34 Å |
| `ChR2_Cre_full` | 6EID | 0 | `LYR` | 28 (14 strong, 9 medium, 5 none) | K257 (Schiff-base, 0 Å) |
| `Chrimson_5ZIH_structure` | 5ZIH | −3 | `LYR` | 28 (16 strong, 7 medium, 5 none) | K296 (Schiff-base, 0 Å) |
| `GtACR1_6CSM_structure` | 6CSM | −3 | `RET` | 25 (14 strong, 6 medium, 5 none) | K235 @ 1.47 Å |

`CaRhGC_full` has no published structure and stays unwired — candidates under
that scaffold fall back to the reason-string rule under graded mode per
spec §7.3.

## Command

```bash
# Legacy (flag off)
cli run --scaffolds configs/calibration_scaffolds.json \
        --calibration configs/calibration_literature.json \
        --target-family bacteriorhodopsin \
        --target-phenotype spectral_tuning \
        --max-mutations 2 --max-combinations-per-scaffold 500 \
        --per-scaffold-cap 50 --per-position-cap 20 \
        --calibration-top-k 20 --out out/legacy

# Graded (flag on, pocket maps present on 4/5 scaffolds)
cli run --scaffolds ... --graded-pocket --out out/graded
```

## Results

| Mode | Useful matched | Disruptive matched | AUROC (useful vs disruptive) | MRR useful | Useful in top 20 |
|---|---:|---:|---:|---:|---:|
| Legacy | 17/17 | 5/5 | **0.4588** | 0.1376 | 4 |
| Graded | 17/17 | 5/5 | **0.4176** | 0.1241 | 3 |

Random AUROC baseline = 0.5; random useful-in-top-20 ≈ 2.18.

## Interpretation — graded is not a win here

Graded scoring is **worse** than legacy on this calibration. The reason is
diagnostic, not a bug:

1. **Useful and disruptive mutations share the pocket.** Every calibration
   entry we curated is at a retinal-pocket residue (that's how they ended up
   worth publishing about). So both classes land in the `strong` band and get
   the same `+2.0` pocket bonus from graded scoring. Same signal on both sides
   → no discrimination.
2. **Mutation-count penalty asymmetry.** 4 of 17 usefuls are doubles (each
   carries `-0.4` mutation penalty); all 5 disruptives are singles (each
   `-0.2`). Graded scoring amplifies pocket signal equally across singles
   and doubles, so the total scores land in roughly the same place and the
   penalty gap pushes usefuls slightly below disruptives — dragging AUROC
   below the legacy baseline.
3. **Structure geometry is position, not chemistry.** `D85N` and `D85A` sit at
   the same position, same distance to retinal, same graded band. The signal
   that separates them is the *identity of the substitution* (N preserves a
   polar contact, A removes the counterion), which graded scoring doesn't
   consider. This is exactly the blind spot a PLM plausibility prior or a
   LigandMPNN substitution score is designed to fill.

## Decision

**`--graded-pocket` default stays OFF.** Flipping it on against the current
calibration would trade one below-random AUROC (0.46) for a slightly-more-
below-random AUROC (0.42) without improving ranking quality. The pocket-map
plumbing is still valuable:

- It surfaces geometric evidence into the review workflow
  (`pocket-annotate`), which is independently useful.
- It's a necessary building block for Stage 3 (PLM plausibility adapter).
  A PLM score conditioned on *the substitution* combined with graded
  *position* score should discriminate D85N from D85A — neither alone will.
- The flag remains available for experimentation.

## Next concrete chunk

**Stage 3 — PLM plausibility adapter.** Per spec roadmap:

- Thin adapter that scores each candidate's mutated sequence under a
  small protein language model (ESM2-150M, or SaProt, or ESM-C distilled).
- Emits per-mutation `plm_log_likelihood_delta` on the `Mutation` object.
- New score key `plm_plausibility` with graded bands.
- Combined scoring: `position_score * chemistry_score` where `position_score`
  is the existing graded band and `chemistry_score` is the PLM-derived
  plausibility. Re-evaluate AUROC on the same literature calibration —
  hypothesis is that the *combination* moves AUROC above 0.5 even though
  neither component does alone.
- Keep the PLM import optional so the pipeline still runs without torch.

Estimate: one chunk of similar size to Stage 2 (plumbing + one test fixture
using a tiny pretrained model or mocked log-likelihoods).
