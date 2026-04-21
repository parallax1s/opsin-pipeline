# Real ESM2 PLM Predictions — Calibration Snapshot

Branch: `stage-3-real-predictions`
Date: 2026-04-21
Model: `esm2_t12_35M_UR50D` (facebook, 35M params, CPU forward pass)

## Setup

Dependencies installed in an isolated venv at `/tmp/esm-venv` — not committed,
not required to run the default 128-test suite. Generated real PLM predictions
for all 5 calibration scaffolds by running `cli plm` (no `--mock`):

| Scaffold | Preds | Median delta |
|---|---:|---:|
| BR_Hsal_BACR | 10 | −2.500 |
| ChR2_Cre_full | 15 | −2.953 |
| Chrimson_5ZIH_structure | 0 | — |
| GtACR1_6CSM_structure | 0 | — |
| CaRhGC_full | 0 | — |

Three scaffolds produce zero predictions because they don't yet carry
`mutable_positions` in the calibration set (noted follow-up from the
earlier `rerun_report.md`).

## Plumbing unchanged — calibration AUROC stays 0.4588

Rerunning the full pipeline with real predictions wired in:

```
calibration: useful 17/17, disruptive 5/5, AUROC 0.4588
```

Identical to the mock-data rerun and to the pre-PLM legacy baseline, as
required by spec §7.5. The real deltas reach the `min_plm_log_likelihood_delta`
CSV column but do not influence scoring — Stage 3's "no-scoring-change"
contract holds.

## Raw per-entry PLM deltas

```
SCAFFOLD           MUT        CATEGORY     PLM_DELTA  LABEL
BR_Hsal_BACR       p98N       useful         -3.695  D85N
BR_Hsal_BACR       p98T       useful         -2.701  D85T
BR_Hsal_BACR       p98S       useful         -2.724  D85S
BR_Hsal_BACR       p98V       useful         -2.299  D85V
BR_Hsal_BACR       p109N      neutral        -1.135  D96N
BR_Hsal_BACR       p95Q       neutral        -3.200  R82Q
BR_Hsal_BACR       p229A      disruptive     -0.169  K216A
BR_Hsal_BACR       p229G      disruptive     -1.046  K216G
ChR2_Cre_full      p134R      useful         -2.953  H134R
ChR2_Cre_full      p123T      useful         -2.961  E123T
ChR2_Cre_full      p123A      useful         -2.709  E123A
ChR2_Cre_full      p128T      useful         -0.967  C128T
ChR2_Cre_full      p128S      useful         -1.456  C128S
ChR2_Cre_full      p128A      useful         -0.612  C128A
ChR2_Cre_full      p156A      useful         -1.517  D156A
ChR2_Cre_full      p132C      useful         -5.471  L132C
ChR2_Cre_full      p159C      useful         -4.097  T159C
ChR2_Cre_full      p134K      neutral        -3.411  H134K
ChR2_Cre_full      p159S      neutral        -3.029  T159S
ChR2_Cre_full      p120A      disruptive     -3.090  R120A
ChR2_Cre_full      p257A      disruptive     -1.599  K257A
ChR2_Cre_full      p90A       disruptive     -1.682  E90A
```

## The negative result worth reading carefully

**PLM delta alone does NOT separate useful from disruptive. On BR, it
anti-correlates.**

Mean PLM delta by category:

| Class | BR | ChR2 | All |
|---|---:|---:|---:|
| useful (n=13) | −2.86 | −2.60 | −2.68 |
| neutral (n=4) | −2.17 | −3.22 | −2.69 |
| disruptive (n=5) | **−0.61** | **−2.12** | **−1.52** |

- **Disruptive deltas are LARGER (less negative)** than useful on BR, and only
  slightly more negative than useful on ChR2.
- Treating "higher delta = more plausible = more useful", AUROC of useful vs
  disruptive would be *below* 0.5 — probably around 0.2–0.3 on BR, around 0.45
  on ChR2.

## Why this happens (mechanistic diagnosis)

ESM2 learns plausibility from natural-sequence co-occurrence. For opsin
calibration positions:

- **D85 in BR** is deeply conserved across the opsin family (counterion
  triad). ESM2 penalises *any* mutation there, including the useful
  spectrum-shifting `D85N/T/S/V` variants. Big negative deltas for useful
  mutations.
- **K216 in BR**, while *functionally* essential (covalent Schiff base),
  is not as strictly *sequence*-conserved in the training-data sense —
  different opsins anchor retinal at different positions or via different
  chemistry, and many related proteins have non-K residues there. ESM2
  doesn't strongly penalise `K216A/G`. Small deltas for disruptive mutations.
- The function that distinguishes useful from disruptive at the *same*
  position (e.g. `D85N` keeps a polar H-bond partner, `D85A` removes the
  counterion entirely) is chemistry — and PLM delta doesn't see chemistry,
  only sequence plausibility.

**This is interesting information, not a broken model.** ESM2 captures
"which residues are sequence-conserved across homologs," not "which
substitutions preserve function." For the Stage 4 scoring proposal, this
changes the hypothesis.

## Implications for Stage 4 scoring

The spec (`docs/plm_adapter_spec.md` §11) anticipated this possibility:

> If "PLM delta alone" already moves AUROC above 0.5 on our calibration,
> the combined signal might not help — we still ship Stage 3 + Stage 4 but
> document the finding. If neither alone nor combined moves AUROC above
> 0.5, Stage 4's merge waits on better calibration (more entries,
> chemically-diverse controls) rather than flipping a default on bad
> evidence.

We are now in the "neither alone moves AUROC above 0.5" regime on this
22-entry calibration set. Two realistic paths:

**Path A — bigger calibration set first.** Expand the literature set with
more chemically-diverse single-mutant entries (conservative substitutions,
small hydrophobic swaps, charge reversals, etc.) so the useful/disruptive
populations are better separated chemically. Then re-evaluate. This is a
curation task, not a code task.

**Path B — ship Stage 4 scoring with the flag off, by convention.** Stage 4
lands a combined scorer (e.g. `alpha * pocket_band + beta * plm_delta`) but
remains opt-in behind `--graded-pocket --use-plm-chemistry`. Nobody relies
on it until someone runs the calibration-expanded evaluation and flips the
default. Low-risk: keeps the architecture in place, defers the evidence
question.

Both paths make Stage 4 a more reflective chunk than "ship a formula that
moves AUROC." The spec's original framing ("Stage 4 merges only if AUROC
moves above 0.5") rules out Path B as written. This report proposes a
spec amendment in the next Stage 4 spec iteration: Stage 4 may ship in
`off-by-default` mode as a documented negative-evidence landing, with the
default-flip contract reassigned to the calibration-expansion chunk.

## What this PR ships

- `configs/plm_predictions/*.json` — real ESM2 predictions replacing the
  mock ones. Same file paths; wired the same way.
- `rerun_report_plm.md` — this document.
- No code changes. All tests (128/128 with 1 opt-in skipped) continue to
  pass because the scorer interface is unchanged and the CSV shape is
  identical to the mock-data case.

## Reproducibility

```bash
python3 -m venv /tmp/esm-venv
/tmp/esm-venv/bin/pip install transformers torch
for s in BR_Hsal_BACR ChR2_Cre_full Chrimson_5ZIH_structure GtACR1_6CSM_structure CaRhGC_full; do
  PYTHONPATH=. /tmp/esm-venv/bin/python -m opsin_pipeline.cli plm \
    --scaffold "$s" \
    --scaffolds configs/calibration_scaffolds.json \
    --model esm2_t12_35M_UR50D \
    --out "configs/plm_predictions/$s.json"
done
```

Each PocketMap PLM prediction carries `sequence_sha256` so divergence from
a future rerun is visible. `model_id` recorded as `esm2_t12_35M_UR50D` on
every entry.
