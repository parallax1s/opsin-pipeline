# RhGC First Design Run — Notes

Status: **exercise, not final candidate selection.** The purpose of this chunk
was Claude's Option 1: use the pipeline on a real target to find out what
works, what breaks, and what's missing — not to ship wet-lab-ready variants.
Read this doc before trusting anything in ``out/rhgc_first_run/``.

## Scaffold correction surfaced during the run

The calibration set's ``CaRhGC_full`` (UniProt A0A1Y2HLQ0, 541 aa) is NOT a
rhodopsin-guanylyl cyclase — it's a cyclase-domain-only paralog from
*Catenaria anguillulae*. Evidence:

- No Schiff-base lysine in the expected TM7 region (300–350).
- Full classic class III cyclase catalytic core (``IGDAYL``) at position ~409.
- Residues 1–300 contain classic microbial-opsin *motifs* but at shifted
  positions (e.g. ``DYAL`` conserved at 258 instead of ~280), suggesting a
  truncated/mis-annotated variant.

Searching UniProt for the actual rhodopsin+cyclase fusion returns
**A0A060H1D7 (BeRhGC1, 626 aa)** from *Blastocladiella emersonii* — this is
the widely-studied RhGC from Scheib 2018 / Trieu 2017 / Butryn 2020. Switched
to this scaffold. The ``CaRhGC_full`` entry should probably be renamed
``A0A1Y2HLQ0_cyclase_paralog`` in the calibration set, or dropped — it's
misleading as "RhGC" when it can't bind retinal.

## Scaffold shipped

``configs/rhgc_design_scaffold.json`` — BeRhGC1_full, 626 aa, UniProt
A0A060H1D7. Key annotations:

- **Protected**: K338 (retinal Schiff-base Lys, confirmed by sequence
  context ``NVEKSGLRW`` in TM7) + positions 499–530 (class III cyclase
  catalytic core, starting at the ``IGDAYL`` motif). 33 protected residues
  total.
- **Mutable**: 4 positions, all flagged ``HEURISTIC`` in the
  ``configs/rhgc_position_review.csv``:
  - `D181` → `E, N` (TM3 counterion candidate, BR D85 analog)
  - `W183` → `F, Y` (TM3 aromatic wall, BR W86 analog)
  - `Y272` → `F, W` (TM5 aromatic, speculative tuning position)
  - `D283` → `E, N` (TM6 NDYAL-motif Asp, BR D96 analog)

**Position selection is sequence-homology-guess only. No alignment tool, no
structure, no Hegemann-lab curation.** Before running anything wet-lab from
this, an opsin-expert pass over the mutable set is required.

## Run command

```bash
PYTHONPATH=. python3 -m opsin_pipeline.cli run \
  --scaffolds configs/rhgc_design_scaffold.json \
  --out out/rhgc_first_run \
  --target-family RhGC --target-phenotype spectral_tuning \
  --max-mutations 2 --per-scaffold-cap 20 --per-position-cap 5
```

`--graded-pocket` was omitted because no RhGC crystal structure is wired into
the scaffold. (BeGC1 has PDB 6SIX from Butryn 2020 — adding it would be a
natural follow-up chunk.)

Real ESM2 predictions were generated for BeRhGC1_full via
``cli plm --scaffold BeRhGC1_full --model esm2_t12_35M_UR50D`` and wired in
through `plm_predictions_path`. 8 per-mutation deltas attached to the
candidates.

## Results

32 candidates generated (8 singles + 24 doubles). Top 10 all tie at total
score 7.8 (singles) or 7.6 (doubles) — **the legacy score function cannot
discriminate them.**

| Rank | Candidate | PLM delta | Comment |
|---:|---|---:|---|
| 1 | BeRhGC1_p283DtoN | −1.87 | Highest PLM plausibility; ESM2 thinks D283 is not strictly conserved |
| 2 | BeRhGC1_p283DtoE | −2.85 | Conservative D→E |
| 3 | BeRhGC1_p272YtoW | −4.44 | Most penalised by PLM — large aromatic swap |
| 4 | BeRhGC1_p272YtoF | −3.14 | Moderate |
| 5 | BeRhGC1_p183WtoY | −2.57 | |
| 6 | BeRhGC1_p183WtoF | −2.16 | |
| 7 | BeRhGC1_p181DtoN | −3.76 | D181 strongly conserved by PLM — supports counterion hypothesis |
| 8 | BeRhGC1_p181DtoE | −3.74 | Same strong conservation |

All 8 get the same legacy total, so the ordering within top-8 is pure
alphabetical tie-break (candidate_id desc). **The ranking is currently not
actionable.**

## What the run actually told us

1. **Scaffold data quality is the first real bottleneck.** The existing
   calibration "CaRhGC" scaffold can't be used for real design; swapping to
   BeRhGC1 took 30 min of UniProt searching + sequence-motif scanning that a
   human-curated scaffold registry should have done upfront.
2. **Legacy scoring is useless at this stage.** 8 candidates, 1 score value.
   Without either a graded pocket map (needs a PDB / AlphaFold) or a scoring
   signal that weights PLM delta, the tool produces candidate *lists* but not
   candidate *rankings*.
3. **Real PLM data is informative as descriptive evidence, not ranking.** The
   deltas above tell a small biological story: `D181` is highly conserved
   (likely real counterion), `D283` is more permissive (possibly secondary).
   A biophysicist could read this column and narrow down by hand — exactly
   the kind of "tool for a scientist, not a replacement for one" use case
   Claude's Option 1 was meant to surface.
4. **Graded pocket scoring would matter here.** With real pocket distances
   from a BeGC1 crystal (PDB 6SIX) or AlphaFold model, positions split into
   strong/medium/none and the top-4 would diverge instead of tying. That's a
   concrete follow-up chunk.

## Failure modes to fix next

In descending priority:

1. **No BeRhGC1 structure wired.** Adding PDB 6SIX (Butryn 2020) or an
   AlphaFold model would unblock graded-pocket scoring on this scaffold,
   which is the only thing that would actually reorder the top candidates.
2. **Mutable positions are sequence-heuristic.** An alignment against BR or
   ChR2 (via MAFFT or hmmer) would put proper homologous landmarks at
   correct positions; right now `Y272` and `D283` are especially soft
   choices.
3. **Scaffold registry has a mis-labelled entry.** The calibration's
   "CaRhGC_full" is a cyclase paralog, not an RhGC. Either rename it or
   drop it.
4. **No target for the `spectral_tuning` objective.** We asked for
   "spectral tuning" but didn't say red or blue. In a real design run this
   would be parametrised (e.g. target Δλ > +20 nm → favour certain
   substitutions).
5. **PLM anti-correlation from PR #8 also applies here.** ESM2 penalises
   conserved positions regardless of whether the substitution is functionally
   useful or disruptive. For the design run, PLM is informative for
   *describing* candidates, not for re-ranking them.

## Sanity checks I did by hand

- K338 really is a K in the sequence: ✓ (sequence[337] == 'K').
- `IGDAYL` really is at position 499: ✓.
- D181 really is at position 181: ✓ (`DIWYGY` motif).
- `NDYAL` motif really is around position 283: ✓ (`NDYAL` at 282–286).
- Cyclase catalytic range 499–530 covers ~30 residues including metal-
  coordinating positions by homology to class III cyclases.

## Deliverables landed on this branch

- `configs/rhgc_design_scaffold.json` — BeRhGC1_full scaffold with reviewed
  heuristic mutables + protected core.
- `configs/rhgc_position_review.csv` — per-position rationale + confidence.
- `configs/plm_predictions/BeRhGC1_full.json` — real ESM2 deltas.
- `out/rhgc_first_run/ranked_candidates.csv` — 32 candidates.
- `out/rhgc_first_run/decision_report.md` — decision report (legacy scoring).
- `docs/rhgc_first_run_notes.md` — this file.

## Next concrete chunk (suggested)

Pick ONE of:

**Option A — wire BeGC1 structure + rerun graded.** Fetch PDB 6SIX, run
`cli pocket --pdb 6SIX.pdb --pdb-offset <tbd>` for BeRhGC1, attach
`pocket_map_path`, rerun with `--graded-pocket`. The top-8 stop tying. This
is the highest-leverage next chunk and probably a single cycle of work.

**Option B — align BeRhGC1 to BR and re-pick mutable positions.** Replace
the 4 heuristic positions with alignment-derived homologs of well-studied BR
residues (D85 → BeRhGC1-equivalent; T89, Y185, W86, W182). Needs MAFFT or
hmmer. Produces a scaffold worth running wet-lab candidates from.

**Option C — both A + B.** Full real design iteration. Probably a weekend.

My recommendation: **Option A first** because it's cheap and immediately
tests whether graded-pocket scoring reorders the top-8 meaningfully. Option
B unlocks scientific validity; Option A unlocks pipeline usefulness.
