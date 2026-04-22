# RhGC Pocket Map Attempt — Negative Result

Status: **Option A from PR #10 attempted; approach failed; documenting for
future reference.**

This chunk was supposed to wire a structure-derived graded pocket map for
BeRhGC1 so the first design run's top-8 ties (PR #10) would break. It
didn't work, for instructive reasons.

## What was tried

**Step 1: Look for a real retinal-bound crystal structure of BeRhGC1.**
None exists. The three Butryn-lab PDBs (6AO9 / 6AOA / 6AOB, Butryn et al.
2017) are the *guanylyl cyclase domain only* — no retinal, no rhodopsin
domain. RCSB full-text search for "rhodopsin guanylyl cyclase" returns
those three plus unrelated GPCR entries (7VH0 etc.). PDB 6SIX that I
initially tried is an AMP ligase, unrelated.

**Step 2: Fall back to AlphaFold prediction for A0A060H1D7.**
AlphaFold DB has a v6 model (global pLDDT 76.56; K338 local pLDDT 83.12).
The model is available as
`https://alphafold.ebi.ac.uk/files/AF-A0A060H1D7-F1-model_v6.pdb`.

**Step 3: Since AF is apo (no retinal), synthesize retinal atoms.**
Two attempts:

- *First attempt:* extend 20 C atoms from K338 NZ along the CE→NZ vector.
  Output: 5 "pocket" residues, all in TM7 around K338 (pos 327–338). The
  retinal line went *out of* the protein core, not into it.
- *Second attempt:* extend 20 C atoms from K338 NZ along the
  NZ→(rhodopsin-domain CA centroid) vector, so retinal points into the
  protein. Output: 30 "pocket" residues spread across positions 133, 199–
  215, 327–341, 391–398 — but **none** of my heuristic mutable positions
  from PR #10 (D181, W183, Y272, D283) appeared.

## Why it fails

1. **AF-apo sidechain orientations around a ligand-binding pocket are
   unreliable.** Retinal binding stabilises a specific Lys Nε orientation
   and reorganises nearby residues. Without the ligand in the training data,
   AF's predicted pocket geometry is a "best guess" from sequence
   plausibility, not a physically refined cavity.
2. **A straight 20-Å line of carbons doesn't approximate retinal geometry.**
   Real retinal has a β-ionone ring, a polyene chain with specific
   cis/trans kinks, and a Schiff base end. Residues close to the ring end
   are different from residues close to the Schiff-base end. A line-based
   proxy can't capture this.
3. **Anchoring on NZ and extending along "toward the protein centroid"
   cuts through multiple helices.** It doesn't represent one helix bundle's
   interior; it's geometry applied to biology.

The 30-residue "pocket" from attempt 2 is biologically meaningless. Wiring
that as a PocketMap would feed graded scoring with noise.

## What would actually work

Four concrete paths, roughly in increasing order of effort:

1. **Proper sequence alignment (MAFFT / hmmer) of BeRhGC1 onto BR.** BR
   1C3W has a real retinal-bound pocket (we mapped it in PR #2). The
   alignment transfers BR's pocket residues to BeRhGC1 positions by
   homology. Output: a `PocketMap` whose distances are *inherited* from
   BR, not computed from a BeRhGC1 structure — but at least every residue
   in the map has a biologically defensible reason to be there.

2. **Dock retinal into the AF model with Chai-1 / AlphaFold3 / AutoDock.**
   Produces a ligand-bound BeRhGC1 model. Modern ML docking is reasonable
   for a known cofactor in a Rossmann-like fold. Output: real per-residue
   distances, computed the same way as for BR 1C3W.

3. **Run ColabFold with retinal as an MSA-biased ligand.** Similar to (2)
   but via ColabFold's hallucination-with-ligand path. Requires a GPU
   session or Google Colab.

4. **Wait for a BeRhGC1 crystal.** Ideal but slow; no current precedent.

Of these, **(1)** is lowest-effort and tool-light: a 30-minute alignment +
a short script that does the transfer. **(2)** is the right answer for a
thesis deliverable but needs either an Anthropic Chai-1 API or a local
install.

## Recommendation for the next attempt

Do **(1)**. Concretely:

- Install `mafft` (or use the EBI web service).
- Align BeRhGC1 A0A060H1D7 against BR P02945.
- Take the BR pocket from `configs/pocket_maps/BR_Hsal_BACR.json` (already
  in main from PR #5).
- For each pocket residue in BR, find its aligned position in BeRhGC1; if
  aligned to a gap, drop it.
- Write out `configs/pocket_maps/BeRhGC1_full_via_BR_alignment.json` with
  the same structure as the real pocket maps plus a clear
  `derivation_method: "alignment_transfer_from_BR_1C3W"` annotation.
- Rerun the RhGC first-run pipeline with `--graded-pocket`. Top-8 should
  reorder; the candidates at positions aligned to strong-band BR residues
  win, the others drop to medium or none.

The `cli pocket` parser from PR #7 can't consume this directly (it
computes distances from a PDB), so either:

- Add a `cli pocket-transfer` subcommand that takes two scaffold names +
  a source PocketMap + an alignment file and emits a target PocketMap.
- Or hand-write the JSON for this one-off demonstration.

The latter is faster for a single chunk; the former is the proper
engineering answer if homology transfer is going to be a recurring
pattern.

## Artifacts on this branch

Nothing committed under `configs/` or `opsin_pipeline/`. The AF model and
synthetic-retinal PDB stay in `data/structures/` (gitignored) for local
inspection only. The generated `configs/pocket_maps/BeRhGC1_full.json`
was **not** wired into `rhgc_design_scaffold.json` — graded scoring is
safer without it than with biologically noisy input.

This doc is the only commit. PR #10 (the first-run output) stays as the
landed artifact; this branch is the "we tried A, here's why it didn't
work" companion.
