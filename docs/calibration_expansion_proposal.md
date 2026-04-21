# Calibration Expansion Proposal

Status: **proposal (no data changes yet)** — landing target for "Path A" in PR #8.

Context: PR #8's real-ESM2 snapshot showed the 22-entry literature calibration
can't discriminate useful from disruptive on PLM delta alone. The data itself,
not the scorer, is the bottleneck — both classes cluster at retinal-pocket
positions with overlapping PLM plausibility distributions.

This doc proposes concrete additions to close that gap **before** Stage 4
scoring ships. If reviewers prefer Path B (ship Stage 4 off-by-default and
defer calibration expansion), this file stays as a reference for the follow-up
chunk; if Path A wins, the list below becomes the curation target.

## Target shape

Current set: 17 useful + 5 disruptive + 4 neutral, all at retinal-pocket
positions, dominated by BR + ChR2.

Target: **≥40 total, ≥10 disruptive, ≥6 neutral**, with entries chosen so at
least one of the following distinguishes useful from disruptive at each
position:
- **Charge reversal** (e.g. `D→R` vs `D→N` at position 85)
- **Size jump** (e.g. `A→W` vs `A→S` at a pocket-wall residue)
- **Non-pocket controls** (a few entries at surface-exposed residues that are
  biologically neutral — absent in the current set, which skews pocket-only)

The goal is entries whose class labels are determined by *substitution
chemistry at the same position*, so PLM delta (or any chemistry-aware scorer)
can discriminate them.

## Proposed additions (to review, verify, and reject as needed)

Notation: scaffold-numbered positions (same convention as
`calibration_literature.json`). Each entry needs PMID/DOI before landing.

### BR_Hsal_BACR (add ~10)

**Disruptive, charge-reversal at the counterion** (currently absent):
- [ ] `D85R` — positive residue at the primary counterion; expected catastrophic.
      Literature: inferred from structure (Braiman 1988 counterion model);
      may not have a direct mutagenesis paper — verify.
- [ ] `D85K` — same class as D85R.
- [ ] `R82D` — inverse charge reversal at the Arg of the counterion triad.

**Disruptive, pocket-lining hydrophobic → polar**:
- [ ] `L93N`, `W86Q` — hydrophobic pocket-wall residues replaced with polar.
      Literature on W86F spectral shifts exists (Lin & Mathies 1989); the
      `W86Q` chemistry-flip is an extrapolation.

**Useful, conservative at D85**:
- [ ] `D85E` — aspartate → glutamate, same charge, one carbon longer.
      Expected small spectral perturbation, retains counterion.
      Literature: Subramaniam et al. 1991 (Δλmax ~ +20 nm, retains pump).

**Neutral, surface-exposed control**:
- [ ] `K40R`, `K41R`, `R7K` — surface Lys/Arg swaps that don't touch the
      pocket or counterion triad. Expected biologically neutral; no published
      phenotype I know of (review needed).

### ChR2_Cre_full (add ~10)

**Disruptive, primary counterion charge-reversal** (currently absent; E123 is
the only ChR2 pocket glutamate in the set):
- [ ] `E123R`, `E123K` — analogous to D85R on BR.

**Disruptive, Schiff-base-adjacent**:
- [ ] `W223A`, `F217A` — aromatic pocket wall → alanine. Published as
      red-/blue-shifts with loss-of-function tendencies (Kato 2012,
      Schneider 2015 review).

**Useful, conservative DC-gate shifts**:
- [ ] `C128F`, `C128W` — aromatic substitutions at the DC gate. C128T/S/A are
      already in the set; adding aromatics at the same position lets PLM see
      chemistry variation *while fixing the position*.

**Neutral, cytoplasmic-domain controls**:
- [ ] `K325R`, `Q442N` — cytoplasmic-loop residues with no known phenotype.
      These break the pocket-only clustering; PLM delta should be near 0.

### Chrimson_5ZIH_structure (add ~5)

Currently no entries. Adding these requires first populating
`mutable_positions` on the scaffold. Suggested curated set:
- [ ] `K176R` — re-verify numbering (the original attempt failed; see
      ``rerun_report.md``). This is the Chrimson DC-gate analog of ChR2 C128.
- [ ] `A298T`, `S267M` — red-shift variants from the ChrimsonR paper
      (Klapoetke 2014) that are well-established.
- [ ] One Schiff-base-Lys disruptive mutation to parallel BR K229A.
- [ ] One cytoplasmic neutral control.

### GtACR1_6CSM_structure (add ~5)

Currently no entries. Most mutagenesis data for GtACRs is in Govorunova 2015,
Kim 2018, or the Sineshchekov lab publications. Suggested:
- [ ] Schiff-base K238 disruptive.
- [ ] One anion-selectivity residue (literature: E68A changes anion vs
      cation).
- [ ] Two pocket chemistry-flip pairs (size or polarity).

### CaRhGC_full (add ~5)

Currently no entries. Less literature than BR/ChR2 but the user's lab
background is strongest here. Suggested:
- [ ] `D122N` equivalent — if CaRhGC has a counterion analog (check
      alignment to BR D85).
- [ ] `K336R` or similar Schiff-base mutant.
- [ ] 2–3 residues in the guanylyl-cyclase domain as non-pocket controls.

## Ingest / landing plan

Once reviewers agree on entries:

1. Update `configs/calibration_scaffolds.json` to populate
   `mutable_positions` for Chrimson / GtACR1 / CaRhGC at the new positions
   (each allowed AA list includes the target residues of the new calibration
   entries).
2. Add the entries to `configs/calibration_literature.json` with populated
   `source` fields.
3. Re-run `cli plm` on each scaffold to refresh
   `configs/plm_predictions/*.json` with deltas covering the new targets.
4. Rerun `cli run` with and without `--graded-pocket` to capture the new
   AUROC baseline.
5. Write `rerun_report_expanded_calibration.md` with the new numbers.

Estimated effort: half a day of careful literature curation; a follow-up
commit of a few hundred JSON lines. No code changes — the ingest / generate
paths already handle the new data as-is.

## Acceptance criteria for Path A

Before Stage 4 scoring is unblocked:

- [ ] ≥40 total calibration entries.
- [ ] ≥10 disruptive entries (currently 5).
- [ ] ≥1 charge-reversal pair per scaffold (at least for BR + ChR2).
- [ ] ≥3 non-pocket "neutral" controls.
- [ ] `cli run` with the expanded set still matches 100 % of calibration
      entries (no stale IDs).
- [ ] At least one of {PLM delta alone, graded pocket alone, graded pocket ×
      PLM} moves AUROC **above 0.5** on the expanded set — otherwise the
      finding is that the current signals don't discriminate at all, and
      Stage 4 merges with the flag off-by-default.

## Non-goals

- This proposal does **not** mandate Path A. Path B (ship Stage 4
  off-by-default, defer calibration expansion) is still a legitimate choice
  and will leave this doc as a reference.
- No existing calibration entries are removed or relabelled.
- No Chrimson/GtACR1/CaRhGC scaffold sequences change; only their
  `mutable_positions` would be populated when entries are added.
- Literature curation is the bottleneck, not code.

## Ownership / review

If Path A is chosen:

- **Codex** is best placed to do the literature curation (familiarity with the
  original `calibration_literature.json` curation, standard references
  already at hand).
- **User's lab context** (Hegemann-lab RhGC spectral tuning) covers the
  CaRhGC entries that neither Codex nor I have strong priors on.
- **I** can re-run the ESM2 and calibration-rerun pipeline once the
  literature JSON lands, and do the AUROC writeup.

This doc is a claimable handoff surface, not a commitment.
