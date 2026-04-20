# Structure-Grounded Pocket Parser — Spec

Status: **draft (pre-implementation)** — lock before coding.

## 1. Motivation

The current score function treats "retinal pocket" as a binary tag lifted from whatever string the scaffold author typed in `reason`. Two problems fall out of that:

1. Pocket membership is annotation-driven, not geometry-driven. A scaffold author could forget the tag on a real pocket residue, or over-apply it, and neither would be caught.
2. The score contribution is the same `+2.0` whether a residue is in direct van-der-Waals contact with retinal or five shells away. The calibration smoke run on `main` already demonstrates the consequence: AUROC = 0.1667 against the synthetic calibration, because `K8F` (useful) and `K8W` (disruptive) both sit at the same pocket position and score identically.

This chunk replaces the annotation with geometric evidence:

- Parse a PDB structure for the scaffold.
- Identify the retinal-like ligand(s).
- Compute each residue's minimum distance to any heavy atom of the ligand.
- Surface those distances into the review workflow and the score function.

## 2. Non-goals (explicitly out of scope for this chunk)

- mmCIF / PDBx parsing. RCSB defaults to mmCIF now; we defer it to a follow-up. Every opsin structure we need is also available in legacy PDB format.
- Automatic structure fetch from RCSB. The user places `.pdb` files into a known directory. No network calls at pipeline runtime.
- Hydrogen-aware distances. PDB rarely ships with Hs; we compute over heavy atoms only.
- Multi-model NMR ensembles. Take model 1.
- Full binding-pocket volumetrics (Voronoi, CASTp-style), pocket druggability, cavity detection. We only measure distances.
- Multi-chain oligomer analysis. We operate on one scaffold = one chain. Symmetry mates and crystal packing contacts are ignored.
- Alignment of the PDB sequence against the scaffold sequence. Numbering reconciliation is an explicit, user-provided offset for this chunk (see §6); pairwise-alignment fallback is a future extension.
- Any automatic setting of `mutable=true` or `protected=true`. The parser emits evidence only; the existing `position_map` human-review layer retains sole authority over what becomes mutable.

## 3. Module layout

```
opsin_pipeline/
  structure/
    __init__.py
    pdb.py            # pure-Python PDB ATOM/HETATM parser (no Biopython)
    ligands.py        # retinal-like HETATM whitelist + LYR handling
    pocket.py         # distance calculation, PocketResidue data model
  position_map.py     # existing; gains optional pocket-annotation step
  score.py            # gains graded pocket-distance signal
  schemas.py          # MutablePosition gains optional distance_to_retinal field

tests/
  fixtures/
    structures/
      synthetic_mini.pdb   # hand-authored 8-residue fragment + RET ligand
      br_1c3w_strip.pdb    # stripped real bacteriorhodopsin (retinal + ~20 nearest residues)
  test_pdb_parser.py
  test_pocket_distance.py
  test_ligands.py
```

No new third-party dependencies. `math.sqrt` and tuple arithmetic only.

## 4. Parser output data model

Two layered types.

### 4.1 `PDBAtom` (internal)

```python
@dataclass(frozen=True)
class PDBAtom:
    record: str                    # "ATOM" | "HETATM"
    serial: int
    atom_name: str                 # " CA ", " N1 ", etc. (stripped)
    alt_loc: str                   # usually ""; kept for alt-conf filtering
    res_name: str                  # "LYS", "RET", "LYR", ...
    chain: str                     # single character, "A", "B", ...
    res_num: int                   # PDB author resnum
    ins_code: str                  # insertion code (usually "")
    x: float
    y: float
    z: float
    occupancy: float
    element: str                   # "C", "N", "O", "H", ...
```

### 4.2 `PocketResidue` (public)

```python
@dataclass(frozen=True)
class PocketResidue:
    pdb_chain: str
    pdb_resnum: int
    pdb_ins_code: str              # "" unless the PDB uses insertion codes
    res_name_three: str            # "LYS"
    res_name_one: str              # "K"
    min_distance_A: float          # to nearest heavy atom of any matched ligand
    closest_retinal_atom: str      # atom name, for debugging
    closest_ligand_id: str         # "RET" | "LYR" | ...
    closest_ligand_resnum: int
    band: str                      # "strong" | "medium" | "none"  (per §7 thresholds)
```

### 4.3 `PocketMap` (scaffold-level public output)

```python
@dataclass(frozen=True)
class PocketMap:
    scaffold_name: str
    pdb_path: str
    pdb_sha256: str                # reproducibility: which file, exactly
    pdb_chain: str                 # the chain we analyzed
    ligand_matches: list[LigandMatch]
    pocket_residues: list[PocketResidue]   # sorted ascending by (chain, resnum, ins_code)
    cutoff_A: float                # the distance cap used (default 6.0; see §7)
    thresholds_A: tuple[float, float]   # (strong_max, medium_max) = (4.0, 5.5)
```

`LigandMatch` records which HETATM group was identified and why:

```python
@dataclass(frozen=True)
class LigandMatch:
    res_name: str                  # "RET" | "LYR" | ...
    chain: str
    res_num: int
    ins_code: str
    heavy_atom_count: int
    matched_by: str                # "whitelist" | "lyr_schiff_base" | "heavy_atom_heuristic"
```

`PocketMap` is what the review layer consumes. It is **not** the scaffold JSON. It is serialized to JSON alongside the scaffold for provenance:

```
out/pocket/<scaffold>.json
```

## 5. Ligand identification

### 5.1 Whitelist (first pass)

HETATM `res_name` matches against a curated set. Initial list:

| Code | Name | Notes |
|---|---|---|
| `RET` | all-trans retinal (free) | most common in crystal opsins |
| `RTN` | retinal (alt code) | seen in some older entries |
| `LYR` | retinal Schiff base to Lys | **covalent**; see §5.3 |
| `A1H` | all-trans retinal (variant) | e.g. ChRmine |
| `A1C` | 13-cis retinal | |
| `13M` | 13-cis retinal (alt) | |
| `9CR` | 9-cis retinal | rhodopsin photocycle intermediates |
| `BCR` | β-carotene | exclude — not retinal; listed for deny-list clarity |
| `RET_A`, `RT0`, `RET2` | vendor/refinement variants | accept; regex `^RE[T0-9A-Z]$` as safety net |

The whitelist is stored as a module constant, imported from `opsin_pipeline/structure/ligands.py`, and easy to extend. No silent expansion via heuristics — unknown HETATM codes are logged and ignored, not assumed to be retinal.

### 5.2 Deny-list (important)

Rhodopsin structures routinely carry other heterogeneous residues that must **not** be treated as the pocket ligand: lipid tails (`PLM`, `CLR`, `PC1`, `POPC`, `OLC`), detergents (`OGA`, `LMT`, `LDA`, `LUT`), fusion-protein ligands (from T4L / BRIL cocrystals), crystallization buffers (`ACT`, `GOL`, `PEG`, `TRS`). The parser hard-ignores any HETATM not on the allow-list.

### 5.3 `LYR` special case

Bacteriorhodopsin and several opsins in the PDB encode the retinal-bound Lys as a single covalent `LYR` residue (Lys + retinal connected through the Schiff base). Treating `LYR` naïvely breaks the parser two ways:

1. It's a HETATM on some entries, an ATOM-tagged residue on others. Must match either record type.
2. Its "protein" atoms (CA, CB, CG, CD, CE, NZ) belong to the Lys half; the Schiff-base N and the retinal polyene belong to the ligand half. For distance calculations, only the retinal-half atoms should count as "ligand heavy atoms".

Resolution: when we detect `LYR`, we split it. All atoms whose names match the Lys set (`N`, `CA`, `C`, `O`, `CB`, `CG`, `CD`, `CE`, `NZ`) are treated as the protein-side retinal-binding lysine (and its position becomes a pocket residue at distance 0). All other atoms (the retinal carbons `C1`–`C20`, `C1'`, etc.) are treated as ligand heavy atoms. The retinal-binding Lys is always output as a pocket residue with band `strong`, regardless of cutoff.

### 5.4 Failure behavior

- Zero matched ligands across the whole file → **hard error**: `NoRetinalLigandError`, with the list of HETATM codes actually observed so the user can update the whitelist if needed.
- Multiple matched ligands (common in oligomer asymmetric units) → **soft warning**. The parser picks a ligand deterministically: (a) if the caller passed `--pdb-chain X`, use the ligand on chain X; (b) else use the first ligand by (chain, resnum) sort order, and emit `LigandMatch` entries for all so the user can see what else was around.
- Heavy-atom count out of range for retinal (retinal is C₂₀; allow 15–25 heavy atoms as a sanity check) → **hard error** unless `--allow-weird-ligand` is passed. Catches mis-whitelisted entries.

## 6. PDB numbering ↔ scaffold sequence numbering

This is the single most error-prone surface and the one the spec must nail down.

### 6.1 Facts of life

- PDB author numbering is whatever the crystallographer typed. For opsins it routinely differs from UniProt because:
  - Signal peptide cleaved (ChR2: UniProt 1–20 absent in mature protein; PDB resnum 1 = UniProt 21, typically).
  - N- or C-terminal fusion tags (BRIL, T4L) inserted in the middle of the sequence, with their own resnum block.
  - Missing-density loops renumber discontinuously.
  - Truncations (Chrimson, GtACR1 N/C termini trimmed).
- Insertion codes (`82A`, `82B`) exist and must survive round-trips.
- PDB sequence (the amino acid actually at that resnum) may differ from our scaffold sequence at that position (engineered mutations in the crystal construct, e.g. `C128T` was crystallized in one ChR2 entry).

### 6.2 Rules

1. **The parser emits PDB author numbering only.** Chain, resnum, insertion code. It never silently produces a `seq_index` for the scaffold.
2. **Reconciliation is a separate, explicit step.** The caller provides one of:
   - `--pdb-offset N`: scalar offset such that `seq_index = pdb_resnum + N` for the selected chain. No insertion codes allowed when using an offset; the presence of an insertion code on any pocket residue becomes a hard error in this mode.
   - `--pdb-mapping path/to/mapping.csv`: explicit table with columns `pdb_chain, pdb_resnum, pdb_ins_code, seq_index`. One row per residue the caller wants mapped. Residues not in the mapping are emitted in PDB numbering and marked `review_needed=true`.
3. **Verification.** Whichever method produced `seq_index`, the parser must check `scaffold.residue_at(seq_index) == pdb_amino_acid`. Any disagreement (e.g., crystal construct has `C128T` but scaffold has `C128`) is flagged: the pocket residue is still emitted, but `review_needed=true` and the mismatch recorded in a `mapping_note` field. Default behavior: **warn and continue**; with `--strict-mapping` the mismatch becomes a hard error.
4. **No fallback to pairwise alignment in this chunk.** If the caller provides neither an offset nor a mapping, the parser emits PDB-numbered `PocketMap` only. The review workflow (§8) can still consume it — it just stays in PDB numbering.

### 6.3 Emitting into the position-map CSV

When a mapping is provided, pocket-annotation writes these fields into the existing `draft-position-map` CSV:

| Column | Source |
|---|---|
| `pdb_chain` | parser |
| `pdb_resnum` | parser |
| `region` | `retinal_pocket` appended if `band != "none"` |
| `notes` | `distance_to_retinal=4.2A; ligand=RET; mapping_note=...` |
| (new) `distance_to_retinal_A` | parser |
| (new) `pocket_band` | `strong` \| `medium` \| `none` |
| (new) `review_needed` | `true` if `mapping_note` was non-empty, else `false` |

`mutable` and `protected` are **never** written by the parser — they remain human-controlled in `apply-position-map`.

## 7. Graded distance scoring

### 7.1 Thresholds

| Band | Range | Default score contribution |
|---|---|---|
| `strong` | ≤ 4.0 Å | `+2.0` |
| `medium` | (4.0 Å, 5.5 Å] | `+1.0` |
| `none` | > 5.5 Å | `0.0` |

Thresholds are module-level constants and overridable via `score_candidate(..., pocket_thresholds=(4.0, 5.5))`. The 4.0 / 5.5 split comes from standard "first shell" (van der Waals contact, typically ≤4 Å) and "second shell" (polar / ordered-water bridge, ~4–6 Å) conventions in the opsin literature.

### 7.2 Candidate-level score

Let `positions` = the positions a candidate mutates. For each position, look up the pocket band (from the scaffold's pocket map, merged at load time — see §8). Candidate score:

```
pocket_score = max(band_points(position) for position in positions_with_pocket_data)
```

Rationale: a single strong pocket contact is the relevant signal. Averaging across mutations would dilute a strong hit; summing would over-reward shotgun combinations that happen to graze the pocket at multiple positions.

### 7.3 Backward-compat fallback

If a candidate has **no** mutation with `distance_to_retinal` data (scaffold lacks a pocket map), fall back to the existing reason-string rule: `+2.0` if any mutation's reason contains `retinal_pocket`. This keeps the existing `configs/example_scaffolds.json` and any user scaffolds without crystal structures working unchanged.

Scaffolds that **do** have a pocket map always prefer the graded signal; the reason-string fallback is ignored to avoid double-counting.

### 7.4 What changes in `score_candidate`

Two new score keys replace the single `retinal_pocket` key:

```
scores["retinal_pocket_strong"]   # {0.0, 2.0}
scores["retinal_pocket_medium"]   # {0.0, 1.0}
```

Using two keys instead of one preserves the transparency the score dict is designed for — a reader can see exactly which band contributed. `total` is the sum as before.

## 8. Integration with the existing workflow

### 8.1 New CLI subcommand

```
python3 -m opsin_pipeline.cli pocket \
  --scaffold BR \
  --pdb data/structures/1c3w.pdb \
  --pdb-chain A \
  --pdb-offset -2 \                # optional
  --cutoff 6.0 \
  --strong-max 4.0 \
  --medium-max 5.5 \
  --out out/pocket/BR.json
```

Writes `PocketMap` JSON. Exits nonzero on hard errors per §5.4 and §6.2.

### 8.2 Extension of `draft-position-map`

```
python3 -m opsin_pipeline.cli draft-position-map \
  --scaffolds configs/real_scaffolds_reviewed.json \
  --pocket-map out/pocket/BR.json \
  --out out/review/position_map_draft.csv
```

When `--pocket-map` is provided, the CSV draft is pre-annotated with the new columns (§6.3). Human review then proceeds as today.

### 8.3 Scaffold → Candidate path

Scaffold JSON gains an optional field:

```json
{
  "name": "BR",
  "sequence": "...",
  "pocket_map_path": "out/pocket/BR.json",
  "mutable_positions": [...]
}
```

At load time, `load_scaffolds` merges pocket distances from the referenced file into the `MutablePosition` objects (new optional `distance_to_retinal` field). Mutations inherit that distance at generation time. Scoring in §7.2 reads it.

Load-time mismatch: if the pocket map references positions that aren't in `mutable_positions`, that's fine (the human curator decided not to mutate there). If `mutable_positions` references positions absent from the pocket map, those mutations score via the reason-string fallback (§7.3).

### 8.4 Existing code unchanged

`generate.py`, `diversify.py`, `calibration.py`, `report.py` require no changes. The pocket data flows through the existing `MutablePosition` → `Mutation` → `Candidate` → `rank_candidates` path transparently.

## 9. Failure modes — explicit table

| Condition | Classification | Action |
|---|---|---|
| PDB file missing | hard | `FileNotFoundError` |
| PDB unparseable (malformed `ATOM`/`HETATM` line) | hard | `ValueError` with line number |
| Zero retinal-like ligands matched | hard | `NoRetinalLigandError` listing HETATM codes seen |
| ≥2 retinal ligands, no `--pdb-chain` | soft | warn, pick first by (chain, resnum); emit all matches |
| Ligand heavy-atom count out of range (<15 or >25) | hard | `LigandShapeError` unless `--allow-weird-ligand` |
| Insertion code present but `--pdb-offset` used | hard | `NumberingError` — force user to switch to `--pdb-mapping` |
| PDB residue AA ≠ scaffold AA at mapped `seq_index` | soft by default | record `mapping_note`, `review_needed=true`; hard with `--strict-mapping` |
| Chain requested via `--pdb-chain X` absent | hard | `MissingChainError` |
| Alternate-location atoms | silent | keep highest occupancy (ties: alt_loc `A`) |
| NMR multi-model file | silent | take model 1, warn in log |
| Hydrogen atoms in file | silent | skip all `H` element atoms |

Everything listed as "hard" exits the CLI with a nonzero status and a message including the failing file path. Nothing silently produces an empty pocket map — if the parser runs to completion with zero pocket residues, that's itself an error (`EmptyPocketError`).

## 10. Test fixtures

Two PDB files, both committed to the repo under `tests/fixtures/structures/`.

### 10.1 `synthetic_mini.pdb`

Hand-authored. Eight ATOM records forming a short polypeptide, plus one HETATM group labeled `RET` with ten atoms at known coordinates, such that:

- Exactly three residues fall within 4 Å of at least one retinal atom (`strong`).
- Two more residues fall within 4–5.5 Å (`medium`).
- The remaining three are >5.5 Å (`none`).

All coordinates are picked so distances are whole-number verifiable by hand. Tests assert the parser finds exactly the expected bands.

### 10.2 `br_1c3w_strip.pdb`

Real bacteriorhodopsin 1C3W, stripped to: retinal (`RET` or `LYR` depending on how we strip it), K216 (retinal-binding lysine), and all residues within 8 Å of the retinal. Typical size: ~25 residues plus the ligand, well under 10 KB. Used for:

- Ligand-detection smoke test on a real opsin (including the `LYR` split if applicable).
- Band-assignment smoke test against literature expectations: T89, D85, Y185, W86, W182 must come out as `strong` or `medium`, never `none`.
- Numbering reconciliation: scaffold = full BR UniProt sequence; `--pdb-offset` known from the PDB REMARK 465 / DBREF lines; verification round-trips residue identities.

We do **not** commit full-size real PDBs. Anyone running against full structures points `--pdb` at their own local files.

### 10.3 Unit tests to write

- `test_pdb_parser.py`
  - parses ATOM and HETATM lines with fixed-column widths correctly, including negative resnums and insertion codes.
  - rejects malformed lines with a line-numbered error.
  - skips hydrogens, respects alt-loc selection.
- `test_ligands.py`
  - whitelist match / deny-list reject.
  - `LYR` split into Lys atoms vs retinal atoms.
  - heavy-atom-count sanity check fires on a `BCR`-sized group even if mislabeled.
  - `NoRetinalLigandError` with helpful observed-HETATM list.
- `test_pocket_distance.py`
  - synthetic fixture: exact band counts, min-distance values match hand calculation.
  - real fixture: residues known to be in BR pocket come out strong/medium.
  - `--pdb-offset` mode: `seq_index = pdb_resnum + offset`, AA verification fires on a crafted mismatch.
  - `--pdb-mapping` mode: residues in the map get `seq_index`; residues outside get `review_needed=true`.
  - insertion code + `--pdb-offset` raises `NumberingError`.
- `test_score_graded_pocket.py`
  - scaffold with pocket map: `K→F` at a `strong` position scores +2, at a `medium` position scores +1.
  - scaffold without pocket map: legacy reason-string rule still applies.
  - two mutations, one `strong`, one `medium`: candidate score uses `max`, i.e., +2.
- `test_cli_pocket.py`
  - end-to-end: `pocket` subcommand produces a `PocketMap` JSON matching fixtures.
  - `draft-position-map --pocket-map` produces the enriched CSV.

Target: add roughly 8–12 new tests, keeping the 23/23 → 30+/30+ trajectory.

## 11. Calibration-set dependency

**Lock this in the spec even though it's a separate deliverable.** The graded scoring change is meaningless without a real calibration set — the only way to tell whether graded > binary is to run `evaluate_ranking` against known literature mutants and see AUROC move. Current synthetic calibration cannot show improvement because every synthetic "useful" and "disruptive" sits at the same position.

Commitment: before merging the pocket parser to `main`, land a calibration set of **≥20 curated literature mutations** in `configs/calibration_literature.json` (BR D85N/T89A/Y185F, ChR C128T/H134R/E123T/L132C, Chrimson variants, RhGC published shifts, a handful of reported loss-of-function controls). The pocket parser PR includes a before/after AUROC comparison on this set as evidence the change helps.

If curating the calibration set drags on, the parser branch still lands on its own merits (geometric evidence in the review workflow is useful independent of scoring), but the graded scoring toggle stays off-by-default until the calibration set is in.

## 12. Open questions / decisions to confirm before coding

1. **Band points.** Spec uses `+2.0 / +1.0 / 0.0`. Drop to `+1.5 / +0.75 / 0.0` to reduce the pocket signal's dominance over family/phenotype signals? Defer to calibration results.
2. **Max vs sum.** §7.2 uses `max` across mutations. For Hamming-3+ candidates we may want sum-with-diminishing-returns (e.g., `band_points_first + 0.5 * band_points_second + ...`). Defer to calibration results.
3. **Strict vs permissive AA mismatch.** Default is warn. Any scaffold where the crystal was of a mutant construct (common) will produce warnings. If the noise is high we flip to strict + a whitelist of known construct mutations. Start permissive; revisit after first real scaffold run.
4. **`LYR` pocket-center residue.** Always emitted with band `strong`, distance 0. Do we want a distinct `band="schiff_base"` to make the covalent linkage explicit in reports? Nice-to-have; tag as a follow-up.
5. **mmCIF timing.** RCSB is mmCIF-first now. Defer for this chunk, but expect a follow-up within 1–2 iterations. The `structure/` module layout leaves room for `cif.py` next to `pdb.py`.
6. **Pairwise-alignment fallback for numbering.** Explicitly deferred (§6.2 rule 4). If this becomes a recurring request, use a minimal Needleman–Wunsch in pure Python (~30 lines) or pull in `parasail`/`biopython` at that point.

## 13. What this chunk does NOT do

- Does not move `mutable`/`protected` decisions away from human review.
- Does not couple the pipeline to a specific PDB — every scaffold either points at a local file or continues working via the reason-string fallback.
- Does not add a scoring tunable that hasn't been calibrated. The 4.0 / 5.5 thresholds and 2.0 / 1.0 points are starting values; the calibration set from §11 is the evidence that any of these numbers are right.
- Does not fetch over the network. No `requests`, no RCSB API calls.
- Does not break the existing 23 tests. Every change is additive or backward-compat.

## 14. Branch plan

- Branch: `structure-pocket-parser` (this branch).
- First commit: this spec, for review.
- Subsequent commits: parser, ligand handling, pocket computation, position-map integration, scoring integration, tests, one real PDB fixture — ideally one commit per subsystem so review stays tractable.
- Merge target: `main`, via PR. Do not merge until the calibration dependency in §11 is either satisfied or explicitly descoped to a follow-up branch.

## 15. Review checklist (fill in before implementation starts)

- [ ] §5 whitelist covers all retinal codes present in the PDB files we care about — anyone reviewing has checked the 10 most-used opsin PDBs against the list.
- [ ] §6 numbering rules are acceptable; in particular, §6.2 rule 1 (parser never silently emits scaffold `seq_index`) is non-negotiable.
- [ ] §7 thresholds and scoring scheme are accepted as starting values; calibration set (§11) will update them.
- [ ] §9 failure table classifies each condition correctly; hard-error cases fail loudly, soft-warning cases produce actionable messages.
- [ ] §10 fixtures plan is sufficient — no need for a third fixture.
- [ ] §11 calibration commitment is acceptable; either a parallel chunk or a blocker before merge.
