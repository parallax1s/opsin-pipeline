# Stage 3 — PLM Plausibility Adapter — Spec

Status: **draft (pre-implementation)** — lock before coding.

## 1. Motivation

Stage 2's graded pocket scoring was a clean negative result on the literature calibration (see `rerun_report_graded.md` on `main`): legacy AUROC 0.4588, graded AUROC 0.4176. Every calibration entry is a pocket residue, so both useful (`D85N`, `T89A`) and disruptive (`D85A`, `K216A`) land in the `strong` band and get identical +2 bonuses. The position signal alone cannot separate them — what changes between useful and disruptive at the same position is the **chemistry of the substitution**: `N` preserves a polar H-bond partner; `A` removes the counterion entirely.

A protein language model (PLM) trained on natural protein sequences captures exactly this distinction: `log P(N | context)` at BR position 85 is much higher than `log P(A | context)` because many homologs have `N` or `D` there while very few have `A`. That's the chemistry signal Stage 2's graded pocket bands are blind to.

Stage 3 adds a thin PLM adapter that emits per-mutation log-likelihood deltas. **No scoring changes in Stage 3** — this is evidence plumbing only, same pattern Stage 1 followed for pocket evidence. Scoring that combines position × chemistry (pocket band × PLM delta) lands in Stage 4, gated behind a flag, and only flips on once calibration AUROC actually moves.

## 2. Non-goals

- **No hard torch dependency.** The pipeline's existing 92 tests must keep passing without torch / transformers installed. PLM adapter imports are lazy and guarded.
- **No scoring changes.** Mutations gain a new optional field; scorer ignores it for now.
- **No model training.** We use pretrained ESM2 weights (HuggingFace Hub) via `transformers`.
- **No fine-tuning for opsins.** We evaluate the stock ESM2 prior. Fine-tuning on opsin-specific data is a later chunk, separate branch.
- **No structure-aware PLMs (SaProt/ESM3).** Interesting but adds complexity; defer. ESM2-sequence-only first.
- **No context-window handling beyond ESM2's native 1024.** Scaffold sequences >1024 aa (e.g. `ChR2_Cre_full` at 737, `CaRhGC_full` at 541) fit comfortably; if a future scaffold exceeds 1024, hard-error.
- **No batch-of-scaffolds optimization.** One scaffold per invocation; good enough for 4–6 scaffolds per calibration run.

## 3. Module layout

```
opsin_pipeline/
  plm/
    __init__.py
    predictions.py     # PLMPrediction, PLMPredictionSet data model + JSON I/O
    scorer.py          # abstract PLMScorer, MockPLMScorer
    esm.py             # ESM2Scorer (imports torch + transformers lazily)
  ingest.py            # gains plm_predictions_path merge (mirrors pocket_map_path)
  schemas.py           # MutablePosition / Mutation gain optional plm_log_likelihood_delta

tests/
  test_plm_predictions.py     # MockPLMScorer, JSON round-trip, ingest merge
  test_plm_cli.py             # `cli plm` subcommand end-to-end using mock scorer
  fixtures/
    plm/
      mock_scores_BR.json     # hand-authored mock predictions for BR test
```

No new *always-required* third-party dependencies. `torch` + `transformers` optional, imported only inside `esm.py` and only when a non-Mock scorer is invoked.

## 4. Data model

### 4.1 `PLMPrediction` (per mutation)

```python
@dataclass(frozen=True)
class PLMPrediction:
    position: int                      # scaffold position (1-indexed)
    from_aa: str                       # one-letter
    to_aa: str                         # one-letter
    log_likelihood_delta: float        # log P(to_aa|ctx) - log P(from_aa|ctx)
    model_id: str                      # "esm2_t12_35M_UR50D", "mock", etc.
```

Positive delta → mutation is *more* plausible than wild-type at that position (e.g. mutation to a commonly-seen homolog residue). Negative → less plausible. Stock ESM2 deltas typically range ~[-15, +5].

### 4.2 `PLMPredictionSet` (per scaffold)

```python
@dataclass(frozen=True)
class PLMPredictionSet:
    scaffold_name: str
    sequence_sha256: str               # so rerunning on a changed scaffold sequence is detectable
    model_id: str
    predictions: list[PLMPrediction]
```

Serialized to `configs/plm_predictions/<scaffold>.json`, mirroring `configs/pocket_maps/`.

### 4.3 Schema additions (spec §8.4 pattern)

```python
@dataclass(frozen=True)
class MutablePosition:
    ...
    plm_log_likelihood_deltas: dict[str, float] | None = None   # { to_aa: delta }

@dataclass(frozen=True)
class Mutation:
    ...
    plm_log_likelihood_delta: float | None = None   # the single delta for this (from, to)
```

`MutablePosition` holds a dict because one position can have N allowed target amino acids; each mutation generated from it picks the specific delta.

### 4.4 `Scaffold` addition

```python
@dataclass(frozen=True)
class Scaffold:
    ...
    plm_predictions_path: str | None = None
```

Ingest loads the JSON, builds `{(position, to_aa): delta}` lookup, and attaches a per-position dict to each `MutablePosition`.

## 5. Model choice

**Default:** `esm2_t12_35M_UR50D` — small ESM2 (~140 MB weights, ~35M params). Runs on CPU in <30s for a 300–700 aa sequence. Sufficient context and vocabulary for opsin-length proteins.

**Test scorer:** `MockPLMScorer` — deterministic pseudo-random deltas keyed on `(position, from_aa, to_aa, seed)`. Same interface as ESM2Scorer. All unit tests use this; no network access, no model download.

**Opt-in alternatives (documented, not wired in Stage 3):**
- `esm2_t30_150M_UR50D` — better quality, ~600 MB
- `esm2_t33_650M_UR50D` — even better, ~2.5 GB, needs GPU for speed
- `facebook/esm-c-600M` — ESM-C (newer, smaller param count for same quality)
- SaProt (structure-aware) — requires residue-level secondary structure tokens; spec §12 deferred

## 6. Computation

### 6.1 Mask-free full-position scoring

For a WT sequence of length L, one forward pass through ESM2 gives logits of shape `[L, vocab_size]`. Softmax per position gives `log P(aa | full WT context)` for all 20 amino acids at every position. Then for any mutation `(position p, from_aa, to_aa)`:

```
delta = log_probs[p][to_aa_index] - log_probs[p][from_aa_index]
```

This is O(L × vocab) space, O(1) per mutation lookup after one forward pass per scaffold. For 20 allowed AAs × ~10 mutable positions × 4 scaffolds = trivial.

**Note on mask-free vs masked pseudo-log-likelihood (PLL):** the standard "ESM pseudo-log-likelihood" masks each position in turn and reads the log prob. That's L forward passes per scaffold — 300× slower. Mask-free uses the full WT context, which overstates confidence (the model already "sees" the WT residue at position p when scoring p), but the *delta* between two AAs at the same position is approximately unchanged because both values share the same overstated confidence. For ranking mutations at the same position — which is exactly what we need — mask-free is fine. The spec mentions this explicitly so anyone reading knows what we're approximating.

### 6.2 Per-mutation deltas we actually compute

Only for mutations the pipeline will generate:

```python
for mutable in scaffold.mutable_positions:
    from_aa = scaffold.residue_at(mutable.position)
    for to_aa in mutable.allowed:
        if to_aa == from_aa: continue          # silent; skipped by generator anyway
        delta = log_probs[mutable.position - 1][aa_to_idx[to_aa]] \
              - log_probs[mutable.position - 1][aa_to_idx[from_aa]]
        ...
```

We don't emit deltas for positions not in `mutable_positions` — keeps the output compact.

### 6.3 WT AA sanity check

Before computing deltas, verify that the tokenized WT sequence matches the scaffold sequence character-for-character. ESM2 uses standard one-letter amino acids plus a small special-token set. Any non-ACDEFGHIKLMNPQRSTVWY character in the scaffold (`X`, `B`, `Z`, `U`, `O`, etc.) is either mapped to `<unk>` by ESM2 or is a construct marker. Treatment:

- Standard 20 AAs → straightforward.
- `X` (unknown) in scaffold → skip any mutation *at* this position; log a warning.
- `U` (selenocysteine), `O` (pyrrolysine) → skip position, log warning. No opsin in our calibration has these but future scaffolds might.

## 7. CLI surface

New subcommand, mirroring `cli pocket`:

```
python3 -m opsin_pipeline.cli plm \
  --scaffold BR_Hsal_BACR \
  --scaffolds configs/calibration_scaffolds.json \
  --model esm2_t12_35M_UR50D \
  --out configs/plm_predictions/BR_Hsal_BACR.json \
  [--device cpu|cuda] \
  [--mock]                    # use MockPLMScorer; no torch required
```

Hard-errors:

- `--model` value not recognized and no `--mock`.
- `torch` / `transformers` not importable when a non-mock model is selected.
- Scaffold name not found in `--scaffolds`.
- Scaffold sequence length > 1024 (ESM2 context limit).
- WT AA outside the standard 20 encountered at a `mutable_positions` entry (hard; alternative is a long-tail of silent skips).

Soft (warn + continue):

- No `mutable_positions` on the scaffold (produces an empty prediction set; callers probably want this).

## 8. Integration with existing workflow

### 8.1 Scaffold JSON gains `plm_predictions_path`

```json
{
  "name": "BR_Hsal_BACR",
  "pocket_map_path": "pocket_maps/BR_Hsal_BACR.json",
  "plm_predictions_path": "plm_predictions/BR_Hsal_BACR.json",
  ...
}
```

### 8.2 Ingest merge

`load_scaffolds` does the same pattern as `pocket_map_path`:

1. Read JSON.
2. Build `{(position, to_aa): delta}` lookup.
3. For each `MutablePosition` at `position`, attach
   `plm_log_likelihood_deltas = {to_aa: delta for to_aa in mutable.allowed}` (only the ones the user declared).

### 8.3 Generate propagation

For each `Mutation` produced: copy the single delta `plm_log_likelihood_deltas[to_aa]` onto `Mutation.plm_log_likelihood_delta`. If the position's dict is `None` (no PLM data for this scaffold) or the to_aa is missing, the `Mutation` keeps `plm_log_likelihood_delta=None`.

### 8.4 CSV column

`ranked_candidates.csv` gains a `plm_log_likelihood_delta` column (blank when no data). Per-candidate value = the *minimum* delta across the candidate's mutations (worst-case plausibility), mirroring how `min_distance_to_retinal_A` is aggregated.

### 8.5 Scoring

**No changes to `score.py` in Stage 3.** Mutations now carry PLM data, but no score key reads it. Stage 4 lands the combined scorer behind `--use-plm-chemistry` (or similar).

## 9. Failure modes

| Condition | Classification | Action |
|---|---|---|
| torch / transformers not installed, `--mock` not passed | hard | `PLMBackendUnavailableError` with install hint |
| `--model` string unrecognized | hard | `PLMModelUnknownError` listing allowed values |
| Scaffold name not in `--scaffolds` | hard | `ScaffoldNotFoundError` |
| Scaffold sequence > 1024 aa | hard | `SequenceTooLongError` |
| WT AA ∉ standard 20 at a `mutable_positions` entry | hard | `UnsupportedResidueError` |
| Downloaded model checksum differs from pinned value | hard | `ModelIntegrityError` |
| `to_aa` in `allowed` is not one of the 20 | hard | `UnsupportedTargetAAError` |
| Scaffold has no `mutable_positions` | soft | warn; write empty prediction set |

Every hard error exits nonzero; none silently produce empty/partial JSON.

## 10. Tests

All in-process, zero network, zero torch dependency:

- `test_plm_predictions.py`
  - `MockPLMScorer` is deterministic for a given seed across runs.
  - JSON round-trip via `predictions_to_dict` / `predictions_from_dict`.
  - `sequence_sha256` changes when the sequence changes.
- `test_plm_cli.py`
  - `cli plm --mock` produces a valid `PLMPredictionSet` JSON.
  - WT AA sanity check fires on a scaffold with `X` at a mutable position.
  - Unrecognized `--model` errors without `--mock`.
- `test_plm_ingest_plumbing.py`
  - `load_scaffolds` with `plm_predictions_path` attaches `plm_log_likelihood_deltas` dict to `MutablePosition`.
  - `generate_candidates` propagates the delta onto the `Mutation`.
  - CSV column `plm_log_likelihood_delta` is populated (minimum across mutations).
- `test_plm_no_torch_path.py`
  - Importing `opsin_pipeline.plm` with torch not installed does not raise.
  - Invoking `MockPLMScorer` with torch absent succeeds.

Optional integration test (skipped by default; runs when `OPSIN_ENABLE_ESM_TESTS=1`):

- `test_plm_esm2_smoke.py` — downloads ESM2-t12 once to a cache, runs a 10-residue dummy scaffold, asserts the output shape and that deltas are finite.

Target: ~12 new tests, keeping the cumulative tests trajectory 92 → 104+.

## 11. Calibration strategy

Stage 3 ships evidence-plumbing only; Stage 4 lands the scoring change. The calibration comparison we care about:

| Scorer | Expected AUROC on literature calibration | Commentary |
|---|---:|---|
| legacy (reason string binary) | 0.4588 (current) | no chemistry |
| graded pocket only | 0.4176 (worse) | position but no chemistry |
| PLM delta alone | **to measure in Stage 4** | chemistry but no position |
| graded pocket × PLM delta | **hypothesis: >0.5** | position × chemistry |

If "PLM delta alone" already moves AUROC above 0.5 on our calibration, the combined signal might not help — we still ship Stage 3 + Stage 4 but document the finding. If neither alone nor combined moves AUROC above 0.5, Stage 4's merge waits on better calibration (more entries, chemically-diverse controls) rather than flipping a default on bad evidence.

## 12. Decisions (locked in this draft)

1. **Default model:** ESM2-t12-35M. Small, fast, good-enough prior for opsin-scale sequences. Upgradeable via `--model` later.
2. **Mask-free vs masked PLL:** mask-free. Fast, and accurate enough for *relative* ranking of mutations at the same position, which is what we're measuring.
3. **Aggregation across mutations in a candidate:** `min` across the candidate's per-mutation deltas (worst-case plausibility). Mirror of how `min_distance_to_retinal_A` is aggregated in the CSV.
4. **Scoring integration:** deferred to Stage 4 behind a flag. Stage 3 is plumbing-only.
5. **Scaffold AA alphabet:** standard 20 only; non-standard residues at mutable positions hard-error.

## 12b. Still open / deferred

- **ESM-C / SaProt alternatives.** Defer to Stage 5+ if ESM2 plateaus.
- **Fine-tuning on opsin MSA.** Separate chunk; would live in its own spec.
- **Long-sequence handling beyond 1024 aa.** Sliding-window or last-layer pooling. No opsin in the calibration needs it; revisit when/if we add a multi-domain scaffold.
- **Combined scorer formula (Stage 4).** Options: `sum(pocket, plm)`, `pocket * sigmoid(plm)`, `max(pocket, plm)`. Decide in Stage 4 based on calibration AUROC.

## 13. What this chunk does NOT do

- Does not change any existing score.
- Does not pull torch into the standard test path.
- Does not require a GPU.
- Does not fetch models over the network in the default test run.
- Does not touch calibration data files (configs/calibration_*).

## 14. Branch plan

- Branch: `stage-3-plm-spec` (this branch).
- First commit: this spec. **No code.**
- Subsequent commits (on a follow-up branch `stage-3-plm-adapter` off this one, or on this branch if spec lands fast):
  - schemas + ingest plumbing + generate propagation + MockPLMScorer + tests
  - `opsin_pipeline/plm/esm.py` with lazy imports + CLI `cli plm` subcommand + `--mock` path
  - real ESM2 smoke test (skipped by default)
  - CSV column + scaffold integration

Merge target: `main`, via PR per chunk.

## 15. Review checklist (fill in before implementation)

- [ ] §5 default model choice (ESM2-t12-35M) is acceptable; upgrade path in §12b OK.
- [ ] §6.1 mask-free vs masked PLL tradeoff explicitly accepted.
- [ ] §8 ingest/generate integration mirrors pocket pattern — no surprises.
- [ ] §9 failure table is strict enough; hard-error conditions all listed.
- [ ] §11 calibration strategy: Stage 4 merges only if AUROC moves; empty-move case explicitly documented.
- [ ] §12.3 `min` aggregation accepted (vs mean or median).
- [ ] §14 branch plan: spec-first, plumbing second, scoring last is the agreed cadence.
