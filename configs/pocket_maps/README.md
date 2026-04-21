# Calibration PocketMaps

Scaffold-numbered `PocketMap` JSONs for the literature calibration scaffolds
in `../calibration_scaffolds.json`. Produced from RCSB PDBs via the
`opsin_pipeline.cli pocket` subcommand with `--pdb-offset` chosen so PDB
author numbering aligns onto scaffold sequence coordinates (AA verification
clean on all four).

## Regenerate

```bash
# Fetch the structures (not committed)
mkdir -p data/structures
for pdb in 1C3W 6EID 5ZIH 6CSM; do
  curl -s -o "data/structures/$pdb.pdb" "https://files.rcsb.org/download/$pdb.pdb"
done

# BR / bacteriorhodopsin — mature-minus-signal-peptide offset
PYTHONPATH=. python3 -m opsin_pipeline.cli pocket \
  --scaffold BR_Hsal_BACR --pdb data/structures/1C3W.pdb \
  --pdb-chain A --pdb-offset 13 \
  --scaffolds configs/calibration_scaffolds.json \
  --out configs/pocket_maps/BR_Hsal_BACR.json

# ChR2 / Chlamydomonas channelrhodopsin-2
PYTHONPATH=. python3 -m opsin_pipeline.cli pocket \
  --scaffold ChR2_Cre_full --pdb data/structures/6EID.pdb \
  --pdb-chain A --pdb-offset 0 \
  --scaffolds configs/calibration_scaffolds.json \
  --out configs/pocket_maps/ChR2_Cre_full.json

# Chrimson — scaffold trimmed 3 residues relative to PDB
PYTHONPATH=. python3 -m opsin_pipeline.cli pocket \
  --scaffold Chrimson_5ZIH_structure --pdb data/structures/5ZIH.pdb \
  --pdb-chain A --pdb-offset -3 \
  --scaffolds configs/calibration_scaffolds.json \
  --out configs/pocket_maps/Chrimson_5ZIH_structure.json

# GtACR1 — same offset pattern as Chrimson
PYTHONPATH=. python3 -m opsin_pipeline.cli pocket \
  --scaffold GtACR1_6CSM_structure --pdb data/structures/6CSM.pdb \
  --pdb-chain A --pdb-offset -3 \
  --scaffolds configs/calibration_scaffolds.json \
  --out configs/pocket_maps/GtACR1_6CSM_structure.json
```

Each PocketMap records `pdb_sha256`, so rerunning on a different version
of a PDB entry will be visibly detectable.

## What's not here

- `CaRhGC_full` has no published crystal structure as of 2026-04; the
  scaffold stays without a `pocket_map_path`, which means candidates on
  CaRhGC fall back to the reason-string rule under `--graded-pocket`
  (spec §7.3 fallback — `pocket_graded_fallback` tag in output).
- Full PDB files are not committed; they're re-fetchable from RCSB and
  the `pdb_sha256` in each PocketMap lets reviewers verify they got the
  same bytes.
