# Opsin Pipeline Decision Report

- Target family: RhGC
- Target phenotype: spectral_tuning
- Scaffolds screened: 1
- Candidates generated: 32
- Protected-residue violations: 0
- Diversity caps: per_scaffold_cap=20, per_position_cap=5
- Pocket signal: legacy (reason-string binary, flag off)

## Scaffolds

| Scaffold | Family | λmax start (nm) | Protected | Mutable |
|---|---|---:|---:|---:|
| BeRhGC1_full | RhGC | 540.0 | 33 | 4 |

## Generation

- Total candidates: 32

## Top 10 Candidates (diversified)

| Rank | Candidate | Scaffold | Hamming | Score | λmax start (nm) | Tags |
|---:|---|---|---:|---:|---:|---|
| 1 | BeRhGC1_full_p283DtoN | BeRhGC1_full | 1 | 7.8 | 540.0 | assay_ready, retinal_pocket, spectral_tuning, target_family, target_phenotype |
| 2 | BeRhGC1_full_p283DtoE | BeRhGC1_full | 1 | 7.8 | 540.0 | assay_ready, retinal_pocket, spectral_tuning, target_family, target_phenotype |
| 3 | BeRhGC1_full_p272YtoW | BeRhGC1_full | 1 | 7.8 | 540.0 | assay_ready, retinal_pocket, spectral_tuning, target_family, target_phenotype |
| 4 | BeRhGC1_full_p272YtoF | BeRhGC1_full | 1 | 7.8 | 540.0 | assay_ready, retinal_pocket, spectral_tuning, target_family, target_phenotype |
| 5 | BeRhGC1_full_p183WtoY | BeRhGC1_full | 1 | 7.8 | 540.0 | assay_ready, retinal_pocket, spectral_tuning, target_family, target_phenotype |
| 6 | BeRhGC1_full_p183WtoF | BeRhGC1_full | 1 | 7.8 | 540.0 | assay_ready, retinal_pocket, spectral_tuning, target_family, target_phenotype |
| 7 | BeRhGC1_full_p181DtoN | BeRhGC1_full | 1 | 7.8 | 540.0 | assay_ready, retinal_pocket, spectral_tuning, target_family, target_phenotype |
| 8 | BeRhGC1_full_p181DtoE | BeRhGC1_full | 1 | 7.8 | 540.0 | assay_ready, retinal_pocket, spectral_tuning, target_family, target_phenotype |
| 9 | BeRhGC1_full_p272YtoW__p283DtoN | BeRhGC1_full | 2 | 7.6 | 540.0 | assay_ready, retinal_pocket, spectral_tuning, target_family, target_phenotype |
| 10 | BeRhGC1_full_p272YtoW__p283DtoE | BeRhGC1_full | 2 | 7.6 | 540.0 | assay_ready, retinal_pocket, spectral_tuning, target_family, target_phenotype |

## Next Validation Layer

- Run structure and retinal-pocket checks on the top candidates.
- Reserve QM/MM or excited-state calculations for a small shortlist.
- Treat this MVP score as triage, not a final biophysical prediction.
