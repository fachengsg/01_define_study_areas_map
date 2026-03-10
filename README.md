## Reproducible workflow (how to run)

This repository contains the R workflow used to (i) define the study domain and regional partitions, (ii) download and clean OBIS occurrences, (iii) harmonise taxonomy with WoRMS, (iv) apply 0.01° grid-cell thinning (spatial deduplication), (v) integrate environmental predictors (Bio-ORACLE v2.0), and (vi) generate maps/models reported in the manuscript and Supplementary Materials.

### Run order (scripts)
Run scripts in numerical order:

1. `01_define_study_areas_map.R` — Define study domain, regional partitions (Mediterranean, Atl-North, Atl-South), and export region polygons.
2. `02_download_obis_occurrences.R` — Download OBIS occurrences for target clades and study domain (via `robis` / OBIS API).
3. `03_taxonomy_worms_standardisation.R` — Standardise species names using WoRMS; remove unmatched names; keep binomial/trinomial only.
4. `04_clean_coordinates_filters.R` — Remove invalid coordinates/out-of-domain records; apply additional spatial filters (e.g., Iceland sector exclusion).
5. `05_grid_thinning_0p01deg.R` — 0.01° × 0.01° thinning: collapse repeated occurrences of the same species within the same grid cell to a single presence per region.
6. `06_environment_biooracle_extraction.R` — Extract Bio-ORACLE v2.0 variables (SST, SFT, bathymetry) for thinned occurrences.
7. `07_trait_merge_ornamentation.R` — Merge species-level ornamentation coding; derive binary ornamentation (ornamented vs non-ornamented).
8. `08_models_and_validation.R` — Fit binomial GLMs (additive + clade interactions), z-score sensitivity, geographic proxies (latitude / depth zones), regional contingency, and spatial diagnostics (Moran’s I; RAC if enabled).
9. `09_figures_tables_export.R` — Recreate main figures and supplementary figures/tables used in the manuscript.

> Tip: Script names may differ slightly between versions; the intended execution order is always **01 → 09**.
