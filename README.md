## Reproducible workflow (how to run)

This repository contains the R workflow used to (i) define the study domain and regional partitions, (ii) download and clean OBIS occurrences, (iii) harmonise taxonomy with WoRMS, (iv) apply 0.01° grid-cell thinning (spatial deduplication), (v) integrate environmental predictors (Bio-ORACLE v2.0), and (vi) generate maps/models reported in the manuscript and Supplementary Materials.

### Run order (scripts)
Run scripts in numerical order:

1. `01_define_study_areas_map.R` — Define study domain, regional partitions (Mediterranean, Atl-North, Atl-South).
2. `02_raw_data_scrape.R` — Download OBIS occurrences for target clades and study domain (via `robis` / OBIS API).
3. `03_data_cleaning.R` — Standardise species names using WoRMS; remove unmatched names; keep binomial/trinomial only, 0.01° × 0.01° thinning, K80.
4. `04_data_integration.R` — Traits + Filtering + Environmental Extraction.
5. `05_Exploratory_Data_Analysis.R` — Summarize taxonomic/spatial coverage, diagnose predictors (SST, SFT, Depth).
6. `06_Statistical modelling.R` — predictors (SST, SFT, Depth[m] via log10(Depth+1)), binomial GLMs and predicted probability.
7. `07_Spatial synthesis.R` — summary mapping.
8. `08_Regional_Contingency_Analysis.R` — validation and robustness.
9. `09_Spatial_Autocorrelation.R`

> Tip: Script names may differ slightly between versions; the intended execution order is always **01 → 09**.
