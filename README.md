## Reproducible workflow (how to run)

This repository contains the complete R workflow used to: (i) define the study domain and regional biogeographic partitions, (ii) download and taxonomically harmonise OBIS occurrence records against WoRMS, (iii) mitigate spatial sampling biases via 0.01° grid-cell deduplication and H3 discrete global grid aggregation, (iv) integrate high-resolution environmental predictors (Bio-ORACLE v2.0), (v) fit assemblage-level Generalized Additive Models (GAMs) and Beta regressions to evaluate environmental drivers, and (vi) rigorously test and control for fine-scale spatial autocorrelation (SAC) using Residual Autocovariate (RAC) models to generate all macroecological maps and statistical results reported in the manuscript.

### Run order (scripts)
Run scripts in numerical order:

1. `01_define_study_areas_map.R` — Define study domain, regional partitions (Mediterranean, Atl-North, Atl-South).
2. `02_raw_data_scrape.R` — Download OBIS occurrences for target clades and study domain (via `robis` / OBIS API).
3. `03_data_cleaning.R` — Standardise species names using WoRMS; remove unmatched names; keep binomial/trinomial only, 0.01° × 0.01° thinning, K80.
4. `04_data_integration.R` — Traits + Filtering + Environmental Extraction.
5. `05_Exploratory_Data_Analysis.R` — Summarise taxonomic and spatial coverage; conduct predictor diagnostics and multicollinearity checks (Spearman's ρ and VIF) for the environmental variables.
6. `06_Spatial_synthesis.R` — Aggregate occurrence records into the H3 discrete global grid system (resolution 4) to map large-scale spatial gradients of ornamentation prevalence.
7. `07_Binomial_GAM_Beta.R` — Fit assemblage-level binomial Generalized Additive Models (GAMs) and complementary Beta regression models; perform deviance partitioning to assess variable importance; compute and visualise marginal effects for temperature and depth.
8. `08_Ornamental_Trend_Morphotype.R` — Analyse raw, unadjusted species- and occurrence-level proportions along environmental gradients; map the continuous ecological turnover of specific bivalve morphotypes (e.g., smooth vs. spinose/cancellate).
9. `09_Spatial_Autocorrelation.R` — Assess fine-scale spatial clustering using global Moran’s I; construct and evaluate Residual Autocovariate (RAC) GAMs to ensure that the core environmental signals are statistically robust to spatial autocorrelation

> Tip: Script names may differ slightly between versions; the intended execution order is always **01 → 09**.
