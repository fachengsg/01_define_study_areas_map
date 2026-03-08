## =============================================================================
## 09_Integrated_Pre_Analysis_and_Zonation  (REVISED)
##
## Goal:
##   - Summarize taxonomic/spatial coverage (Region x Group)
##   - Diagnose predictors (SST, SFT, Depth): range + Pearson correlations
##   - Define depth zonation schemes (Primary, Sens_A, Sens_B)
##   - Sensitivity testing: binomial GLMs (ornamented vs non-ornamented) + AIC
##   - Visualize robustness of ornamentation patterns across zonations
##
## Notes (alignment with manuscript):
##   - Response is binary: ornamented (1) vs non-ornamented (0)
##   - "a" is treated as non-ornamented (smooth/commarginal); all other states = ornamented
##   - Analyses are on 0.01° spatially thinned grid-cell occupancies (record-level)
## =============================================================================

library(data.table)
library(ggplot2)

# Ensure dt_final is a data.table
setDT(dt_final)

# Optional: enforce target regions/groups if you want strict consistency
TARGET_REGIONS <- c("Atl-North", "Atl-South", "Mediterranean")
TARGET_GROUPS  <- c("Bivalvia", "Brachiopoda")
dt_final <- dt_final[Region %chin% TARGET_REGIONS & Group %chin% TARGET_GROUPS]

# -----------------------------------------------------------------------------
# STEP 1: Taxonomic & Spatial Coverage (Representativeness)
# -----------------------------------------------------------------------------
cat("--- 1. Taxonomic & Spatial Coverage (Region x Group) ---\n")
coverage_summary <- dt_final[, .(
  n_records = .N,
  n_species = uniqueN(species)
), by = .(Group, Region)][order(Group, Region)]

print(coverage_summary)

# -----------------------------------------------------------------------------
# STEP 2: Trait & Environmental Diagnostics
# -----------------------------------------------------------------------------
cat("\n--- 2. Ornamentation Trait Distribution (raw categories) ---\n")
trait_summary <- dt_final[, .(Count = .N), by = .(Group, Ornament)]
trait_summary[, Proportion := round(Count / sum(Count), 3), by = Group]
print(trait_summary[order(Group, -Count)])

# Binary ornamentation (required for binomial GLMs)
# Assumption consistent with manuscript coding: "a" = smooth/non-ornamented; others = ornamented
dt_final[, Orn_bin := fifelse(Ornament == "a", 0L, 1L)]

cat("\n--- 2b. Binary Ornamentation Distribution (0=non-ornamented, 1=ornamented) ---\n")
trait_bin_summary <- dt_final[, .(Count = .N), by = .(Group, Orn_bin)]
trait_bin_summary[, Proportion := round(Count / sum(Count), 3), by = Group]
print(trait_bin_summary[order(Group, Orn_bin)])

cat("\n--- 3. Environmental Variable Summary (SST, SFT, Depth) ---\n")
env_vars <- c("SST", "SFT", "Depth")
missing_env <- setdiff(env_vars, names(dt_final))
if (length(missing_env) > 0) stop("[STOP] Missing environmental columns in dt_final: ", paste(missing_env, collapse = ", "))

print(summary(dt_final[, ..env_vars]))

cat("\n--- 4. Collinearity Check (Pearson Correlation) ---\n")
cor_matrix <- cor(dt_final[, ..env_vars], use = "complete.obs", method = "pearson")
print(round(cor_matrix, 2))

# Temperature gradient (optional diagnostic)
dt_final[, Temp_Gradient := SST - SFT]

# -----------------------------------------------------------------------------
# STEP 3: Depth Zonation & Sensitivity Setup
# -----------------------------------------------------------------------------
cat("\n--- 5. Depth Zonation Assignment ---\n")

dt_final[, `:=`(
  Zone_Primary = fcase(
    Depth <= 20,                "Coastal (0-20m)",
    Depth > 20  & Depth <= 80,  "Inner Shelf (20-80m)",
    Depth > 80  & Depth <= 200, "Open Shelf (80-200m)",
    Depth > 200,                "Beyond Shelf (>200m)",
    default = NA_character_
  ),
  Zone_Sens_A = fcase(
    Depth <= 20,                "Coastal (0-20m)",
    Depth > 20  & Depth <= 60,  "Inner Shelf (20-60m)",
    Depth > 60  & Depth <= 200, "Open Shelf (60-200m)",
    Depth > 200,                "Beyond Shelf (>200m)",
    default = NA_character_
  ),
  Zone_Sens_B = fcase(
    Depth <= 50,                "Zone 1 (0-50m)",
    Depth > 50  & Depth <= 100, "Zone 2 (50-100m)",
    Depth > 100 & Depth <= 200, "Zone 3 (100-200m)",
    Depth > 200,                "Beyond Shelf (>200m)",
    default = NA_character_
  )
)]

# Drop any NA zones (should be rare if Depth is present)
dt_final <- dt_final[!is.na(Zone_Primary) & !is.na(Zone_Sens_A) & !is.na(Zone_Sens_B)]

# Set ordered factor for consistent plotting (Primary only; others are optional)
zone_levels_primary <- c("Coastal (0-20m)", "Inner Shelf (20-80m)", "Open Shelf (80-200m)", "Beyond Shelf (>200m)")
dt_final[, Zone_Primary := factor(Zone_Primary, levels = zone_levels_primary)]

# -----------------------------------------------------------------------------
# STEP 4: Sensitivity Testing (GLMs + AIC)
# -----------------------------------------------------------------------------
cat("\n--- 6. Sensitivity Testing: GLMs (Orn_bin ~ Depth Zonation) ---\n")

model_primary <- glm(Orn_bin ~ Zone_Primary, data = dt_final, family = binomial)
model_sens_a  <- glm(Orn_bin ~ Zone_Sens_A,  data = dt_final, family = binomial)
model_sens_b  <- glm(Orn_bin ~ Zone_Sens_B,  data = dt_final, family = binomial)

aic_results <- data.frame(
  Zonation = c("Primary (20-80-200)", "Sens_A (20-60-200)", "Sens_B (50-100-200)"),
  AIC_Value = c(AIC(model_primary), AIC(model_sens_a), AIC(model_sens_b))
)

cat("\n--- Model Explanatory Power (Lower AIC is better) ---\n")
print(aic_results[order(aic_results$AIC_Value), ])

# -----------------------------------------------------------------------------
# STEP 5: Visualization of Robustness (Proportion ornamented by zone)
# -----------------------------------------------------------------------------
cat("\n--- 7. Plot: Proportion ornamented by depth zone (three schemes) ---\n")

sens_summary <- rbind(
  dt_final[, .(Prop = mean(Orn_bin), Def = "Primary"), by = .(Zone = as.character(Zone_Primary))],
  dt_final[, .(Prop = mean(Orn_bin), Def = "Sens_A"),  by = .(Zone = as.character(Zone_Sens_A))],
  dt_final[, .(Prop = mean(Orn_bin), Def = "Sens_B"),  by = .(Zone = as.character(Zone_Sens_B))],
  fill = TRUE
)

# -----------------------------------------------------------------------------
# STEP 6: Save cleaned dataset for downstream modelling
# -----------------------------------------------------------------------------
saveRDS(dt_final, "09_final_research_cleaned.rds")
cat("\nSaved: 09_final_research_cleaned.rds\n")












## =============================================================================
## 09_Integrated_Pre_Analysis_and_Zonation
##
## Goal: Clean spatial outliers, summarize data, and define depth-based zones.
## =============================================================================

library(data.table)
library(ggplot2)

# Ensure dt_final is a data.table
setDT(dt_final)


# -----------------------------------------------------------------------------
# STEP 1: Taxonomic & Spatial Coverage (Representativeness)
# -----------------------------------------------------------------------------
cat("--- 1. Taxonomic & Spatial Coverage ---\n")
coverage_summary <- dt_final[, .(
  Total_Records = .N,
  Unique_Species = uniqueN(species)
), by = .(Group, Region)][order(Group, Region)]

print(coverage_summary)

# -----------------------------------------------------------------------------
# STEP 2: Trait & Environmental Diagnostics
# -----------------------------------------------------------------------------
cat("\n--- 2. Ornamentation Trait Distribution ---\n")
trait_summary <- dt_final[, .(Count = .N), by = .(Group, Ornament)]
trait_summary[, Proportion := round(Count / sum(Count), 3), by = Group]
print(trait_summary[order(Group, -Count)])

cat("\n--- 3. Environmental Variable Summary ---\n")
env_vars <- c("SST", "SFT", "Depth")
print(summary(dt_final[, ..env_vars]))
cor_matrix <- cor(dt_final[, ..env_vars], use = "complete.obs", method = "pearson")
print(round(cor_matrix, 2))

# Calculate Temperature Gradient (Stratification indicator)
dt_final[, Temp_Gradient := SST - SFT]

# -----------------------------------------------------------------------------
# STEP 3: Depth Zonation & Sensitivity Setup
# -----------------------------------------------------------------------------
# Define Primary (20-80-200m) and Alternative (Sensitivity) thresholds
dt_final[, `:=`(
  Zone_Primary = fcase(
    Depth <= 20,                "Coastal (0-20m)",
    Depth > 20  & Depth <= 80,  "Inner Shelf (20-80m)",
    Depth > 80  & Depth <= 200, "Open Shelf (80-200m)",
    Depth > 200,                "Beyond Shelf (>200m)"
  ),
  Zone_Sens_A = fcase(
    Depth <= 20,                "Coastal (0-20m)",
    Depth > 20  & Depth <= 60,  "Inner Shelf (20-60m)",
    Depth > 60  & Depth <= 200, "Open Shelf (60-200m)",
    Depth > 200,                "Beyond Shelf (>200m)"
  ),
  Zone_Sens_B = fcase(
    Depth <= 50,                "Coastal/Inner (0-50m)",
    Depth > 50  & Depth <= 100, "Mid Shelf (50-100m)",
    Depth > 100 & Depth <= 200, "Outer Shelf (100-200m)",
    Depth > 200,                "Beyond Shelf (>200m)"
  )
)]

# Set as ordered factor for consistent plotting
zone_levels <- c("Coastal (0-20m)", "Inner Shelf (20-80m)", "Open Shelf (80-200m)", "Beyond Shelf (>200m)")
dt_final[, Zone_Primary := factor(Zone_Primary, levels = zone_levels)]

# -----------------------------------------------------------------------------
# STEP 4: Sensitivity Testing (AIC & Robustness)
# -----------------------------------------------------------------------------
# Test explanatory power (AIC) of different zonations
model_primary <- glm(as.factor(Ornament) ~ Zone_Primary, data = dt_final, family = binomial)
model_sens_a  <- glm(as.factor(Ornament) ~ Zone_Sens_A,  data = dt_final, family = binomial)
model_sens_b  <- glm(as.factor(Ornament) ~ Zone_Sens_B,  data = dt_final, family = binomial)

aic_results <- data.frame(
  Zonation = c("Primary (20-80-200)", "Sens_A (20-60-200)", "Sens_B (50-100-200)"),
  AIC_Value = c(AIC(model_primary), AIC(model_sens_a), AIC(model_sens_b))
)

cat("\n--- Model Explanatory Power (Lower AIC is better) ---\n")
print(aic_results[order(aic_results$AIC_Value), ])

# Visualize Trait Trends across definitions
sens_summary <- rbind(
  dt_final[, .(Prop = sum(Ornament == "a")/.N, Def = "Primary"), by = .(Zone = Zone_Primary)],
  dt_final[, .(Prop = sum(Ornament == "a")/.N, Def = "Sens_A"),  by = .(Zone = Zone_Sens_A)],
  dt_final[, .(Prop = sum(Ornament == "a")/.N, Def = "Sens_B"),  by = .(Zone = Zone_Sens_B)],
  fill = TRUE
)


# Save clean dataset
saveRDS(dt_final, "09_final_research_cleaned.rds")


### prepare Table 1 for the TEXT
tab_table1 <- dt_final[, .(
  n_records = .N,
  n_species = uniqueN(species),
  prop_ornamented = round(mean(Orn_bin), 3),
  median_depth = round(median(Depth), 1),
  median_SFT = round(median(SFT), 2)
), by = .(Region, Group)][order(Region, Group)]

print(tab_table1)
fwrite(tab_table1, file.path("outputs", "Table1_final_sample_size_region_group.csv"))