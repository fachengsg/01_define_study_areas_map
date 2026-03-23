## =============================================================================
## 05 Integrated Pre-Analysis and Zonation
##
## Goal:
##   - Summarize taxonomic/spatial coverage
##   - Diagnose predictors and correlations
##   - Define depth zonation schemes (Primary + Sensitivity A/B)
##   - Run sensitivity testing (Binomial GLMs + AIC)
## =============================================================================

library(data.table)
library(ggplot2)

# Ensure dt_final is a data.table and filter for consistency
setDT(dt_final)
TARGET_REGIONS <- c("Atl-North", "Atl-South", "Mediterranean")
TARGET_GROUPS  <- c("Bivalvia", "Brachiopoda")
dt_final <- dt_final[Region %chin% TARGET_REGIONS & Group %chin% TARGET_GROUPS]

# -----------------------------------------------------------------------------
# 1. Taxonomic & Spatial Coverage
# -----------------------------------------------------------------------------
cat("--- 1. Taxonomic & Spatial Coverage (Region x Group) ---\n")
# **Refinement**: Consolidated the two coverage summary blocks into one efficient table
coverage_summary <- dt_final[, .(
  n_records = .N,
  n_species = uniqueN(species)
), by = .(Group, Region)][order(Group, Region)]
print(coverage_summary)

# -----------------------------------------------------------------------------
# 2. Trait & Environmental Diagnostics
# -----------------------------------------------------------------------------
# **Improvement**: Created binary trait (Orn_bin) immediately to avoid redundant calculations
# Coding: "a" = smooth/non-ornamented (0); all others = ornamented (1)
dt_final[, Orn_bin := fcase(Ornament == "a", 0L, default = 1L)]

cat("\n--- 2. Binary Ornamentation Distribution (0=non-orn, 1=orn) ---\n")
trait_summary <- dt_final[, .(Count = .N), by = .(Group, Orn_bin)]
trait_summary[, Proportion := round(Count / sum(Count), 3), by = Group]
print(trait_summary[order(Group, Orn_bin)])

cat("\n--- 3. Environmental Variables & Collinearity ---\n")
env_vars <- c("SST", "SFT", "Depth")
# Calculate Temperature Gradient (Stratification indicator)
dt_final[, Temp_Gradient := SST - SFT]

print(summary(dt_final[, ..env_vars]))
cor_matrix <- cor(dt_final[, ..env_vars], use = "complete.obs")
print(round(cor_matrix, 2))

# -----------------------------------------------------------------------------
# 3. Depth Zonation Assignment
# -----------------------------------------------------------------------------
# **Optimization**: Combined factor conversion with fcase to ensure clean formatting
cat("\n--- 4. Depth Zonation Assignment ---\n")

dt_final[, `:=`(
  Zone_Primary = fcase(
    Depth <= 20, "Coastal (0-20m)",
    Depth <= 80, "Inner Shelf (20-80m)",
    Depth <= 200, "Open Shelf (80-200m)",
    Depth > 200, "Beyond Shelf (>200m)"
  ),
  Zone_Sens_A = fcase(
    Depth <= 20, "Coastal (0-20m)",
    Depth <= 60, "Inner Shelf (20-60m)",
    Depth <= 200, "Open Shelf (60-200m)",
    Depth > 200, "Beyond Shelf (>200m)"
  ),
  Zone_Sens_B = fcase(
    Depth <= 50, "Zone 1 (0-50m)",
    Depth <= 100, "Zone 2 (50-100m)",
    Depth <= 200, "Zone 3 (100-200m)",
    Depth > 200, "Beyond Shelf (>200m)"
  )
)]

# Set factors for plotting consistency
zone_levels <- c("Coastal (0-20m)", "Inner Shelf (20-80m)", "Open Shelf (80-200m)", "Beyond Shelf (>200m)")
dt_final[, Zone_Primary := factor(Zone_Primary, levels = zone_levels)]

# -----------------------------------------------------------------------------
# 4. Sensitivity Testing (GLMs + AIC)
# -----------------------------------------------------------------------------
# **Refinement**: Automated the model comparison to be more concise
cat("\n--- 5. Sensitivity Testing: AIC Comparison ---\n")

models <- list(
  "Primary (20-80-200)" = glm(Orn_bin ~ Zone_Primary, data = dt_final, family = binomial),
  "Sens_A (20-60-200)"  = glm(Orn_bin ~ Zone_Sens_A,  data = dt_final, family = binomial),
  "Sens_B (50-100-200)" = glm(Orn_bin ~ Zone_Sens_B,  data = dt_final, family = binomial)
)

aic_results <- data.table(
  Zonation = names(models),
  AIC_Value = sapply(models, AIC)
)[order(AIC_Value)]
print(aic_results)

# -----------------------------------------------------------------------------
# 5. Output Generation (Table 1 & Clean Data)
# -----------------------------------------------------------------------------
# **Consolidated Export**: Grouped all summary exports into one clear section
cat("\n--- 6. Exporting Final Summaries & Dataset ---\n")

tab_table1 <- dt_final[, .(
  n_records = .N,
  n_species = uniqueN(species),
  prop_ornamented = round(mean(Orn_bin), 3),
  median_depth = round(median(Depth), 1),
  median_SFT = round(median(SFT), 2)
), by = .(Region, Group)][order(Region, Group)]

print(tab_table1)
fwrite(tab_table1, file.path("outputs", "Table1_final_summary.csv"))

# Final Species counts for Manuscript text
cat("Total unique species:", uniqueN(dt_final$species), "\n")
print(dt_final[, .(n_species = uniqueN(species)), by = Group])

# Save final processed dataset
saveRDS(dt_final, "09_final_research_cleaned.rds")
cat("\nProcess Complete. Saved: 09_final_research_cleaned.rds\n")

