## =============================================================================
## 04-07 Integrated Pipeline (SHELLFUTURE)
## Traits + Filtering + Environmental Extraction + QA + Coverage checks
##
## Goal:
##   (1) Join species-level ornamentation to 0.01°-thinned occurrences (dedup_dt)
##   (2) Extract SST/SFT + Depth (Bio-ORACLE)
##   (3) QA: coverage + final counts (Region x Group) for modelling dataset
## =============================================================================

library(data.table)
library(sdmpredictors)
library(terra)

# -----------------------------------------------------------------------------
# 0. Settings
# -----------------------------------------------------------------------------
TARGET_REGIONS <- c("Atl-North", "Atl-South", "Mediterranean")
TARGET_GROUPS  <- c("Bivalvia", "Brachiopoda")
OUT_DIR <- "outputs"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Helper: summary by Region x Group
summ_rg <- function(x, step_name) {
  x[Region %chin% TARGET_REGIONS & Group %chin% TARGET_GROUPS,
    .(step = step_name, n_records = .N, n_species = uniqueN(species)),
    by = .(Region, Group)
  ][order(Region, Group)]
}

# -----------------------------------------------------------------------------
# 1. Load Occurrence & Trait Data
# -----------------------------------------------------------------------------
# Load the 0.01°-thinned spatial data (dedup_dt)
if (!exists("dedup_dt")) {
  dedup_dt <- readRDS("final_occ_dedup_grid_0.01deg.rds")
}

# Ensure we only work with the 3 target regions + 2 target groups
dedup_dt <- dedup_dt[Region %chin% TARGET_REGIONS & Group %chin% TARGET_GROUPS]

# Load ornamentation trait files
# NOTE: both CSVs MUST contain: species, Ornament
brach_orna <- fread("brachiopoda_orna.csv")
biv_orna   <- fread("Bivalviak80new_orna.csv")

# Basic trait file checks
for (nm in c("species","Ornament")) {
  if (!nm %in% names(brach_orna)) stop("[STOP] brachiopoda_orna.csv missing column: ", nm)
  if (!nm %in% names(biv_orna))   stop("[STOP] Bivalviak80new_orna.csv missing column: ", nm)
}

# -----------------------------------------------------------------------------
# 2. Build the Trait Dictionary
# -----------------------------------------------------------------------------
trait_dt <- rbindlist(list(brach_orna, biv_orna), use.names = TRUE, fill = TRUE)
trait_dt[, species := trimws(gsub("\\s+", " ", as.character(species)))]
trait_dt <- trait_dt[!is.na(species) & nzchar(species)]
trait_dt <- trait_dt[, .(species, Ornament)]
trait_dt <- unique(trait_dt, by = "species")  # enforce uniqueness

# -----------------------------------------------------------------------------
# 2.5 Coverage check + Summary AFTER adding Ornament (post-merge)
# -----------------------------------------------------------------------------
# (A) Species-level coverage: which species in dedup_dt have no Ornament coding?
dedup_species <- unique(dedup_dt[Region %chin% TARGET_REGIONS & Group %chin% TARGET_GROUPS, .(Region, Group, species)])
trait_species <- unique(trait_dt[, .(species)])

missing_traits <- dedup_species[!species %chin% trait_species$species]
cat("Species in dedup_dt missing Ornament coding:", nrow(missing_traits), "\n")
fwrite(missing_traits, file.path(OUT_DIR, "species_missing_ornament.csv"))

# (B) Record/species counts AFTER adding Ornament (i.e., post-merge + Ornament != NA)
# NOTE: n_records here = 0.01° thinned occupancies (species×grid-cells), not raw OBIS observations.
dt_integrated <- merge(dedup_dt, trait_dt, by = "species", all.x = TRUE)  # ensure this exists at this point

dt_integrated[, Ornament := trimws(as.character(Ornament))]

dt_after_orn <- dt_integrated[
  Region %chin% TARGET_REGIONS &
    Group  %chin% TARGET_GROUPS  &
    !is.na(Ornament) &
    nzchar(trimws(as.character(Ornament)))
]

summary_after_orn <- dt_after_orn[, .(
  n_records = .N,
  n_species = uniqueN(species)
), by = .(Region, Group)][order(Region, Group)]

cat("\n--- After adding Ornament (records retained): Region x Group summary ---\n")
print(summary_after_orn)

fwrite(summary_after_orn, file.path(OUT_DIR, "summary_after_ornament_region_group.csv"))

# -----------------------------------------------------------------------------
# 3. Merge and Core Filtering (Traits)
# -----------------------------------------------------------------------------
dt_integrated <- merge(dedup_dt, trait_dt, by = "species", all.x = TRUE)

# Keep only records with Ornament (and only target groups already filtered)
dt_traits <- dt_integrated[!is.na(Ornament)]

# Keep only necessary columns to optimize memory
dt_traits <- dt_traits[, .(species, decimalLatitude, decimalLongitude, Region, Group, Ornament)]


# -----------------------------------------------------------------------------
# 4. Temperature Extraction (Bio-ORACLE)
# -----------------------------------------------------------------------------
cat("Loading Bio-ORACLE temperature layers...\n")
temp_layers <- c("BO22_tempmean_ss", "BO21_tempmean_bdmean")
temp_stack  <- load_layers(temp_layers, datadir = "./env_data")
temp_raster <- rast(temp_stack)

cat("Extracting SST and SFT for", nrow(dt_traits), "points...\n")
coords <- dt_traits[, .(decimalLongitude, decimalLatitude)]
extracted_temp <- as.data.table(terra::extract(temp_raster, as.matrix(coords)))

# Safe assignment using column names
if ("BO22_tempmean_ss" %in% names(extracted_temp)) {
  dt_traits[, SST := extracted_temp$BO22_tempmean_ss]
} else {
  stop("[STOP] Temperature layer missing in extraction output: BO22_tempmean_ss")
}

if ("BO21_tempmean_bdmean" %in% names(extracted_temp)) {
  dt_traits[, SFT := extracted_temp$BO21_tempmean_bdmean]
} else {
  stop("[STOP] Temperature layer missing in extraction output: BO21_tempmean_bdmean")
}

# -----------------------------------------------------------------------------
# 5. Marine Masking / Final Cleaning (Temperature)
# -----------------------------------------------------------------------------
# IMPORTANT: require BOTH SST and SFT
dt_temp <- dt_traits[!is.na(SST) & !is.na(SFT)]

# Save intermediate (optional but useful)
saveRDS(dt_temp, "06_final_integrated_research_data.rds")
fwrite(dt_temp,  "06_final_integrated_research_data.csv")

cat("\n--- Temperature Pipeline Complete ---\n")
cat("Records after traits + temp:", nrow(dt_temp), "\n")
cat("File saved: 06_final_integrated_research_data.rds\n")

# -----------------------------------------------------------------------------
# 7. Water Depth Extraction (Bio-ORACLE bathymetry)
# -----------------------------------------------------------------------------
cat("Loading Bio-ORACLE Bathymetry layer...\n")
depth_layer  <- load_layers("BO_bathymean", datadir = "./env_data")
depth_raster <- rast(depth_layer)

cat("Extracting depth for", nrow(dt_temp), "points...\n")
coords <- dt_temp[, .(decimalLongitude, decimalLatitude)]
extracted_depth <- as.data.table(terra::extract(depth_raster, as.matrix(coords)))

depth_col_name <- names(depth_raster)  # usually "BO_bathymean"
if (!(depth_col_name %in% names(extracted_depth))) {
  stop("[STOP] Depth column missing in extraction output: ", depth_col_name)
}

# Bathymetry often negative; convert to positive Depth (m)
dt_temp[, Depth := abs(extracted_depth[[depth_col_name]])]

# Final QC: drop NA depth
dt_final <- dt_temp[!is.na(Depth)]

cat("\n--- Depth Extraction Summary ---\n")
print(summary(dt_final$Depth))

# -----------------------------------------------------------------------------
# STEP 0: Surgical Filter (Excluding Eastern Iceland Outliers)
# -----------------------------------------------------------------------------
# Remove records between -16 and -10 Longitude and above 62 Latitude
# to focus on the contiguous NE Atlantic–Mediterranean system.
n_before_iceland <- nrow(dt_final)

dt_final <- dt_final[!(decimalLongitude > -16 &
                         decimalLongitude < -10 &
                         decimalLatitude  >  62)]

cat(sprintf("\nIceland outlier filter removed %s records (%.3f%%).\n",
            format(n_before_iceland - nrow(dt_final), big.mark = ","),
            100 * (n_before_iceland - nrow(dt_final)) / n_before_iceland))

# -----------------------------------------------------------------------------
# 8. FINAL QA + Summary counts for modelling dataset (Region x Group)
# -----------------------------------------------------------------------------
tab_final <- summ_rg(dt_final, "after_traits+temp+depth_FINAL")
print(tab_final)
fwrite(tab_final, file.path(OUT_DIR, "summary_FINAL_counts_region_group.csv"))

# Cross-check manuscript-level species counts (final dataset)
species_counts_final <- dt_final[, .(n_species = uniqueN(species)), by = Group]
print(species_counts_final)
fwrite(species_counts_final, file.path(OUT_DIR, "final_species_counts_by_group.csv"))

# Optional: warn if far from expected manuscript values
biv_n <- species_counts_final[Group=="Bivalvia", n_species]
bra_n <- species_counts_final[Group=="Brachiopoda", n_species]

if (length(biv_n)==1 && (biv_n < 240 || biv_n > 320)) {
  warning(sprintf("Bivalvia species in FINAL dataset = %d (manuscript says ~270). Check trait list / K80 definition / name matching.", biv_n))
}
if (length(bra_n)==1 && (bra_n < 30 || bra_n > 50)) {
  warning(sprintf("Brachiopoda species in FINAL dataset = %d (manuscript says 38). Check Ye et al. (2023) list / name matching.", bra_n))
}

# -----------------------------------------------------------------------------
# 9. Save FINAL dataset
# -----------------------------------------------------------------------------
saveRDS(dt_final, "07_database_with_temp_and_depth.rds")
fwrite(dt_final,  "07_database_with_temp_and_depth.csv")
