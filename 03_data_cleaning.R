## =============================================================================
## OBIS (0.01° thinning) + Region-specific K80/K90/K95 for Bivalvia
## + ALL Brachiopoda retained
## + Summaries after coord/region filter and after thinning
## =============================================================================

library(data.table)

# -----------------------------------------------------------------------------
# 1) Parameters
# -----------------------------------------------------------------------------
GRID_RES_DEG <- 0.01
K_LEVELS <- c(0.80, 0.90, 0.95)

TARGET_REGIONS <- c("Atl-North", "Atl-South", "Mediterranean")
TARGET_GROUPS  <- c("Bivalvia", "Brachiopoda")

dir.create("outputs", showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 2) Load raw data (OBIS download output)
# -----------------------------------------------------------------------------
dt <- fread("Raw_Occurrence_Data_v1.csv")

# Keep only target regions/groups early
dt <- dt[Region %chin% TARGET_REGIONS & Group %chin% TARGET_GROUPS]

# A) Basic coordinate filter
dt <- dt[!is.na(decimalLatitude) & !is.na(decimalLongitude)]
dt <- dt[between(decimalLatitude, -90, 90) & between(decimalLongitude, -180, 180)]

# B) Safe de-duplication
has_id <- "id" %in% names(dt)
if (has_id) {
  # if id exists but mostly missing, fallback to a composite key
  if (dt[!is.na(id) & nzchar(as.character(id)), .N] > 0) {
    dt <- unique(dt, by = "id")
  } else {
    dt <- unique(dt, by = c("scientificName","decimalLatitude","decimalLongitude","Region","Group"))
  }
} else {
  dt <- unique(dt, by = c("scientificName","decimalLatitude","decimalLongitude","Region","Group"))
}

# C) Species name standardization (IMPORTANT FIX: do NOT prioritize 'species' column)
dt[, species_std := as.character(scientificName)]
dt[, species_std := trimws(gsub("\\s+", " ", species_std))]

# remove obvious uncertainty tags (align with Methods)
dt <- dt[!grepl("\\b(sp\\.|indet\\.|cf\\.|aff\\.)\\b", species_std, ignore.case = TRUE)]

# remove trailing "(Author, year)" or ", 1758" if present
dt[, species_std := sub("\\s*\\(.*\\)$", "", species_std)]
dt[, species_std := sub(",\\s*\\d{4}$", "", species_std)]
dt[, species_std := trimws(gsub("\\s+", " ", species_std))]

# keep only binomial/trinomial and strip extra tails (e.g., "Genus species something...")
dt[, species_std := sub("^([A-Z][A-Za-z-]+\\s+[a-z][A-Za-z0-9-]+(?:\\s+[a-z][A-Za-z0-9-]+)?).*$", "\\1", species_std)]
dt <- dt[grepl("^[A-Z][A-Za-z-]+\\s+[a-z][A-Za-z0-9-]+(\\s+[a-z][A-Za-z0-9-]+)?$", species_std)]

# -----------------------------------------------------------------------------
# 2.5) Summary AFTER coord/region filtering (Region x Group)
# -----------------------------------------------------------------------------
summary_after_filter <- dt[, .(
  n_records = .N,
  n_species = uniqueN(species_std)
), by = .(Region, Group)][order(Region, Group)]

print(summary_after_filter)
write.csv(summary_after_filter, "summary_after_filter.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 3) 0.01° thinning (grid deduplication)
# -----------------------------------------------------------------------------
cat("Performing 0.01° thinning...\n")

dt[, lat_bin := floor((decimalLatitude  + 90)  / GRID_RES_DEG)]
dt[, lon_bin := floor((decimalLongitude + 180) / GRID_RES_DEG)]

# one record per (Region, Group, species, grid cell)
dedup_dt <- dt[, .(
  decimalLatitude  = median(decimalLatitude),
  decimalLongitude = median(decimalLongitude),
  basisOfRecord    = paste(unique(basisOfRecord), collapse = ";"),
  id               = if (has_id) as.character(id[1]) else NA_character_,
  n_raw_records    = .N
), by = .(Region, Group, species = species_std, lat_bin, lon_bin)]

saveRDS(dedup_dt, sprintf("final_occ_dedup_grid_%sdeg.rds", GRID_RES_DEG))
fwrite(dedup_dt,  sprintf("final_occ_dedup_grid_%sdeg.csv", GRID_RES_DEG))

# -----------------------------------------------------------------------------
# 3.5) Summary AFTER thinning (Region x Group)
# -----------------------------------------------------------------------------
summary_after_thinning <- dedup_dt[, .(
  n_records = .N,              # = occupied species×cells
  n_species = uniqueN(species) # unique species after cleaning
), by = .(Region, Group)][order(Region, Group)]

print(summary_after_thinning)
write.csv(summary_after_thinning, "summary_after_thinning.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 4) Region-specific occupancy (Bivalvia): cells per species per region
#    IMPORTANT: Everything below is based on 0.01° thinning (dedup_dt).
# -----------------------------------------------------------------------------
ss_reg <- dedup_dt[Group == "Bivalvia", .(cells = .N), by = .(Region, species)]
setorder(ss_reg, Region, -cells)

# total cells per region (for shares)
region_tot <- ss_reg[, .(cells_total = sum(cells)), by = Region]

# helper: get region-specific K-threshold selection
get_region_kset <- function(ss_reg_dt, k) {
  ss_reg_dt[, {
    x <- .SD[order(-cells)]
    x[, `:=`(
      rank = seq_len(.N),
      cum_share = cumsum(cells) / sum(cells)
    )]
    k_cut <- which.max(x$cum_share >= k)
    x[1:k_cut]
  }, by = Region]
}

# -----------------------------------------------------------------------------
# 4.1) Build K80/K90/K95 sets + summaries
# -----------------------------------------------------------------------------
k_sets <- list()
k_summ_region <- list()
k_summ_global <- list()

for (k in K_LEVELS) {
  
  k_dt <- get_region_kset(ss_reg, k)
  k_dt[, K := k]
  
  # region-level summary
  reg_sum <- k_dt[, .(
    n_species_selected = uniqueN(species),
    n_rows_region_species = .N,          # Region×species rows (your "before dedup" ~350 idea)
    cells_selected = sum(cells)
  ), by = .(K, Region)]
  
  reg_sum <- merge(reg_sum, region_tot, by = "Region", all.x = TRUE)
  reg_sum[, cells_share := cells_selected / cells_total]
  setorder(reg_sum, K, Region)
  
  # global summary across regions
  global_sum <- k_dt[, .(
    K = k,
    n_rows_region_species = .N,          # ~350 for K80 (depends on data)
    n_species_unique = uniqueN(species), # ~270–280 for K80 (depends on data)
    cells_selected_total = sum(cells)
  )]
  
  k_sets[[as.character(k)]] <- k_dt
  k_summ_region[[as.character(k)]] <- reg_sum
  k_summ_global[[as.character(k)]] <- global_sum
}

summary_k_region <- rbindlist(k_summ_region, use.names = TRUE, fill = TRUE)
summary_k_global <- rbindlist(k_summ_global, use.names = TRUE, fill = TRUE)

cat("\n--- Region-specific K summaries (Bivalvia) ---\n")
print(summary_k_region)

cat("\n--- Global K summaries (Bivalvia; union across regions) ---\n")
print(summary_k_global)

write.csv(summary_k_region, "summary_bivalvia_region_specific_K80_K90_K95_by_region.csv", row.names = FALSE)
write.csv(summary_k_global, "summary_bivalvia_region_specific_K80_K90_K95_global.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 4.2) Export candidate lists (Region-specific K80) + unique union list (for trait coding)
# -----------------------------------------------------------------------------
k80_dt <- k_sets[[as.character(0.80)]]

# (a) Region×species list (this is the "350 lines" style)
biv_k80_by_region <- unique(k80_dt[, .(Region, species)])
biv_k80_by_region[, `:=`(orn_type = NA_character_, orn_notes = NA_character_)]

write.csv(biv_k80_by_region, "bivalvia_k80_species_for_ornamentation_by_region.csv", row.names = FALSE)

# (b) Unique species union across regions (this is the "270 unique species" style)
biv_k80_unique <- unique(k80_dt[, .(species)])
biv_k80_unique[, `:=`(orn_type = NA_character_, orn_notes = NA_character_)]
write.csv(biv_k80_unique, "bivalvia_k80_species_for_ornamentation_unique.csv", row.names = FALSE)

# Also export per-region files (optional)
for (r in unique(biv_k80_by_region$Region)) {
  region_fn <- sprintf("bivalvia_k80_species_%s.csv", gsub("[^A-Za-z0-9]+", "_", r))
  fwrite(biv_k80_by_region[Region == r], region_fn)
}

# -----------------------------------------------------------------------------
# 4.3) Brachiopoda: retain ALL species (based on 0.01° thinning)
# -----------------------------------------------------------------------------
bra_all <- dedup_dt[Group == "Brachiopoda", .(cells = .N), by = .(Region, species)]
setorder(bra_all, Region, -cells)

# unique brachiopod species list for coding
fwrite(unique(bra_all[, .(species)]), "brachiopoda_species_names_unique.csv")

# -----------------------------------------------------------------------------
# 4.4) Combined “core list” for downstream (Bivalvia region-K80 + Brachiopoda all)
# -----------------------------------------------------------------------------
biv_k80_core <- k80_dt[, .(
  Region, Group = "Bivalvia", species,
  cells, rank, cum_share
)]

bra_core <- bra_all[, .(
  Region, Group = "Brachiopoda", species,
  cells, rank = NA_integer_, cum_share = NA_real_
)]

core_list_final <- rbind(biv_k80_core, bra_core, fill = TRUE)
saveRDS(core_list_final, "candidate_species_BivalviaRegionK80_BrachiopodaALL.rds")
fwrite(core_list_final,  "candidate_species_BivalviaRegionK80_BrachiopodaALL.csv")

# -----------------------------------------------------------------------------
# 4.5) SUMMARY after applying region-specific K80 (Bivalvia) + ALL Brachiopoda
#      Report n_records (cells) and n_species by Region x Group
# -----------------------------------------------------------------------------
k80_summary_region_group <- core_list_final[, .(
  n_species = uniqueN(species),
  n_records = sum(cells)   # cells = number of occupied 0.01° grid cells per species
), by = .(Region, Group)][order(Region, Group)]

cat("\n--- After K80 (Bivalvia) + ALL Brachiopoda: Region x Group summary ---\n")
print(k80_summary_region_group)

# save
write.csv(k80_summary_region_group,
          "summary_after_K80_corelist_region_group.csv",
          row.names = FALSE)
