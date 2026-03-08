## ============================================================
## OBIS occurrence download + WoRMS AphiaID matching (SHELLFUTURE)
## ============================================================

# >>> load packages only once (and in a consistent order)
library(robis)
library(dplyr)
library(sf)
library(worrms)
library(purrr)

# >>> basic input checks to avoid silent errors
if (!exists("study_areas")) {
  stop("[STOP] Object 'study_areas' not found. Please load it before running this script.")
}
if (!inherits(study_areas, "sf")) {
  stop("[STOP] 'study_areas' must be an sf object.")
}
if (!("Region" %in% names(study_areas))) {
  stop("[STOP] Column 'Region' not found in 'study_areas'.")
}

# >>> define target groups
target_groups <- c("Bivalvia", "Brachiopoda")

# Create an empty list to store results
all_species_data <- list()

# ------------------------------------------------------------
# 1) Loop over regions and taxa to fetch OBIS occurrences
# ------------------------------------------------------------
for (reg_name in unique(study_areas$Region)) {
  
  # >>> 
  message(sprintf("Fetching region: %s ...", reg_name))
  
  # >>> safer geometry handling:
  # - subset sf
  # - union in case the same Region has multiple polygons
  # - convert to a single WKT string
  region_sf <- study_areas %>% dplyr::filter(Region == reg_name)
  if (nrow(region_sf) == 0) next
  
  wkt_geom <- st_as_text(st_union(st_geometry(region_sf)))
  wkt_geom <- as.character(wkt_geom)[1]
  
  for (taxon in target_groups) {
    
    message(sprintf("  -- Fetching taxon/group: %s", taxon))
    
    # >>> keep tryCatch, but also set verbose=FALSE explicitly
    # NOTE: robis::occurrence can return large data; consider limiting fields to what you really need.
    occ <- tryCatch({
      occurrence(
        scientificname = taxon,
        geometry       = wkt_geom,
        fields         = c("scientificName", "species", "decimalLatitude", "decimalLongitude", "basisOfRecord"),
        verbose        = FALSE
      )
    }, error = function(e) {
      message(sprintf("     !!! Failed or no data returned for: %s", taxon))
      return(data.frame())
    })
    
    if (!is.null(occ) && nrow(occ) > 0) {
      occ$Region <- reg_name
      occ$Group  <- taxon
      all_species_data[[paste(reg_name, taxon, sep = "_")]] <- occ
      
      # >>> 
      message(sprintf("     Retrieved records: %d", nrow(occ)))
    }
  }
}

# ------------------------------------------------------------
# 2) Combine results + quick summary
# ------------------------------------------------------------
# >>> handle the "no data" case BEFORE using final_occ_df later
if (length(all_species_data) > 0) {
  final_occ_df <- bind_rows(all_species_data)
  message("--- Finished fetching all regions ---")
  print(table(final_occ_df$Region, final_occ_df$Group))
} else {
  stop("[STOP] No data retrieved. Please check your internet connection or your region geometries.")
}

# ------------------------------------------------------------
# >>> write CSV first (raw), then a smaller RDS for macroecology
write.csv(final_occ_df, "Raw_Occurrence_Data_v1.csv", row.names = FALSE, na = "")
saveRDS(final_occ_df, "Raw_Occurrence_Data_v1.rds")

# ------------------------------------------------------------
# 3) Species list for the next step (trait/taxonomic queries)
# ------------------------------------------------------------
# >>> ensure final_occ_df exists and is non-empty (already guaranteed by stop() above)

species_list <- final_occ_df %>%
  filter(!is.na(scientificName) & scientificName != "") %>%
  mutate(scientificName = trimws(scientificName)) %>%
  filter(!grepl("\\b(sp\\.|indet\\.|cf\\.|aff\\.)\\b", scientificName, ignore.case = TRUE)) %>%
  distinct(scientificName, Group, .keep_all = TRUE) %>%
  select(scientificName, Group)
  
# ------------------------------------------------------------
# 4) Save raw outputs 
write.csv(species_list, "Species_List_for_Traits_v1.csv", row.names = FALSE)


# --- Future Load Command ---
# >>> correct object name in the comment
# Raw_Occurrence_Data <- readRDS("Raw_Occurrence_Data_v1.rds")

# ------------------------------------------------------------
# Load raw OBIS data
# ------------------------------------------------------------
# 
rds_file <- "Raw_Occurrence_Data_v1.rds"  
if (!file.exists(rds_file)) stop("[STOP] Missing raw RDS file: ", rds_file)  
Raw_Occurrence_Data <- readRDS(rds_file)  

# basic column checks (fail fast)
need_cols <- c("scientificName", "decimalLatitude", "decimalLongitude", "basisOfRecord", "Region", "Group")  
miss_cols <- setdiff(need_cols, names(Raw_Occurrence_Data))  
if (length(miss_cols) > 0) stop("[STOP] Raw data missing columns: ", paste(miss_cols, collapse = ", "))  


## ============================================================
## WoRMS matching (AphiaID) with checkpoint/resume  (SHELLFUTURE)
## ============================================================

# ----------------------------
# Packages
# ----------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(worrms)
  library(purrr)
})

# ------------------------------------------------------------
# Load raw OBIS data
# ------------------------------------------------------------
rds_file <- "Raw_Occurrence_Data_v1.rds"
if (!file.exists(rds_file)) stop("[STOP] Missing raw RDS file: ", rds_file)
Raw_Occurrence_Data <- readRDS(rds_file)

# basic column checks (fail fast)
need_cols <- c("scientificName", "decimalLatitude", "decimalLongitude", "basisOfRecord", "Region", "Group")
miss_cols <- setdiff(need_cols, names(Raw_Occurrence_Data))
if (length(miss_cols) > 0) stop("[STOP] Raw data missing columns: ", paste(miss_cols, collapse = ", "))

## ============================================================
## 5.1 Build species list from raw data (reproducible)
## ============================================================

# clean scientific names (trim, remove author tails, drop sp./indet./cf./aff.)
clean_sciname <- function(x) {
  x <- trimws(x)
  x <- gsub("\\s+", " ", x)
  x <- sub("\\s*\\(.*\\)$", "", x)   # remove trailing "(Author, year)" if present
  x <- sub(",\\s*\\d{4}$", "", x)    # remove trailing ", 1758" if present
  x
}

species_list <- Raw_Occurrence_Data %>%
  mutate(scientificName = clean_sciname(scientificName)) %>%
  filter(!is.na(scientificName) & scientificName != "") %>%
  filter(!grepl("\\b(sp\\.|indet\\.|cf\\.|aff\\.)\\b", scientificName, ignore.case = TRUE)) %>%
  distinct(scientificName, Group, .keep_all = FALSE) %>%
  select(scientificName, Group)

# optional save for trait coding / reproducibility
write.csv(species_list, "Species_List_for_WoRMS_v1.csv", row.names = FALSE)

# ------------------------------------------------------------
# 5.2 Unique names to match in WoRMS
# ------------------------------------------------------------
species_to_match <- species_list %>%
  distinct(scientificName) %>%
  pull(scientificName)

message(sprintf("Total unique scientific names to match in WoRMS: %s",
                format(length(species_to_match), big.mark = ",")))

# Configuration for checkpoints
output_file <- "worms_matches_checkpoint.csv"
batch_size  <- 50  # save progress every 50 names (you can increase to 200 to reduce I/O)

# Helper: safely extract a field from possibly inconsistent WoRMS return
get_first_nonnull <- function(x, keys) {
  for (k in keys) {
    if (!is.null(x[[k]]) && length(x[[k]]) > 0 && !is.na(x[[k]][1])) return(x[[k]][1])
  }
  return(NA)
}

# ------------------------------------------------------------
# Initialize or load progress (with migration/fill for old checkpoint)
# ------------------------------------------------------------
if (file.exists(output_file)) {
  
  # (optional) back up checkpoint once
  tryCatch({
    file.copy(
      output_file,
      paste0("backup_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", output_file),
      overwrite = FALSE
    )
  }, error = function(e) {})
  
  matched_data <- read.csv(output_file, stringsAsFactors = FALSE)
  
  # defensive checks
  if (!("query_name" %in% names(matched_data))) {
    stop("[STOP] Checkpoint file exists but has no 'query_name' column: ", output_file)
  }
  
  # ensure expected columns exist (for older checkpoint files)
  for (nm in c("aphiaid","scientificname","status","match_type",
               "matched_aphiaid","valid_aphiaid","matched_name","valid_name","error_msg")) {
    if (!nm %in% names(matched_data)) matched_data[[nm]] <- NA
  }
  
  # --- MIGRATION / FILL (critical fix) ---
  # fill new fields from old fields where missing
  matched_data <- matched_data %>%
    mutate(
      aphiaid         = suppressWarnings(as.numeric(aphiaid)),
      matched_aphiaid = suppressWarnings(as.numeric(matched_aphiaid)),
      valid_aphiaid   = suppressWarnings(as.numeric(valid_aphiaid)),
      matched_aphiaid = coalesce(matched_aphiaid, aphiaid),
      valid_aphiaid   = coalesce(valid_aphiaid, matched_aphiaid),
      matched_name    = coalesce(as.character(matched_name), as.character(scientificname)),
      valid_name      = coalesce(as.character(valid_name), as.character(matched_name))
    )
  
  # write back migrated checkpoint (so next run is consistent)
  write.csv(matched_data, output_file, row.names = FALSE)
  
  done_species <- unique(matched_data$query_name)
  
  # names not in checkpoint yet
  remaining_species <- setdiff(species_to_match, done_species)
  
  # plus names present but still missing AphiaID OR previously not_found/error
  needs_rerun <- matched_data %>%
    filter(
      is.na(valid_aphiaid) |
        status %in% c("not_found", "not_found_or_error", "empty_query")
    ) %>%
    pull(query_name) %>%
    unique()
  
  needs_rerun <- intersect(needs_rerun, species_to_match)
  remaining_species <- unique(c(remaining_species, needs_rerun))
  
  message(sprintf(
    "Resuming... %d names in checkpoint. %d new + %d to rerun (missing AphiaID / error) = %d total to process.",
    length(done_species),
    length(setdiff(species_to_match, done_species)),
    length(needs_rerun),
    length(remaining_species)
  ))
  
} else {
  
  matched_data <- data.frame(
    query_name      = character(),
    aphiaid         = numeric(),
    scientificname  = character(),
    status          = character(),
    match_type      = character(),
    matched_aphiaid = numeric(),
    valid_aphiaid   = numeric(),
    matched_name    = character(),
    valid_name      = character(),
    error_msg       = character(),
    stringsAsFactors = FALSE
  )
  
  remaining_species <- species_to_match
  message(sprintf("Starting a new matching process for %d names.", length(remaining_species)))
}

# ------------------------------------------------------------
# Core matching function with error handling
# ------------------------------------------------------------
match_worms_safe <- function(name) {
  name <- trimws(name)
  if (is.na(name) || name == "") {
    return(data.frame(
      query_name       = name,
      matched_aphiaid   = NA_real_,
      valid_aphiaid     = NA_real_,
      matched_name      = NA_character_,
      valid_name        = NA_character_,
      status            = "empty_query",
      match_type        = "none",
      error_msg         = NA_character_,
      stringsAsFactors  = FALSE
    ))
  }
  
  # polite API usage (avoid 429 bursts)
  Sys.sleep(runif(1, 0.15, 0.30))
  
  tryCatch({
    res_obj <- worrms::wm_records_taxamatch(name)
    res <- if (is.list(res_obj)) res_obj[[1]] else res_obj
    
    if (is.null(res) || nrow(res) == 0) {
      return(data.frame(
        query_name       = name,
        matched_aphiaid   = NA_real_,
        valid_aphiaid     = NA_real_,
        matched_name      = NA_character_,
        valid_name        = NA_character_,
        status            = "not_found",
        match_type        = "none",
        error_msg         = NA_character_,
        stringsAsFactors  = FALSE
      ))
    }
    
    # take the first/best match row
    res1 <- res[1, , drop = FALSE]
    
    matched_aphia <- get_first_nonnull(res1, c("AphiaID", "AphiaId", "aphiaID", "aphiaid"))
    matched_name  <- get_first_nonnull(res1, c("scientificname", "scientificName"))
    stat          <- get_first_nonnull(res1, c("status", "taxonomicStatus"))
    
    # accepted/valid info when present
    valid_aphia <- get_first_nonnull(res1, c("valid_AphiaID", "validAphiaID", "valid_aphiaid"))
    valid_name  <- get_first_nonnull(res1, c("valid_name", "validName"))
    
    # fall back
    if (is.na(valid_aphia)) valid_aphia <- matched_aphia
    if (is.na(valid_name))  valid_name  <- matched_name
    
    data.frame(
      query_name       = name,
      matched_aphiaid   = suppressWarnings(as.numeric(matched_aphia)),
      valid_aphiaid     = suppressWarnings(as.numeric(valid_aphia)),
      matched_name      = as.character(matched_name),
      valid_name        = as.character(valid_name),
      status            = as.character(stat),
      match_type        = "taxamatch",
      error_msg         = NA_character_,
      stringsAsFactors  = FALSE
    )
  }, error = function(e) {
    data.frame(
      query_name       = name,
      matched_aphiaid   = NA_real_,
      valid_aphiaid     = NA_real_,
      matched_name      = NA_character_,
      valid_name        = NA_character_,
      status            = "not_found_or_error",
      match_type        = "none",
      error_msg         = conditionMessage(e),
      stringsAsFactors  = FALSE
    )
  })
}

# ------------------------------------------------------------
# Execution loop with batch saving
# ------------------------------------------------------------
if (length(remaining_species) > 0) {
  
  for (i in seq(1, length(remaining_species), by = batch_size)) {
    
    end_idx <- min(i + batch_size - 1, length(remaining_species))
    current_batch <- remaining_species[i:end_idx]
    
    message(sprintf("Processing batch %d to %d ...", i, end_idx))
    
    batch_results <- purrr::map_dfr(current_batch, match_worms_safe)
    
    # bind + de-duplicate; prefer records with valid_aphiaid and valid_name
    matched_data <- bind_rows(matched_data, batch_results) %>%
      arrange(
        query_name,
        desc(!is.na(valid_aphiaid)),
        desc(!is.na(valid_name))
      ) %>%
      distinct(query_name, .keep_all = TRUE) %>%
      # keep backward-compatible columns filled where possible
      mutate(
        aphiaid        = coalesce(suppressWarnings(as.numeric(aphiaid)), matched_aphiaid),
        scientificname = coalesce(as.character(scientificname), matched_name)
      )
    
    write.csv(matched_data, output_file, row.names = FALSE)
    message("   Batch saved to checkpoint.")
  }
}

message("Taxonomic matching complete!")

# ------------------------------------------------------------
# Summary of results
# ------------------------------------------------------------
summary_stats <- matched_data %>%
  summarise(
    Total   = n(),
    Success = sum(!is.na(valid_aphiaid)),
    Failed  = sum(is.na(valid_aphiaid)),
    Success_rate = round(100 * Success / Total, 2),
    Matched_rate_by_status = round(100 * mean(status != "not_found_or_error"), 2)
  )
print(summary_stats)

# ------------------------------------------------------------
# Join AphiaID back to your species list and save
# ------------------------------------------------------------
species_with_aphia <- species_list %>%
  left_join(
    matched_data %>%
      select(query_name, matched_aphiaid, valid_aphiaid, matched_name, valid_name, status),
    by = c("scientificName" = "query_name")
  )

write.csv(species_with_aphia, "species_worms_matched_v1.csv", row.names = FALSE)
message("Saved joined species-WoRMS table: species_worms_matched_v1.csv")



### prepare data scrape summary table
library(dplyr)

# ---- 0) safety checks ----
stopifnot(exists("Raw_Occurrence_Data"))
stopifnot(exists("matched_data"))
stopifnot("scientificName" %in% names(Raw_Occurrence_Data))
stopifnot("query_name" %in% names(matched_data))

# ---- 1) clean function (use your existing one if already defined) ----
if (!exists("clean_sciname")) {
  clean_sciname <- function(x) {
    x <- trimws(x)
    x <- gsub("\\s+", " ", x)
    x <- sub("\\s*\\(.*\\)$", "", x)
    x <- sub(",\\s*\\d{4}$", "", x)
    x
  }
}

# ---- 2) RAW counts ----
raw_df <- Raw_Occurrence_Data %>%
  mutate(scientificName_clean = clean_sciname(scientificName))

raw_n_records <- nrow(raw_df)

raw_n_unique_scientificName <- raw_df %>%
  summarise(n = n_distinct(scientificName_clean)) %>% pull(n)

# species-level (binomial/trinomial) and exclude sp./cf./aff./indet. (align with Methods)
raw_specieslevel <- raw_df %>%
  filter(!is.na(scientificName_clean) & scientificName_clean != "") %>%
  filter(!grepl("\\b(sp\\.|indet\\.|cf\\.|aff\\.)\\b", scientificName_clean, ignore.case = TRUE)) %>%
  filter(grepl("^[A-Z][A-Za-z-]+\\s+[a-z-]+(\\s+[a-z-]+)?$", scientificName_clean))

raw_n_specieslevel <- raw_specieslevel %>%
  summarise(n = n_distinct(scientificName_clean)) %>% pull(n)

# ---- 3) WoRMS match counts (species-list level) ----
# species list is what you actually sent to WoRMS (cleaned + filtered)
species_to_match <- raw_specieslevel %>%
  distinct(scientificName_clean) %>%
  pull(scientificName_clean)

worms_total_queries <- length(species_to_match)

# matched_data may include mixed generations; use valid_aphiaid if present, else fall back to aphiaid
id_field <- dplyr::case_when(
  "valid_aphiaid" %in% names(matched_data) && sum(!is.na(matched_data$valid_aphiaid)) > 0 ~ "valid_aphiaid",
  "aphiaid"       %in% names(matched_data) && sum(!is.na(matched_data$aphiaid))       > 0 ~ "aphiaid",
  TRUE ~ NA_character_
)
if (is.na(id_field)) stop("[STOP] No usable AphiaID column found in matched_data.")

# restrict to names in your current species_to_match list (avoid old leftovers)
md <- matched_data %>%
  filter(query_name %in% species_to_match)

worms_n_matched_names <- md %>%
  summarise(n = sum(!is.na(.data[[id_field]]))) %>% pull(n)

worms_n_unmatched_names <- md %>%
  summarise(n = sum(is.na(.data[[id_field]]))) %>% pull(n)

# standardized species count after WoRMS:
# prefer valid_name if available, else scientificname
name_field <- dplyr::case_when(
  "valid_name" %in% names(md) && sum(!is.na(md$valid_name)) > 0 ~ "valid_name",
  "scientificname" %in% names(md) && sum(!is.na(md$scientificname)) > 0 ~ "scientificname",
  TRUE ~ NA_character_
)
if (is.na(name_field)) stop("[STOP] No usable standardized name column found in matched_data.")

worms_n_standardized_species <- md %>%
  filter(!is.na(.data[[id_field]])) %>%
  summarise(n = n_distinct(.data[[name_field]])) %>% pull(n)

# ---- 4) WoRMS match counts (occurrence-record level; drop unmatched) ----
occ_worms <- raw_df %>%
  mutate(scientificName_clean = clean_sciname(scientificName)) %>%
  left_join(md, by = c("scientificName_clean" = "query_name"))

worms_occ_n_records_kept <- occ_worms %>%
  filter(!is.na(.data[[id_field]])) %>%
  summarise(n = n()) %>% pull(n)

worms_occ_n_species_kept <- occ_worms %>%
  filter(!is.na(.data[[id_field]])) %>%
  summarise(n = n_distinct(.data[[name_field]])) %>% pull(n)

# ---- 5) print a compact summary table ----
summary_table <- tibble::tibble(
  Stage = c("RAW (all records)", "RAW (species-level binomial/trinomial)", "WoRMS match (species list)", "WoRMS match (occurrences kept)"),
  n_records = c(raw_n_records, nrow(raw_specieslevel), NA, worms_occ_n_records_kept),
  n_species = c(raw_n_unique_scientificName, raw_n_specieslevel, worms_n_standardized_species, worms_occ_n_species_kept),
  notes = c(
    "All downloaded OBIS rows; species count = unique cleaned scientificName (may include sp./cf./aff.)",
    "Filtered to binomial/trinomial; excluded sp./cf./aff./indet.",
    paste0("Total queries=", worms_total_queries,
           "; matched names=", worms_n_matched_names,
           "; unmatched names=", worms_n_unmatched_names,
           "; species standardized by WoRMS"),
    "Occurrences with successful WoRMS match retained; species count uses standardized WoRMS names"
  )
)

print(summary_table)


