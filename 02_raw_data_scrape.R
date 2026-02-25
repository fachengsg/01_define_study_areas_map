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
if (!exists("study_areas_final_v3")) {
  stop("[STOP] Object 'study_areas_final_v3' not found. Please load it before running this script.")
}
if (!inherits(study_areas_final_v3, "sf")) {
  stop("[STOP] 'study_areas_final_v3' must be an sf object.")
}
if (!("Region" %in% names(study_areas_final_v3))) {
  stop("[STOP] Column 'Region' not found in 'study_areas_final_v3'.")
}

# >>> define target groups
target_groups <- c("Bivalvia", "Gastropoda", "Brachiopoda")

# Create an empty list to store results
all_species_data <- list()

# ------------------------------------------------------------
# 1) Loop over regions and taxa to fetch OBIS occurrences
# ------------------------------------------------------------
for (reg_name in unique(study_areas_final_v3$Region)) {
  
  # >>> 
  message(sprintf("Fetching region: %s ...", reg_name))
  
  # >>> safer geometry handling:
  # - subset sf
  # - union in case the same Region has multiple polygons
  # - convert to a single WKT string
  region_sf <- study_areas_final_v3 %>% dplyr::filter(Region == reg_name)
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
# 3) Species list for the next step (trait/taxonomic queries)
# ------------------------------------------------------------
# >>> ensure final_occ_df exists and is non-empty (already guaranteed by stop() above)
species_list <- final_occ_df %>%
  filter(!is.na(species) & species != "") %>%
  distinct(species, .keep_all = TRUE) %>%
  select(species, Group)

# ------------------------------------------------------------
# 4) Save raw outputs
# ------------------------------------------------------------
# >>> UPDATED: write CSV first (raw), then a smaller RDS for macroecology
write.csv(final_occ_df, "Raw_Occurrence_Data_v1.csv", row.names = FALSE)

final_occ_cleaned <- final_occ_df %>%
  select(scientificName, species, Group, Region, decimalLatitude, decimalLongitude, basisOfRecord)

saveRDS(final_occ_cleaned, "SHELLFUTURE_Occurrences_Raw_v1.rds")

# >>> remove hard-coded "2.1 million" and report the real number
message(sprintf("Successfully saved %s records to RDS!",
                format(nrow(final_occ_cleaned), big.mark = ",")))

# --- Future Load Command ---
# >>> correct object name in the comment (cleaned data is what you saved)
# final_occ_cleaned <- readRDS("SHELLFUTURE_Occurrences_Raw_v1.rds")


## ============================================================
## 5) WoRMS matching (AphiaID) with checkpoint/resume
## ============================================================

# >>> prepare unique species vector from species_list (already cleaned)
species_to_match <- species_list %>%
  distinct(species) %>%
  pull(species)

message(sprintf("Total unique species to match in WoRMS: %s",
                format(length(species_to_match), big.mark = ",")))

# Configuration for checkpoints
output_file <- "worms_matches_checkpoint.csv"
batch_size  <- 50  # save progress every 50 species

# >>> helper to safely extract a field from a possibly inconsistent WoRMS return
get_first_nonnull <- function(x, keys) {
  for (k in keys) {
    if (!is.null(x[[k]]) && length(x[[k]]) > 0 && !is.na(x[[k]][1])) return(x[[k]][1])
  }
  return(NA)
}

# Initialize or load progress
if (file.exists(output_file)) {
  matched_data <- read.csv(output_file, stringsAsFactors = FALSE)
  
  # >>> be defensive about column existence (older checkpoints, manual edits, etc.)
  if (!("query_name" %in% names(matched_data))) {
    stop("[STOP] Checkpoint file exists but has no 'query_name' column: ", output_file)
  }
  
  done_species <- unique(matched_data$query_name)
  remaining_species <- setdiff(species_to_match, done_species)
  
  message(sprintf("Resuming... %d species already matched. %d remaining.",
                  length(done_species), length(remaining_species)))
} else {
  matched_data <- data.frame(
    query_name      = character(),
    aphiaid         = numeric(),
    scientificname  = character(),
    status          = character(),
    match_type      = character(),
    stringsAsFactors = FALSE
  )
  remaining_species <- species_to_match
  message(sprintf("Starting a new matching process for %d species.", length(remaining_species)))
}

# Core matching function with error handling
match_worms_safe <- function(name) {
  Sys.sleep(0.2)  # polite API usage (avoid 429 Too Many Requests)
  
  tryCatch({
    res_list <- wm_records_taxamatch(name)
    res <- res_list[[1]]
    
    if (is.null(res) || nrow(res) == 0) {
      return(data.frame(
        query_name     = name,
        aphiaid        = NA,
        scientificname = NA,
        status         = "not_found",
        match_type     = "none",
        stringsAsFactors = FALSE
      ))
    }
    
    # >>> take the first/best match row
    res1 <- res[1, ]
    
    aphia_id <- get_first_nonnull(res1, c("AphiaID", "AphiaId", "aphiaID", "aphiaid"))
    sci_name <- get_first_nonnull(res1, c("scientificname", "scientificName", "valid_name"))
    stat     <- get_first_nonnull(res1, c("status", "taxonomicStatus"))
    
    return(data.frame(
      query_name     = name,
      aphiaid        = suppressWarnings(as.numeric(aphia_id)),
      scientificname = as.character(sci_name),
      status         = as.character(stat),
      match_type     = "exact/fuzzy",
      stringsAsFactors = FALSE
    ))
  }, error = function(e) {
    return(data.frame(
      query_name     = name,
      aphiaid        = NA,
      scientificname = NA,
      status         = "not_found_or_error",
      match_type     = "none",
      stringsAsFactors = FALSE
    ))
  })
}

# Execution loop with batch saving
if (length(remaining_species) > 0) {
  for (i in seq(1, length(remaining_species), by = batch_size)) {
    
    end_idx <- min(i + batch_size - 1, length(remaining_species))
    current_batch <- remaining_species[i:end_idx]
    
    message(sprintf("Processing batch %d to %d ...", i, end_idx))
    
    batch_results <- map_dfr(current_batch, match_worms_safe)
    
    # >>> UPDATED: bind + de-duplicate by query_name (important when resuming)
    matched_data <- bind_rows(matched_data, batch_results) %>%
      distinct(query_name, .keep_all = TRUE)
    
    write.csv(matched_data, output_file, row.names = FALSE)
    message("   Batch saved to checkpoint.")
  }
}

message("Taxonomic matching complete!")

# Summary of results
summary_stats <- matched_data %>%
  summarise(
    Total   = n(),
    Success = sum(!is.na(aphiaid)),
    Failed  = sum(is.na(aphiaid))
  )
print(summary_stats)

# >>> (optional but useful): join AphiaID back to your species list and save
species_with_aphia <- species_list %>%
  left_join(matched_data %>% select(query_name, aphiaid, scientificname, status),
            by = c("species" = "query_name"))

write.csv(species_with_aphia, "species_worms_matched_v1.csv", row.names = FALSE)
message("Saved joined species-WoRMS table: species_worms_matched_v1.csv")
