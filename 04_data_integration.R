## =============================================================================
## 04 Integrated Pipeline (SHELLFUTURE)
## Traits + Filtering + Environmental Extraction + QA + Coverage checks
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
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. Load & Clean Trait Data
# -----------------------------------------------------------------------------
# Optimization: Combine reading and cleaning to minimize intermediate variables
# NOTE: CSVs must contain 'species' and 'Ornament'
trait_files <- list(
  Brach = "Supplementary_Material_S2_Brachiopoda.csv",
  Bival = "Supplementary_Material_S1_Bivalviak80.csv"
)

# Use lapply to loop through files and validate column names
trait_list <- lapply(trait_files, function(f) {
  d <- fread(f)
  if (!all(c("species", "Ornament") %in% names(d))) stop("[STOP] Missing columns in: ", f)
  return(d[, .(species, Ornament)])
})

trait_dt <- rbindlist(trait_list, use.names = TRUE)

# **Streamlined Cleaning**: Uniform handling of whitespace, null values, and uniqueness
trait_dt[, `:=`(
  species = trimws(gsub("\\s+", " ", as.character(species))),
  Ornament = trimws(as.character(Ornament))
)]
trait_dt <- trait_dt[nzchar(species) & !is.na(Ornament) & nzchar(Ornament)]
trait_dt <- unique(trait_dt, by = "species") 

# =============================================================================
# NEW: Trait Processing - Mapping Morphotypes to Complexity Indices
# =============================================================================

# Standardize 'Ornament' codes to lowercase to ensure robust matching
trait_dt[, orn_clean := tolower(trimws(as.character(Ornament)))]

# -----------------------------------------------------------------------------
# Define Trait Complexity Index (TCI)
# -----------------------------------------------------------------------------
# Mapping rationale:
# 0: Baseline (No reinforcement)
# 1: 1D Structural Reinforcement (Radial or Concentric ribs)
# 2: 2D Grid Reinforcement (Cancellate or Reticulate mesh)
# 3: 3D Architectural Relief (Spines, Nodes, or Tubercles)
# -----------------------------------------------------------------------------
trait_dt[, TCI := fcase(
  orn_clean == "a", 0L,  # smooth or commarginal lines only
  orn_clean == "b", 1L,  # radial
  orn_clean == "c", 1L,  # concentric
  orn_clean == "e", 2L,  # cancellate, reticulate, or other ornaments
  orn_clean == "d", 3L,  # spines, nodes, or tubercles
  default = NA_integer_
)]

# -----------------------------------------------------------------------------
# Define Binary Trait Variable
# -----------------------------------------------------------------------------
# 0 = non-ornamented; 1 = ornamented (any structural relief)
trait_dt[, Binary_Trait := fifelse(TCI > 0, 1, 0)]

# -----------------------------------------------------------------------------
# Quality Assurance (QA)
# -----------------------------------------------------------------------------
# Check for any unrecognized codes that failed the mapping
if(any(is.na(trait_dt$TCI))) {
  unrecognized_codes <- trait_dt[is.na(TCI), unique(Ornament)]
  warning(paste("Mapping Alert: Unrecognized ornamentation codes found:", 
                paste(unrecognized_codes, collapse=", ")))
}

# Remove temporary cleaning column to keep the dataset tidy
trait_dt[, orn_clean := NULL]

# Proceed to Merge with Occurrence Data...

# -----------------------------------------------------------------------------
# 2. Load Occurrence & Merge
# -----------------------------------------------------------------------------
if (!exists("dedup_dt")) {
  dedup_dt <- readRDS("final_occ_dedup_grid_0.01deg.rds")
}

# **Optimization**: Pre-filter regions and groups to reduce computational load during Merge
dt_integrated <- dedup_dt[Region %chin% TARGET_REGIONS & Group %chin% TARGET_GROUPS]

# QA: Identify species missing Trait coding
missing_traits <- dt_integrated[!species %chin% trait_dt$species, unique(.(Region, Group, species))]
cat("Species missing Ornament coding:", nrow(missing_traits), "\n")
fwrite(missing_traits, file.path(OUT_DIR, "species_missing_ornament.csv"))

# **Core Merge**: Retain only records with valid trait data (Inner Join logic)
dt_final <- merge(dt_integrated, trait_dt, by = "species")

# Now dt_final contains 'Ornament', 'TCI', and 'Binary_Trait' columns ready for GAMs!

# -----------------------------------------------------------------------------
# 3. Environmental Extraction (Bio-ORACLE)
# -----------------------------------------------------------------------------
# **Optimization**: Combine SST, SFT, and Depth extraction to avoid redundant coordinate processing
cat("Loading Bio-ORACLE layers (SST, SFT, Depth)...\n")
env_layers <- c("BO22_tempmean_ss", "BO21_tempmean_bdmean", "BO_bathymean")
env_stack  <- rast(load_layers(env_layers, datadir = "./env_data"))

cat("Extracting environmental data for", nrow(dt_final), "points...\n")
coords <- as.matrix(dt_final[, .(decimalLongitude, decimalLatitude)])
extracted_vals <- as.data.table(terra::extract(env_stack, coords))

# **Safe Assignment and Renaming**
dt_final[, `:=`(
  SST   = extracted_vals[[env_layers[1]]],
  SFT   = extracted_vals[[env_layers[2]]],
  Depth = abs(extracted_vals[[env_layers[3]]]) # Convert to positive depth (m)
)]

# Remove records with missing environmental data
dt_final <- na.omit(dt_final, cols = c("SST", "SFT", "Depth"))

# -----------------------------------------------------------------------------
# 4. Surgical Filter (Excluding Eastern Iceland Outliers)
# -----------------------------------------------------------------------------
n_before <- nrow(dt_final)
dt_final <- dt_final[!(decimalLongitude > -16 & decimalLongitude < -10 & decimalLatitude > 62)]

cat(sprintf("\nIceland filter removed %d records.\n", n_before - nrow(dt_final)))

# -----------------------------------------------------------------------------
# 5. Final QA & Export
# -----------------------------------------------------------------------------
# **Optimization**: Use a single function for all summary reports
summary_report <- function(dt) {
  res <- dt[, .(n_records = .N, n_species = uniqueN(species)), by = .(Region, Group)]
  return(res[order(Region, Group)])
}

final_summary <- summary_report(dt_final)
print(final_summary)
fwrite(final_summary, file.path(OUT_DIR, "summary_FINAL_counts.csv"))

# **Species Count Alert System**
check_counts <- function(dt) {
  counts <- dt[, .(n = uniqueN(species)), by = Group]
  # Define expected ranges (List)
  ref <- list("Bivalvia" = c(240, 320), "Brachiopoda" = c(30, 50))
  
  for (g in names(ref)) {
    val <- counts[Group == g, n]
    if (length(val) == 0 || val < ref[[g]][1] || val > ref[[g]][2]) {
      warning(sprintf("%s species count (%d) is outside expected range!", g, val))
    }
  }
}
check_counts(dt_final)

# **Final Export**
saveRDS(dt_final, "07_database_final_integrated.rds")
fwrite(dt_final,  "07_database_final_integrated.csv")

cat("\n--- Pipeline Complete ---\n")
