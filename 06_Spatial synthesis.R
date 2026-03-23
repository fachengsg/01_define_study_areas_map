################################################################################
## Spatial Analysis: Binary Presence vs. Trait Complexity (TCI)
## 
## Goals:
## - Map Binary Ornamentation (%) and Mean TCI across H3 Hexagonal Grids
## - Generate a 2x2 Main Figure for the Manuscript (Presence + Complexity)
################################################################################

library(data.table)
library(h3jsr)
library(sf)
library(ggplot2)
library(viridis)
library(patchwork)
library(scales)
library(rnaturalearth)

# -----------------------------------------------------------------------------
# 1. Initialization & Settings
# -----------------------------------------------------------------------------
OUT_DIR <- "outputs"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

TARGET_REGIONS <- c("Atl-North", "Atl-South", "Mediterranean")
TARGET_GROUPS  <- c("Bivalvia", "Brachiopoda")
H3_RESOLUTION  <- 4  

xlim_use <- c(-32, 38)
ylim_use <- c(8, 70)
x_breaks <- seq(-30, 40, by = 10)
y_breaks <- seq(10, 70, by = 10)

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

# -----------------------------------------------------------------------------
# 2. Data Preparation
# -----------------------------------------------------------------------------
if (!exists("dt_final")) {
  dt_final <- readRDS("09_final_research_cleaned.rds")
}
setDT(dt_final)

# Filter domain and ensure required columns (Binary_Trait and TCI) exist
dt <- dt_final[Region %chin% TARGET_REGIONS & Group %chin% TARGET_GROUPS]

# Compute H3 Index if missing
if (!("h3" %in% names(dt))) {
  dt[, h3 := point_to_cell(cbind(decimalLongitude, decimalLatitude), res = H3_RESOLUTION)]
}

# -----------------------------------------------------------------------------
# 3. Calculate Hexagon Statistics (Updated to include TCI for Occurrences)
# -----------------------------------------------------------------------------
# Species-weighted stats (Main Text logic)
hex_sp <- unique(dt[, .(h3, Group, species, Binary_Trait, TCI)])
hex_stats_sp <- hex_sp[, .(
  prop_orn_sp = mean(Binary_Trait), 
  mean_tci_sp = mean(TCI), 
  n_species = .N
), by = .(h3, Group)]

# Occurrence-weighted stats (Supplementary logic)
# Note: Here we calculate the mean directly from 'dt' (every record counts)
hex_stats_occ <- dt[, .(
  prop_orn_occ = mean(Binary_Trait), 
  mean_tci_occ = mean(TCI), 
  n_records = .N
), by = .(h3, Group)]

# Merge stats and generate spatial SF object
hex_stats <- merge(hex_stats_sp, hex_stats_occ, by = c("h3", "Group"), all = TRUE)
hex_geom  <- cell_to_polygon(unique(hex_stats$h3), simple = FALSE)
hex_sf    <- st_as_sf(data.frame(h3 = unique(hex_stats$h3), geometry = hex_geom))
hex_sf    <- merge(hex_sf, hex_stats, by = "h3")

# -----------------------------------------------------------------------------
# 4. Universal Plotting Function (Flexible for Binary & TCI)
# -----------------------------------------------------------------------------
theme_map_clean <- theme_void(base_size = 12) +
  theme(
    axis.text = element_text(size = 8, color = "grey30"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.margin = margin(2, 2, 2, 2)
  )

plot_trait_map <- function(data, group_name, var_name, title_text, legend_title, 
                           palette_opt = "magma", scale_limits = c(0, 1), is_percent = TRUE) {
  
  df_sub <- data[data$Group == group_name, ]
  
  p <- ggplot() +
    geom_sf(data = world, fill = "grey92", color = "white", linewidth = 0.15) +
    geom_sf(data = df_sub, aes(fill = .data[[var_name]]), color = NA) +
    coord_sf(xlim = xlim_use, ylim = ylim_use, expand = FALSE) +
    scale_x_continuous(breaks = x_breaks, labels = function(x) paste0(x, "°")) +
    scale_y_continuous(breaks = y_breaks, labels = function(y) paste0(y, "°")) +
    labs(title = title_text) +
    theme_map_clean
  
  if (is_percent) {
    p <- p + scale_fill_viridis_c(
      option = palette_opt, limits = scale_limits, name = legend_title,
      labels = percent_format(accuracy = 1),
      guide = guide_colorbar(barheight = unit(1.5, "in"), barwidth = unit(0.1, "in"))
    )
  } else {
    p <- p + scale_fill_viridis_c(
      option = palette_opt, limits = scale_limits, name = legend_title,
      guide = guide_colorbar(barheight = unit(1.5, "in"), barwidth = unit(0.1, "in"))
    )
  }
  return(p)
}

# -----------------------------------------------------------------------------
# 5. Generate FIGURE 1 (2x2 Matrix: Binary + TCI)
# -----------------------------------------------------------------------------

# Row A: Ornamented Species (%) - Binary result
p1_bin_biv <- plot_trait_map(hex_sf, "Bivalvia", "prop_orn_sp", "Bivalvia", "Ornamented\nSpecies (%)")
p1_bin_brach <- plot_trait_map(hex_sf, "Brachiopoda", "prop_orn_sp", "Brachiopoda", "Ornamented\nSpecies (%)")

# Row B: Mean Trait Complexity (TCI)
# Note: Use a different color scale (e.g., "viridis" or "cividis") to differentiate from Binary
p1_tci_biv <- plot_trait_map(hex_sf, "Bivalvia", "mean_tci_sp", "", "Mean Trait\nComplexity", 
                             palette_opt = "viridis", scale_limits = c(0, 2), is_percent = FALSE)
p1_tci_brach <- plot_trait_map(hex_sf, "Brachiopoda", "mean_tci_sp", "", "Mean Trait\nComplexity", 
                               palette_opt = "viridis", scale_limits = c(0, 2), is_percent = FALSE)

# Combine using patchwork
fig1_main_revised <- (p1_bin_biv | p1_bin_brach) / (p1_tci_biv | p1_tci_brach) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 16, face = "bold"))

ggsave(file.path(OUT_DIR, "Figure_1_Revised_Presence_Complexity.png"), fig1_main_revised, width = 11, height = 9, dpi = 600)

# -----------------------------------------------------------------------------
# 6. Generate FIGURE S1 (Supplementary: Occurrence-Weighted 2x2 Matrix)
# -----------------------------------------------------------------------------

# Row A: Ornamented Occurrences (%)
ps1_bin_biv <- plot_trait_map(hex_sf, "Bivalvia", "prop_orn_occ", "Bivalvia", "Ornamented\nOccurrences (%)")
ps1_bin_brach <- plot_trait_map(hex_sf, "Brachiopoda", "prop_orn_occ", "Brachiopoda", "Ornamented\nOccurrences (%)")

# Row B: Mean Trait Complexity (Occurrence-weighted)
ps1_tci_biv <- plot_trait_map(hex_sf, "Bivalvia", "mean_tci_occ", "", "Mean Trait\nComplexity", 
                              palette_opt = "viridis", scale_limits = c(0, 2), is_percent = FALSE)
ps1_tci_brach <- plot_trait_map(hex_sf, "Brachiopoda", "mean_tci_occ", "", "Mean Trait\nComplexity", 
                                palette_opt = "viridis", scale_limits = c(0, 2), is_percent = FALSE)

figS1_revised <- (ps1_bin_biv | ps1_bin_brach) / (ps1_tci_biv | ps1_tci_brach) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") & 
  theme(plot.tag = element_text(size = 16, face = "bold"))

ggsave(file.path(OUT_DIR, "Fig_S1_Occurrence_Weighted_Full.png"), figS1_revised, width = 11, height = 9, dpi = 600)

# -----------------------------------------------------------------------------
# 7. Morphotype Dominance Maps (Supplementary Figure S2 Combined)
# -----------------------------------------------------------------------------

# --- [A] Define the custom plotting function first ---
# This defines the function so R can find it later
plot_morph <- function(data_sf, var_name, title_text) {
  ggplot() +
    geom_sf(data = world, fill = "grey92", color = "white", linewidth = 0.15) +
    geom_sf(data = data_sf, aes(fill = .data[[var_name]]), color = NA) +
    coord_sf(xlim = xlim_use, ylim = ylim_use, expand = FALSE) +
    scale_x_continuous(breaks = x_breaks, labels = function(x) paste0(x, "°")) +
    scale_y_continuous(breaks = y_breaks, labels = function(y) paste0(y, "°")) +
    scale_fill_manual(values = biv_palette, name = "Dominant Morphotype", drop = FALSE) +
    labs(title = title_text) +
    theme_map_clean # Use the same clean theme as Figure 1
}

# --- [B] Data Preparation for Bivalve Morphotypes ---
biv_labels <- c("a" = "Smooth/Commarginal", 
                "b" = "Radial", 
                "c" = "Concentric", 
                "d" = "Spines/Nodes", 
                "e" = "Cancellate/Others")

tie_priority <- c("d", "e", "b", "c", "a")

get_dominant <- function(x) {
  tb <- table(x)
  winners <- names(tb)[tb == max(tb)]
  winners[order(match(winners, tie_priority))][1]
}

# Use the dt object from previous steps
dt_biv <- dt[Group == "Bivalvia" & Ornament %chin% names(biv_labels)]

# Calculate dominance statistics
dom_occ <- dt_biv[, .(dom_occ = biv_labels[get_dominant(Ornament)]), by = h3]
dom_sp  <- unique(dt_biv[, .(h3, species, Ornament)])[, .(dom_sp = biv_labels[get_dominant(Ornament)]), by = h3]
dom_table <- merge(dom_occ, dom_sp, by = "h3", all = TRUE)

# Factorize for consistent legend ordering
dom_table[, dom_occ := factor(dom_occ, levels = unname(biv_labels[tie_priority]))]
dom_table[, dom_sp  := factor(dom_sp, levels = unname(biv_labels[tie_priority]))]

# --- [C] FIX: Create a proper SF object with H3 column for merging ---
# This resolves the 'by' must specify uniquely valid columns error
hex_geom_sf <- st_as_sf(data.frame(
  h3 = unique(hex_stats$h3), 
  geometry = cell_to_polygon(unique(hex_stats$h3), simple = FALSE)
))

# Final merge for mapping
hex_sf_biv <- merge(hex_geom_sf, dom_table, by = "h3")

# --- [D] Generate the Plot ---
biv_palette <- c("Smooth/Commarginal" = "#999999", "Radial" = "#377EB8", 
                 "Concentric" = "#4DAF4A", "Spines/Nodes" = "#E41A1C", "Cancellate/Others" = "#984EA3")

p_morph_sp  <- plot_morph(hex_sf_biv, "dom_sp", "Species-weighted")
p_morph_occ <- plot_morph(hex_sf_biv, "dom_occ", "Occurrence-weighted")

# Combine using patchwork
figS2_morph_combined <- (p_morph_sp | p_morph_occ) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") & 
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        plot.tag = element_text(size = 16, face = "bold"))

# Save Figure S2
ggsave(file.path(OUT_DIR, "Fig_S2_Morphotype_Combined.png"), figS2_morph_combined, 
       width = 12, height = 6.5, dpi = 600, bg = "white")

cat("\n--- Spatial Mapping Pipeline with TCI Complete ---\n")
