################################################################################
## Spatial Analysis: Shell Ornamentation Mapping (H3 Hexagonal Grids)
## 
## Output: 
## - Figure_1_Main_Species_Weighted.png (Main Text)
## - Fig_S1_Occurrence_Weighted.png     (Supplementary)
## - Fig_S2_Morphotype_Species.png      (Supplementary)
## - Fig_S3_Morphotype_Occurrence.png   (Supplementary)
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
# Load data
if (!exists("dt_final")) {
  dt_final <- readRDS("09_final_research_cleaned.rds")
}
setDT(dt_final)

# Filter domain & clean Ornamentation
dt <- dt_final[Region %chin% TARGET_REGIONS & Group %chin% TARGET_GROUPS]
dt[, Ornament := trimws(as.character(Ornament))]
dt <- dt[!is.na(Ornament) & nzchar(Ornament)]

# Binary ornamentation (1=ornamented, 0=smooth/commarginal)
dt[, Orn_bin := fifelse(Ornament == "a", 0L, 1L)]

# Compute H3 Index
if (!("h3" %in% names(dt))) {
  dt[, h3 := point_to_cell(cbind(decimalLongitude, decimalLatitude), res = H3_RESOLUTION)]
}

# -----------------------------------------------------------------------------
# 3. Calculate Hexagon Statistics (Species & Occurrence Weighted)
# -----------------------------------------------------------------------------
# Species-weighted (Main Text focus - removes sampling bias)
hex_sp <- unique(dt[, .(h3, Group, species, Orn_bin)])
hex_stats_sp <- hex_sp[, .(prop_orn_sp = mean(Orn_bin), n_species = .N), by = .(h3, Group)]

# Occurrence-weighted (Supplementary)
hex_stats_occ <- dt[, .(prop_orn_occ = mean(Orn_bin), n_records = .N), by = .(h3, Group)]

# Merge stats and generate spatial SF object
hex_stats <- merge(hex_stats_sp, hex_stats_occ, by = c("h3", "Group"), all = TRUE)
fwrite(hex_stats, file.path(OUT_DIR, "hex_stats_H3_res4.csv"))

hex_polys <- cell_to_polygon(unique(hex_stats$h3), simple = FALSE)
hex_geom <- st_as_sf(data.frame(h3 = unique(hex_stats$h3), geometry = hex_polys))
hex_sf <- merge(hex_geom, hex_stats, by = "h3")

# -----------------------------------------------------------------------------
# 4. Universal Plotting Theme & Function
# -----------------------------------------------------------------------------
theme_map_geb <- theme_void(base_size = 12) +
  theme(
    axis.text = element_text(size = 9, color = "grey20"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.margin = margin(5, 5, 5, 5)
  )

plot_h3_map <- function(data, group_name, var_name, title_text, legend_title) {
  df_sub <- data[data$Group == group_name, ]
  
  ggplot() +
    geom_sf(data = world, fill = "grey90", color = "white", linewidth = 0.2) +
    geom_sf(data = df_sub, aes(fill = .data[[var_name]]), color = NA) +
    coord_sf(xlim = xlim_use, ylim = ylim_use, expand = FALSE) +
    scale_x_continuous(breaks = x_breaks, labels = function(x) paste0(x, "°")) +
    scale_y_continuous(breaks = y_breaks, labels = function(y) paste0(y, "°")) +
    scale_fill_viridis_c(
      option = "magma", limits = c(0, 1), breaks = seq(0, 1, 0.25),
      labels = percent_format(accuracy = 1), name = legend_title,
      guide = guide_colorbar(barheight = unit(3, "in"), barwidth = unit(0.15, "in"), ticks = FALSE)
    ) +
    labs(title = title_text) +
    theme_map_geb
}

# -----------------------------------------------------------------------------
# 5. Generate FIGURE 1 (Main Text: Species-Weighted)
# -----------------------------------------------------------------------------
p1a <- plot_h3_map(hex_sf, "Bivalvia", "prop_orn_sp", "Bivalvia", "Ornamented\nSpecies (%)")
p1b <- plot_h3_map(hex_sf, "Brachiopoda", "prop_orn_sp", "Brachiopoda", "Ornamented\nSpecies (%)")

fig1_main <- (p1a | p1b) + 
  plot_layout(guides = "collect") + 
  plot_annotation(tag_levels = "a") & 
  theme(plot.tag = element_text(size = 16, face = "bold"))

ggsave(file.path(OUT_DIR, "Figure_1_Main_Species_Weighted.png"), fig1_main, width = 12, height = 6, dpi = 600, bg = "white")

# -----------------------------------------------------------------------------
# 6. Generate FIGURE S1 (Supplementary: Occurrence-Weighted)
# -----------------------------------------------------------------------------
ps1a <- plot_h3_map(hex_sf, "Bivalvia", "prop_orn_occ", "Bivalvia", "Ornamented\nOccurrences (%)")
ps1b <- plot_h3_map(hex_sf, "Brachiopoda", "prop_orn_occ", "Brachiopoda", "Ornamented\nOccurrences (%)")

figS1_supp <- (ps1a | ps1b) + 
  plot_layout(guides = "collect") + 
  plot_annotation(tag_levels = "a") & 
  theme(plot.tag = element_text(size = 16, face = "bold"))

ggsave(file.path(OUT_DIR, "Fig_S1_Occurrence_Weighted.png"), figS1_supp, width = 12, height = 6, dpi = 600, bg = "white")

# -----------------------------------------------------------------------------
# 7. Morphotype Dominance Maps (Supplementary Figure S2 Combined)
# -----------------------------------------------------------------------------
biv_labels <- c("a" = "Smooth/Commarginal", "b" = "Radial", "c" = "Concentric", 
                "d" = "Spines/Nodes", "e" = "Cancellate/Others")
tie_priority <- c("d", "e", "b", "c", "a")

get_dominant <- function(x) {
  tb <- table(x)
  winners <- names(tb)[tb == max(tb)]
  winners[order(match(winners, tie_priority))][1]
}

dt_biv <- dt[Group == "Bivalvia" & Ornament %chin% names(biv_labels)]

# Calculate dominance
dom_occ <- dt_biv[, .(dom_occ = biv_labels[get_dominant(Ornament)]), by = h3]
dom_sp  <- unique(dt_biv[, .(h3, species, Ornament)])[, .(dom_sp = biv_labels[get_dominant(Ornament)]), by = h3]
dom_table <- merge(dom_occ, dom_sp, by = "h3", all = TRUE)

dom_table[, dom_occ := factor(dom_occ, levels = unname(biv_labels[tie_priority]))]
dom_table[, dom_sp  := factor(dom_sp, levels = unname(biv_labels[tie_priority]))]
hex_sf_biv <- merge(hex_geom, dom_table, by = "h3")

# Plotting Morphotypes
biv_palette <- c("Smooth/Commarginal" = "#999999", "Radial" = "#377EB8", 
                 "Concentric" = "#4DAF4A", "Spines/Nodes" = "#E41A1C", "Cancellate/Others" = "#984EA3")

plot_morph <- function(data_sf, var_name, title_text) {
  ggplot() +
    geom_sf(data = world, fill = "grey90", color = "white", linewidth = 0.2) +
    geom_sf(data = data_sf, aes(fill = .data[[var_name]]), color = NA) +
    coord_sf(xlim = xlim_use, ylim = ylim_use, expand = FALSE) +
    scale_x_continuous(breaks = x_breaks, labels = function(x) paste0(x, "°")) +
    scale_y_continuous(breaks = y_breaks, labels = function(y) paste0(y, "°")) +
    scale_fill_manual(values = biv_palette, name = "Dominant Morphotype", drop = FALSE) +
    labs(title = title_text) +
    theme_map_geb
}

# Generate individual panels with simplified titles
p_morph_sp  <- plot_morph(hex_sf_biv, "dom_sp", "Species-weighted")
p_morph_occ <- plot_morph(hex_sf_biv, "dom_occ", "Occurrence-weighted")

# Combine using patchwork, collect legends to the bottom
figS2_morph_combined <- (p_morph_sp | p_morph_occ) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") & 
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    plot.tag = element_text(size = 16, face = "bold")
  )

# Save the combined figure
ggsave(file.path(OUT_DIR, "Fig_S2_Morphotype_Combined.png"), figS2_morph_combined, 
       width = 12, height = 6.5, dpi = 600, bg = "white")

cat("\n--- Spatial Mapping Pipeline Complete ---\n")
