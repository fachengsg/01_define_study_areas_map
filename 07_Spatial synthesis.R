## =============================================================================
## H3 ornamentation maps (GEB-ready) 
## Input: dt_final (from 07_database_with_temp_and_depth.rds)
## Output:
##   - Figure2_Ornamentation_occ_weighted.png
##   - FigS_species_weighted_ornamentation.png
## =============================================================================

library(data.table)
library(h3jsr)
library(sf)
library(ggplot2)
library(viridis)
library(patchwork)
library(scales)
library(rnaturalearth)

OUT_DIR <- "outputs"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 0) Load data
# -----------------------------------------------------------------------------
if (!exists("dt_final")) {
  dt_final <- readRDS("07_database_with_temp_and_depth.rds")
}
setDT(dt_final)

TARGET_REGIONS <- c("Atl-North", "Atl-South", "Mediterranean")
TARGET_GROUPS  <- c("Bivalvia", "Brachiopoda")

dt <- dt_final[Region %chin% TARGET_REGIONS & Group %chin% TARGET_GROUPS]

# clean Ornament and build binary ornamentation (1 = ornamented)
dt[, Ornament := trimws(as.character(Ornament))]
dt <- dt[!is.na(Ornament) & nzchar(Ornament)]
dt[, Orn_bin := fifelse(Ornament == "a", 0L, 1L)]

# -----------------------------------------------------------------------------
# 1) MANUAL map extent (based on your study domain)
# -----------------------------------------------------------------------------
xlim_use <- c(-32, 38)
ylim_use <- c(8, 70)

x_breaks <- seq(-30, 40, by = 10)
y_breaks <- seq(10, 70, by = 10)

# -----------------------------------------------------------------------------
# 2) H3 indexing + hex statistics
# -----------------------------------------------------------------------------
H3_RESOLUTION <- 4
coords <- cbind(dt$decimalLongitude, dt$decimalLatitude)
dt[, h3 := point_to_cell(coords, res = H3_RESOLUTION)]

# occupancy-weighted proportion (across 0.01° records)
hex_occ <- dt[, .(
  prop_orn_occ = mean(Orn_bin)
), by = .(h3, Group)]

# species-weighted proportion (unique species per hex)
hex_sp <- unique(dt[, .(h3, Group, species, Orn_bin)])
hex_sp_stats <- hex_sp[, .(
  prop_orn_sp = mean(Orn_bin)
), by = .(h3, Group)]

hex_stats <- merge(hex_occ, hex_sp_stats, by = c("h3", "Group"), all.x = TRUE)
fwrite(hex_stats, file.path(OUT_DIR, "hex_stats_H3_res4.csv"))

# -----------------------------------------------------------------------------
# 3) Geometry (Natural Earth + H3 polygons)
# -----------------------------------------------------------------------------
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

hex_ids <- unique(hex_stats$h3)
hex_polys <- cell_to_polygon(hex_ids, simple = FALSE)
hex_geom <- st_as_sf(data.frame(h3 = hex_ids, geometry = hex_polys))
hex_sf <- merge(hex_geom, hex_stats, by = "h3")

# -----------------------------------------------------------------------------
# 4) Panel map function (A  Bivalvia / B  Brachiopoda) + coordinates + right legend
# -----------------------------------------------------------------------------
map_panel_coords <- function(sf_data, grp, var, header_text, legend_title) {
  d <- sf_data[sf_data$Group == grp, ]
  
  # header position inside plot
  hx <- xlim_use[1] + 0.8
  hy <- ylim_use[2] - 0.8
  
  ggplot() +
    geom_sf(data = world, fill = "grey95", color = "grey80", linewidth = 0.2) +
    geom_sf(data = d, aes(fill = .data[[var]]), color = NA) +
    coord_sf(xlim = xlim_use, ylim = ylim_use, expand = FALSE) +
    scale_x_continuous(breaks = x_breaks, labels = function(x) paste0(x, "°")) +
    scale_y_continuous(breaks = y_breaks, labels = function(y) paste0(y, "°")) +
    annotate("text", x = hx, y = hy, label = header_text,
             hjust = 0, vjust = 1, fontface = "bold", size = 4.2) +
    scale_fill_viridis_c(
      option = "magma", limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      labels = scales::percent_format(accuracy = 1),
      name = legend_title,
      guide = guide_colorbar(
        direction = "vertical",
        barheight = grid::unit(90, "pt"),
        barwidth  = grid::unit(10, "pt"),
        ticks = FALSE
      )
    ) +
    theme_void(base_size = 11) +
    theme(
      axis.text = element_text(size = 9, color = "grey20"),
      axis.ticks = element_line(linewidth = 0.3, color = "grey30"),
      axis.ticks.length = grid::unit(2, "pt"),
      axis.title = element_blank(),
      panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.3),
      legend.position = "right",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      plot.margin = margin(5.5, 5.5, 5.5, 5.5)
    )
}

# -----------------------------------------------------------------------------
# 5) MAIN figure: occ-weighted
# -----------------------------------------------------------------------------
p_occ_A <- map_panel_coords(
  hex_sf, "Bivalvia", "prop_orn_occ",
  header_text = "A  Bivalvia",
  legend_title = "Ornamented (%)\nOcc-weighted"
)

p_occ_B <- map_panel_coords(
  hex_sf, "Brachiopoda", "prop_orn_occ",
  header_text = "B  Brachiopoda",
  legend_title = "Ornamented (%)\nOcc-weighted"
) + theme(legend.position = "none")

fig2_main <- (
  (p_occ_A | p_occ_B) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = "Spatial variation in shell ornamentation (H3 DGGS, resolution 4)",
      subtitle = "Occupancy-weighted ornamented proportion across 0.01° grid-cell records."
    )
) & theme(
  plot.title = element_text(face = "bold", size = 13),
  plot.subtitle = element_text(size = 10),
  legend.position = "right"
)

ggsave(file.path(OUT_DIR, "Figure2_Ornamentation_occ_weighted.png"),
       fig2_main, width = 14, height = 6.5, dpi = 600, bg = "white")

# -----------------------------------------------------------------------------
# 6) SUPPLEMENTARY figure: species-weighted
# -----------------------------------------------------------------------------
p_sp_A <- map_panel_coords(
  hex_sf, "Bivalvia", "prop_orn_sp",
  header_text = "A  Bivalvia",
  legend_title = "Ornamented (%)\nSpecies-weighted"
)

p_sp_B <- map_panel_coords(
  hex_sf, "Brachiopoda", "prop_orn_sp",
  header_text = "B  Brachiopoda",
  legend_title = "Ornamented (%)\nSpecies-weighted"
) + theme(legend.position = "none")

fig2_sup <- (
  (p_sp_A | p_sp_B) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = "Supplementary Fig. Sx: Species-weighted ornamentation patterns",
      subtitle = "Each species contributes equally within each hexagon."
    )
) & theme(
  plot.title = element_text(face = "bold", size = 13),
  plot.subtitle = element_text(size = 10),
  legend.position = "right"
)

ggsave(file.path(OUT_DIR, "FigS_species_weighted_ornamentation.png"),
       fig2_sup, width = 14, height = 6.5, dpi = 600, bg = "white")



## =============================================================================
## 23_Spatial_Mapping_of_Bivalve_Morphotypes (GEB-ready, TWO versions)
##
## Outputs (in OUT_DIR):
##  - Figure5_Bivalvia_Dominant_Morphotype_occ_weighted_H3res4.png
##  - FigS_Bivalvia_Dominant_Morphotype_species_weighted_H3res4.png
##
## Definitions:
##  - occ-weighted dominance: mode across 0.01° grid-cell occupancy records within each H3 hex
##  - species-weighted dominance: mode across unique species (each species counts once) within each H3 hex
## =============================================================================

OUT_DIR <- "outputs"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 0) Load data
# -----------------------------------------------------------------------------
if (!exists("dt_final")) {
  dt_final <- readRDS("07_database_with_temp_and_depth.rds")
}
setDT(dt_final)

TARGET_REGIONS <- c("Atl-North", "Atl-South", "Mediterranean")

# Manual map extent (same as latest H3 ornamentation maps)
xlim_use <- c(-32, 38)
ylim_use <- c(8, 70)
x_breaks <- seq(-30, 40, by = 10)
y_breaks <- seq(10, 70, by = 10)

# -----------------------------------------------------------------------------
# 1) Filter to Bivalvia + clean Ornament
# -----------------------------------------------------------------------------
dt_biv <- dt_final[
  Region %chin% TARGET_REGIONS &
    Group == "Bivalvia"
]

dt_biv[, Ornament := trimws(as.character(Ornament))]
dt_biv <- dt_biv[!is.na(Ornament) & nzchar(Ornament)]
dt_biv <- dt_biv[Ornament %chin% c("a","b","c","d","e")]

# -----------------------------------------------------------------------------
# 2) H3 indexing
# -----------------------------------------------------------------------------
H3_RESOLUTION <- 4
coords <- cbind(dt_biv$decimalLongitude, dt_biv$decimalLatitude)
dt_biv[, h3 := point_to_cell(coords, res = H3_RESOLUTION)]

# -----------------------------------------------------------------------------
# 3) Dominant morphotype helpers
# -----------------------------------------------------------------------------
biv_labels <- c(
  "a" = "Smooth / Commarginal",
  "b" = "Radial",
  "c" = "Concentric",
  "d" = "Spines / Nodes / Tubercles",
  "e" = "Cancellate / Reticulate / Others"
)

# Tie-break priority (if equal frequency): pick more "ornamented/complex" first
tie_priority <- c("d", "e", "b", "c", "a")

dominant_mode <- function(x) {
  tb <- table(x)
  maxv <- max(tb)
  winners <- names(tb)[tb == maxv]
  if (length(winners) == 1) return(winners)
  winners[order(match(winners, tie_priority))][1]
}

# -----------------------------------------------------------------------------
# 4A) Occ-weighted dominance (mode across records)
# -----------------------------------------------------------------------------
hex_dom_occ <- dt_biv[, .(
  Dominant_Code = dominant_mode(Ornament),
  n_records = .N,
  n_species = uniqueN(species)
), by = .(h3)]

hex_dom_occ[, Dominant_Morph := biv_labels[Dominant_Code]]
hex_dom_occ[, Dominant_Morph := factor(Dominant_Morph, levels = unname(biv_labels[tie_priority]))]

# -----------------------------------------------------------------------------
# 4B) Species-weighted dominance (mode across unique species)
# -----------------------------------------------------------------------------
# Each species counts once per hex (species-level trait, so one row per species is enough)
hex_species <- unique(dt_biv[, .(h3, species, Ornament)])

hex_dom_sp <- hex_species[, .(
  Dominant_Code = dominant_mode(Ornament),
  n_species = .N
), by = .(h3)]

hex_dom_sp[, Dominant_Morph := biv_labels[Dominant_Code]]
hex_dom_sp[, Dominant_Morph := factor(Dominant_Morph, levels = unname(biv_labels[tie_priority]))]

# -----------------------------------------------------------------------------
# 5) Geometry (Natural Earth + H3 polygons)
# -----------------------------------------------------------------------------
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

make_hex_sf <- function(hex_dt) {
  hex_ids <- unique(hex_dt$h3)
  hex_polys <- cell_to_polygon(hex_ids, simple = FALSE)
  hex_geom <- st_as_sf(data.frame(h3 = hex_ids, geometry = hex_polys))
  merge(hex_geom, hex_dt, by = "h3")
}

hex_sf_occ <- make_hex_sf(hex_dom_occ)
hex_sf_sp  <- make_hex_sf(hex_dom_sp)

# -----------------------------------------------------------------------------
# 6) Plot styling (consistent with your GEB-ready H3 maps)
# -----------------------------------------------------------------------------
biv_palette <- c(
  "Smooth / Commarginal"             = "#999999",
  "Radial"                           = "#377EB8",
  "Concentric"                       = "#4DAF4A",
  "Spines / Nodes / Tubercles"       = "#E41A1C",
  "Cancellate / Reticulate / Others" = "#984EA3"
)

theme_map <- theme_void(base_size = 11) +
  theme(
    axis.text = element_text(size = 9, color = "grey20"),
    axis.ticks = element_line(linewidth = 0.3, color = "grey30"),
    axis.ticks.length = grid::unit(2, "pt"),
    axis.title = element_blank(),
    panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.3),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    plot.title = element_text(face = "bold", size = 13, hjust = 0),
    plot.subtitle = element_text(size = 10, hjust = 0),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  )

plot_morph_map <- function(hex_sf, header, title, subtitle, outfile) {
  hx <- xlim_use[1] + 0.8
  hy <- ylim_use[2] - 0.8
  
  p <- ggplot() +
    geom_sf(data = world, fill = "grey95", color = "grey80", linewidth = 0.2) +
    geom_sf(data = hex_sf, aes(fill = Dominant_Morph), color = NA) +
    coord_sf(xlim = xlim_use, ylim = ylim_use, expand = FALSE) +
    scale_x_continuous(breaks = x_breaks, labels = function(x) paste0(x, "°")) +
    scale_y_continuous(breaks = y_breaks, labels = function(y) paste0(y, "°")) +
    annotate("text", x = hx, y = hy, label = header,
             hjust = 0, vjust = 1, fontface = "bold", size = 4.2) +
    scale_fill_manual(values = biv_palette, name = "Dominant morphotype") +
    labs(title = title, subtitle = subtitle) +
    theme_map
  
  ggsave(file.path(OUT_DIR, outfile), p, width = 14, height = 6.5, dpi = 600, bg = "white")
  p
}

# -----------------------------------------------------------------------------
# 7) Export maps
# -----------------------------------------------------------------------------
p_occ <- plot_morph_map(
  hex_sf_occ,
  header   = "Bivalvia",
  title    = "Geographic dominance of bivalve ornamentation morphotypes",
  subtitle = "Occ-weighted: dominant category per H3 hexagon (resolution 4), based on 0.01° grid-cell occupancies",
  outfile  = "Figure5_Bivalvia_Dominant_Morphotype_occ_weighted_H3res4.png"
)

p_sp <- plot_morph_map(
  hex_sf_sp,
  header   = "Bivalvia",
  title    = "Geographic dominance of bivalve ornamentation morphotypes",
  subtitle = "Species-weighted: dominant category per H3 hexagon (resolution 4), each species counted once per hexagon",
  outfile  = "FigS_Bivalvia_Dominant_Morphotype_species_weighted_H3res4.png"
)

print(p_occ)
print(p_sp)




