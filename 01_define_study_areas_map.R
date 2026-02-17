# ============================================================
# Define study-area polygons and plot an overview map
# Project: SHELLFUTURE
# Author: <Facheng Ye>
# Date: <2026-02>
#
# Notes:
# - Coordinates are in lon/lat (WGS84; EPSG:4326).
# - rnaturalearthhires is OPTIONAL; only required for scale = "large".
# ============================================================

# ----------------------------
# 0) Packages
# ----------------------------
required_pkgs <- c("sf", "ggplot2", "rnaturalearth")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing packages: ", paste(missing_pkgs, collapse = ", "),
    "\nInstall them first, e.g. install.packages(c(", paste0('"', missing_pkgs, '"', collapse = ", "), "))."
  )
}

library(sf)
library(ggplot2)
library(rnaturalearth)

# Optional: high-resolution Natural Earth (needed for scale = "large")
# Recommended for GitHub: keep installation OUTSIDE the script, e.g. in README.
# If you want to use it, install once:
#   install.packages("pak")
#   pak::pak("ropensci/rnaturalearthhires")
has_hires <- requireNamespace("rnaturalearthhires", quietly = TRUE)

# ----------------------------
# 1) Define anchor points (lon, lat)
#    (Do NOT repeat the first point at the end; the helper will close the ring.)
# ----------------------------

# 1. Mediterranean: tightly constrained at the eastern entrance of the Strait of Gibraltar
med_pts <- matrix(c(
  -5.6, 35.5,  # Gibraltar east (south end, Morocco side)
  -5.6, 36.2,  # Gibraltar east (north end, Spain side)
  3.0, 43.5,
  15.0, 46.0,
  36.0, 40.0,
  36.0, 30.0,
  19.0, 30.0,  # includes Marsa al Brega coverage
  3.0, 34.0
), ncol = 2, byrow = TRUE)

# 2. Eastern Atlantic (North)
atl_n_pts <- matrix(c(
  -15.0, 36.0,
  -5.7, 36.0,  # aligned latitude with Atl-South
  -8.5, 42.0,
  -1.0, 43.0,
  1.5, 46.5,  # strengthened Royan coverage
  -1.0, 48.5,
  5.0, 50.0,
  12.0, 58.0,
  12.0, 68.0,
  -15.0, 68.0
), ncol = 2, byrow = TRUE)

# 3. Eastern Atlantic (South): NE corner shifted north to 36.0N to align with Atl-North
atl_s_pts <- matrix(c(
  -20.0, 12.0,
  -5.7, 12.0,
  -5.7, 36.0,  # north anchor aligned with the Gibraltar axis
  -10.0, 37.0,
  -20.0, 37.0
), ncol = 2, byrow = TRUE)

# ----------------------------
# 2) Helper: build sf polygons
# ----------------------------
create_region <- function(pts, name, crs = 4326) {
  stopifnot(is.matrix(pts), ncol(pts) == 2, nrow(pts) >= 3)
  
  # Close the polygon ring explicitly
  pts_closed <- rbind(pts, pts[1, , drop = FALSE])
  
  poly <- st_polygon(list(pts_closed))
  st_sf(Region = name, geometry = st_sfc(poly, crs = crs))
}

study_areas <- rbind(
  create_region(med_pts,   "Mediterranean"),
  create_region(atl_s_pts, "Atl-South"),
  create_region(atl_n_pts, "Atl-North")
)

# ----------------------------
# 3) Basemap (Natural Earth)
# ----------------------------
# Use "large" if hires is installed; otherwise fallback to "medium"
world_plot <- if (has_hires) {
  rnaturalearth::ne_countries(scale = "large", returnclass = "sf")
} else {
  rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
}

# ----------------------------
# 4) Plot
# ----------------------------
p <- ggplot() +
  geom_sf(data = world_plot, fill = "#f2f2f2", color = "#d9d9d9") +
  geom_sf(
    data = study_areas,
    aes(fill = Region),
    alpha = 0.5,
    color = "black",
    linewidth = 0.4
  ) +
  coord_sf(xlim = c(-30, 45), ylim = c(5, 75), expand = FALSE) +
  scale_fill_manual(values = c(
    "Mediterranean" = "#e41a1c",
    "Atl-South"     = "#377eb8",
    "Atl-North"     = "#4daf4a"
  )) +
  theme_minimal() +
  labs(
    title = "SHELLFUTURE study areas (v3 final)",
    fill  = "Region"
  ) +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_line(color = "#ebebeb", linewidth = 0.2)
  )

print(p)

# Optional: save figure (recommended for GitHub repo outputs/)
# ggsave("outputs/study_areas_v3_final.png", p, width = 10, height = 6, dpi = 300)
