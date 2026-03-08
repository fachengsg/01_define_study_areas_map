## =============================================================================
## 2.7.5 Spatial autocorrelation check (Moran's I; H3 res=4; kNN k=6)
##
## Goal:
##   Verify that residual spatial clustering does not compromise inference.
## Method:
##   - Deviance residuals from the convergence (interaction) model (m_int)
##   - Aggregate residuals to H3 hexagons (resolution 4): mean residual per cell
##   - Global Moran's I with centroid-based kNN (k=6) weights in projected space
## Outputs (OUT_DIR):
##   - FigureS1_Morans_I_Diagnostic.png
##   - MoranI_results.txt
## =============================================================================

library(data.table)
library(sf)
library(spdep)
library(ggplot2)
library(h3jsr)
library(scales)

OUT_DIR <- "outputs"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

H3_RESOLUTION <- 4
K_NEIGH <- 6

# -----------------------------------------------------------------------------
# 0) Load dt_refined + convergence model (preferred workflow)
# -----------------------------------------------------------------------------
if (!exists("dt_refined")) {
  f_dt <- file.path(OUT_DIR, "dt_refined_for_validation.rds")
  if (!file.exists(f_dt)) stop("[STOP] Missing: ", f_dt, " (Run Section 2.6 first to generate it).")
  dt_refined <- readRDS(f_dt)
}
setDT(dt_refined)

if (!exists("m_int")) {
  f_m <- file.path(OUT_DIR, "m_int.rds")
  if (!file.exists(f_m)) stop("[STOP] Missing: ", f_m, " (Run Section 2.6 first to generate it).")
  m_int <- readRDS(f_m)
}

need_cols <- c("decimalLongitude","decimalLatitude")
miss <- setdiff(need_cols, names(dt_refined))
if (length(miss) > 0) stop("[STOP] dt_refined missing columns: ", paste(miss, collapse = ", "))

# Work on a copy to avoid polluting dt_refined in memory
dt_tmp <- copy(dt_refined)

# -----------------------------------------------------------------------------
# 1) Deviance residuals from convergence model (interaction model)
# -----------------------------------------------------------------------------
dt_tmp[, dev_resid := residuals(m_int, type = "deviance")]

# -----------------------------------------------------------------------------
# 2) H3 indexing (ensure consistent with maps: res=4)
# -----------------------------------------------------------------------------
if (!("h3" %in% names(dt_tmp))) {
  coords <- cbind(dt_tmp$decimalLongitude, dt_tmp$decimalLatitude)
  dt_tmp[, h3 := h3jsr::point_to_cell(coords, res = H3_RESOLUTION)]
}

# -----------------------------------------------------------------------------
# 3) Aggregate residuals to H3 cells
# -----------------------------------------------------------------------------
hex_residuals <- dt_tmp[, .(
  Mean_Residual = mean(dev_resid, na.rm = TRUE),
  Longitude     = mean(decimalLongitude, na.rm = TRUE),
  Latitude      = mean(decimalLatitude, na.rm = TRUE),
  Records       = .N
), by = .(h3)]

if (nrow(hex_residuals) < 10) stop("[STOP] Too few H3 cells for Moran's I (n<10).")

# -----------------------------------------------------------------------------
# 4) Spatial weights: centroid-based kNN in projected space (EPSG:3857)
# -----------------------------------------------------------------------------
hex_sf_wgs <- st_as_sf(hex_residuals, coords = c("Longitude","Latitude"), crs = 4326)
hex_sf_m   <- st_transform(hex_sf_wgs, 3857)
coords_m   <- st_coordinates(hex_sf_m)

k_use <- min(K_NEIGH, nrow(hex_residuals) - 1)
if (k_use < 1) stop("[STOP] Not enough cells to define neighbors.")

nb_k <- knn2nb(knearneigh(coords_m, k = k_use))
lw_k <- nb2listw(nb_k, style = "W", zero.policy = TRUE)

# -----------------------------------------------------------------------------
# 5) Global Moran's I
# -----------------------------------------------------------------------------
moran_result <- moran.test(hex_sf_wgs$Mean_Residual, lw_k, zero.policy = TRUE)

cat("\n--- GLOBAL MORAN'S I (H3-aggregated deviance residuals) ---\n")
cat("H3 resolution:", H3_RESOLUTION, "\n")
cat("Neighbors: kNN k =", k_use, "(centroid-based; EPSG:3857)\n")
print(moran_result)

# Save reproducibility log
sink(file.path(OUT_DIR, "MoranI_results.txt"))
cat("--- GLOBAL MORAN'S I (H3 res=", H3_RESOLUTION, "; kNN k=", k_use, ") ---\n", sep = "")
print(moran_result)
sink()

# -----------------------------------------------------------------------------
# 6) Moran scatterplot (Supplementary Fig. S1)
# -----------------------------------------------------------------------------
hex_sf_wgs$Lagged_Residual <- lag.listw(lw_k, hex_sf_wgs$Mean_Residual, zero.policy = TRUE)

theme_geb <- theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10),
    axis.title = element_text(size = 11),
    axis.text  = element_text(size = 10),
    panel.grid = element_blank()
  )

p_moran <- ggplot(hex_sf_wgs, aes(x = Mean_Residual, y = Lagged_Residual)) +
  geom_point(aes(size = Records), alpha = 0.35, color = "#0072B2") +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.6) +
  scale_size_continuous(labels = scales::comma) +
  labs(
    title = "Supplementary Fig. S1. Spatial autocorrelation diagnostic",
    subtitle = paste0(
      "Global Moran's I = ", round(moran_result$estimate[1], 3),
      " (p = ", signif(moran_result$p.value, 3), "); H3 res=", H3_RESOLUTION, "; kNN k=", k_use
    ),
    x = "Mean deviance residual (H3 cell)",
    y = "Spatially lagged residual",
    size = "Records per H3 cell"
  ) +
  theme_geb

ggsave(file.path(OUT_DIR, "FigureS1_Morans_I_Diagnostic.png"),
       p_moran, width = 8, height = 6, dpi = 600, bg = "white")
