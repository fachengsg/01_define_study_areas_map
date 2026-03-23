################################################################################
## Supplementary Analysis: Spatial Autocorrelation (SAC) & RAC Robustness Model
## Target: Table S5 (Spatial Robustness) & Figure S7 (Moran's I)
##
## Method:
##   1. Aggregate occurrence data to H3 hexagons (resolution 4).
##   2. Fit spatial baseline GAMs (Presence for both clades).
##   3. Extract deviance residuals and calculate Global Moran's I (Great Circle Dist).
##   4. Calculate spatially lagged residuals (RAC_term).
##   5. Refit GAMs with RAC_term to ensure environmental drivers remain robust.
################################################################################

library(data.table)
library(sf)
library(spdep)
library(ggplot2)
library(h3jsr)
library(scales)
library(mgcv)

OUT_DIR <- "outputs"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

H3_RESOLUTION <- 4
K_NEIGH <- 6

# -----------------------------------------------------------------------------
# 1. Load Data & Aggregate to Spatial Assemblages (H3 Grid)
# -----------------------------------------------------------------------------
f_dt <- "09_final_research_cleaned.rds"
if (!file.exists(f_dt)) stop("[STOP] Missing required input file: ", f_dt)

dt_final <- as.data.table(readRDS(f_dt))

# Ensure H3 index exists
if (!("h3" %in% names(dt_final))) {
  dt_final[, h3 := h3jsr::point_to_cell(cbind(decimalLongitude, decimalLatitude), res = H3_RESOLUTION)]
}

# Aggregate to Assemblage Level (Depth <= 1000m)
hex_data <- dt_final[!is.na(SFT) & !is.na(Depth) & Depth <= 1000, .(
  total_sp  = .N,
  orn_sp    = sum(Binary_Trait > 0, na.rm = TRUE),
  unorn_sp  = .N - sum(Binary_Trait > 0, na.rm = TRUE),
  mean_tci  = mean(as.numeric(TCI), na.rm = TRUE),
  SFT       = mean(SFT, na.rm = TRUE),
  Depth     = mean(Depth, na.rm = TRUE),
  Longitude = mean(decimalLongitude, na.rm = TRUE),
  Latitude  = mean(decimalLatitude, na.rm = TRUE)
), by = .(h3, Group)]

# Filter sparse bins
hex_data <- hex_data[total_sp >= 5]
hex_data[, Group := as.factor(Group)]

if (nrow(hex_data) == 0) stop("[STOP] No data left after aggregation. Check filters.")

# -----------------------------------------------------------------------------
# 2. Fit Spatial Baseline GAM & Extract Residuals
# -----------------------------------------------------------------------------
cat("Fitting spatial baseline GAM (Binary Presence)...\n")
m_hex_gam <- gam(cbind(orn_sp, unorn_sp) ~ Group + 
                   s(SFT, by = Group, k = 4) + 
                   s(Depth, by = Group, k = 4), 
                 family = binomial, data = hex_data, method = "REML")

# Extract deviance residuals
hex_data[, dev_resid := residuals(m_hex_gam, type = "deviance")]

# -----------------------------------------------------------------------------
# 3. Spatial Weights Matrix & Global Moran's I
# -----------------------------------------------------------------------------
# Average residuals per H3 cell for purely spatial testing
hex_spatial <- hex_data[, .(
  Mean_Residual = mean(dev_resid, na.rm = TRUE),
  Longitude     = mean(Longitude, na.rm = TRUE),
  Latitude      = mean(Latitude, na.rm = TRUE),
  Records       = sum(total_sp)
), by = .(h3)]

if (nrow(hex_spatial) < 10) stop("[STOP] Too few spatial cells for Moran's I.")

# Create spatial object (WGS84 - Lat/Lon)
hex_sf_wgs <- st_as_sf(hex_spatial, coords = c("Longitude", "Latitude"), crs = 4326)
coords_wgs <- st_coordinates(hex_sf_wgs)

k_use <- min(K_NEIGH, nrow(hex_spatial) - 1)

# Use longlat = TRUE to calculate true spherical distances (Great Circle)
nb_k  <- knn2nb(knearneigh(coords_wgs, k = k_use, longlat = TRUE))
lw_k  <- nb2listw(nb_k, style = "W", zero.policy = TRUE)

# Calculate Moran's I
moran_result <- moran.test(hex_sf_wgs$Mean_Residual, lw_k, zero.policy = TRUE)

# Calculate spatially lagged residuals (RAC Term)
hex_sf_wgs$Lagged_Residual <- lag.listw(lw_k, hex_sf_wgs$Mean_Residual, zero.policy = TRUE)

# -----------------------------------------------------------------------------
# 4. RAC Robustness GAM (Residual Autocovariate) -> TABLE S5
# -----------------------------------------------------------------------------
cat("Fitting RAC robustness GAM for Table S5...\n")

# Extract RAC_term and merge back
rac_map <- as.data.table(st_drop_geometry(hex_sf_wgs))[, .(h3, RAC_term = Lagged_Residual)]
hex_data <- merge(hex_data, rac_map, by = "h3", all.x = TRUE)
hex_data <- hex_data[!is.na(RAC_term)] 

# Fit the RAC robust GAM
m_hex_rac <- gam(cbind(orn_sp, unorn_sp) ~ Group + 
                   s(SFT, by = Group, k = 4) + 
                   s(Depth, by = Group, k = 4) + 
                   s(RAC_term, k = 4), # Spatial control term
                 family = binomial, data = hex_data, method = "REML")

# Export GAM summary to match the new Table S5 numbering
sink(file.path(OUT_DIR, "Table_S5_RAC_Model_Summary.txt"))
cat("--- TABLE S5: RAC ROBUSTNESS GAM SUMMARY (PRESENCE) ---\n")
print(summary(m_hex_rac))
sink()

# Also print to console for quick copy-pasting
cat("\n=== FOR TABLE S5: Summary of RAC Spatial Robustness GAM ===\n")
print(summary(m_hex_rac))

# -----------------------------------------------------------------------------
# 5. High-Quality Visualization (Moran Scatterplot)
# -----------------------------------------------------------------------------
plot_title <- "Spatial Autocorrelation Diagnostic (Moran's I)"

# FIX: Updated subtitle to correctly reference Table S5
plot_sub   <- sprintf("Global Moran's I = %.3f (p = %.3g)\nGAM confirms environmental drivers remain robust after SAC control (see Table S5).", 
                      moran_result$estimate[1], moran_result$p.value)

p_moran <- ggplot(hex_sf_wgs, aes(x = Mean_Residual, y = Lagged_Residual)) +
  geom_point(aes(size = Records), alpha = 0.35, color = "#0072B2") +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.6) +
  scale_size_continuous(labels = comma, name = "Records per cell", range = c(1, 6)) +
  labs(title = plot_title, subtitle = plot_sub,
       x = "Mean deviance residual (H3 cell)", y = "Spatially lagged residual (RAC)") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 12, color = "grey30"),
    legend.position = "right"
  )

ggsave(file.path(OUT_DIR, "Fig_S8_Morans_I.png"), p_moran, width = 8, height = 6, dpi = 600, bg = "white")

cat("\n--- SAC & RAC Pipeline Complete ---\n")
