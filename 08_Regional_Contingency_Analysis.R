## =============================================================================
## 2.7.3–2.7.4 Model validation and robustness (Fig.6 + Fig.Sx; Fig.7 + Table)
##
## 2.7.3 Geographic robustness check (Latitude proxy + depth zones)
##   MAIN:  Orn_bin ~ (Lat_c + Zone) * Group   (Latitude in °N, centered)
##   SUPP:  Orn_bin ~ (Lat_z + Zone) * Group   (Latitude z-scored)
##
## 2.7.4 Regional contingency model (3-level Region)
##   Orn_bin ~ SFT_z + Depth_z + Region * Group   (controls for env; tests residual region signal)
##
## Inputs:
##   - outputs/dt_refined_for_validation.rds  (preferred, produced by Section 2.6)
##   OR
##   - 07_database_with_temp_and_depth.rds
##
## Outputs (OUT_DIR):
##   - Figure6_GeoRobustness_rawLat.png
##   - FigS_GeoRobustness_LatZ.png
##   - Table_Geo_Robustness_Model_rawLat.csv
##   - Table_Geo_Robustness_Model_LatZ.csv
##   - Figure7_Regional_Analysis.png
##   - Table_Regional_Contingency_Model.csv
## =============================================================================

library(data.table)
library(ggplot2)
library(patchwork)
library(scales)

OUT_DIR <- "outputs"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

TARGET_REGIONS <- c("Atl-North","Atl-South","Mediterranean")
TARGET_GROUPS  <- c("Bivalvia","Brachiopoda")

# -----------------------------------------------------------------------------
# Helper 0: load dt_refined robustly
# -----------------------------------------------------------------------------
load_dt_refined <- function(out_dir = "outputs") {
  f_dt <- file.path(out_dir, "dt_refined_for_validation.rds")
  if (file.exists(f_dt)) {
    x <- readRDS(f_dt)
    setDT(x)
    return(x)
  }
  # fallback
  x <- readRDS("07_database_with_temp_and_depth.rds")
  setDT(x)
  x <- x[Region %chin% TARGET_REGIONS & Group %chin% TARGET_GROUPS]
  return(x)
}

# -----------------------------------------------------------------------------
# Helper 1: minimal cleaning shared by all validation steps
# -----------------------------------------------------------------------------
prep_common <- function(x) {
  setDT(x)
  
  # clean Ornament and enforce non-empty
  x[, Ornament := trimws(as.character(Ornament))]
  x <- x[!is.na(Ornament) & nzchar(Ornament)]
  
  # binary response (1=ornamented)
  if (!("Orn_bin" %in% names(x))) {
    x[, Orn_bin := fifelse(Ornament == "a", 0L, 1L)]
  }
  
  # factors
  x[, Group := factor(Group, levels = c("Bivalvia","Brachiopoda"))]
  
  # ensure numeric Depth & Latitude
  x[, Depth_num := as.numeric(Depth)]
  x[, Lat_deg   := as.numeric(decimalLatitude)]
  
  x
}

# -----------------------------------------------------------------------------
# Helper 2: depth zones (0–20, 20–80, 80–200, >200 m)  ✅ UPDATED
# -----------------------------------------------------------------------------
add_depth_zones <- function(x) {
  zone_levels <- c(
    "Coastal (0–20 m)",
    "Inner shelf (20–80 m)",
    "Open shelf (80–200 m)",
    "Beyond shelf (>200 m)"
  )
  
  x[, Zone := fcase(
    Depth_num <= 20,                     zone_levels[1],
    Depth_num >  20 & Depth_num <= 80,   zone_levels[2],
    Depth_num >  80 & Depth_num <= 200,  zone_levels[3],
    Depth_num >  200,                    zone_levels[4],
    default = NA_character_
  )]
  
  x <- x[!is.na(Zone)]
  x[, Zone := factor(Zone, levels = zone_levels, ordered = FALSE)]
  list(dt = x, zone_levels = zone_levels)
}

# -----------------------------------------------------------------------------
# Helper 3: predict with 95% CI (link-scale) — safe for data.table
# -----------------------------------------------------------------------------
predict_with_ci <- function(model, newdt) {
  pred_link <- predict(model, newdata = newdt, type = "link", se.fit = TRUE)
  newdt[, `:=`(fit = as.numeric(pred_link$fit), se = as.numeric(pred_link$se.fit))]
  newdt[, `:=`(
    pred = plogis(fit),
    lo   = plogis(fit - 1.96 * se),
    hi   = plogis(fit + 1.96 * se)
  )]
  newdt
}

# -----------------------------------------------------------------------------
# Load and prepare
# -----------------------------------------------------------------------------
dt_refined <- load_dt_refined(OUT_DIR)
dt_refined <- prep_common(dt_refined)

tmp <- add_depth_zones(dt_refined)
dt_refined <- tmp$dt
zone_levels <- tmp$zone_levels

# Palette and theme consistent with Fig.3
pal <- c("Bivalvia" = "#0072B2", "Brachiopoda" = "#D55E00")

theme_geb <- theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10),
    axis.title = element_text(size = 11),
    axis.text  = element_text(size = 10),
    legend.title = element_blank(),
    panel.grid = element_blank()
  )

# =============================================================================
# 2.7.3 Geographic robustness check (Fig.6 main + Fig.Sx)
# =============================================================================

# Lat variables
lat_mean <- mean(dt_refined$Lat_deg, na.rm = TRUE)
lat_sd   <- sd(dt_refined$Lat_deg, na.rm = TRUE)

dt_refined[, Lat_c := Lat_deg - lat_mean]                  # MAIN (centered degrees)
dt_refined[, Lat_z := (Lat_deg - lat_mean) / lat_sd]       # SUPP (z-score)

# Models
geo_model_raw <- glm(Orn_bin ~ (Lat_c + Zone) * Group, data = dt_refined, family = binomial("logit"))
geo_model_z   <- glm(Orn_bin ~ (Lat_z + Zone) * Group, data = dt_refined, family = binomial("logit"))

# Export coefficient tables
tab_raw <- as.data.table(coef(summary(geo_model_raw)), keep.rownames = "Term")
tab_z   <- as.data.table(coef(summary(geo_model_z)),   keep.rownames = "Term")
fwrite(tab_raw, file.path(OUT_DIR, "Table_Geo_Robustness_Model_rawLat.csv"))
fwrite(tab_z,   file.path(OUT_DIR, "Table_Geo_Robustness_Model_LatZ.csv"))

# Observed bins (2°) for points
dt_refined[, Lat_Bin := round(Lat_deg / 2) * 2]
lat_obs <- dt_refined[, .(prop_orn = mean(Orn_bin), n_records = .N), by = .(Lat_Bin, Group)][n_records >= 50]

# modal Zone for conditional prediction curves
zone_ref <- as.character(dt_refined[, .N, by = Zone][order(-N)][1, Zone])

# Shared Panel B: by-zone observed proportions
zone_obs <- dt_refined[, .(prop_orn = mean(Orn_bin), n_records = .N), by = .(Zone, Group)]

p_zone <- ggplot(zone_obs, aes(x = Zone, y = prop_orn, fill = Group)) +
  geom_col(position = position_dodge(width = 0.8), color = "black", linewidth = 0.25, alpha = 0.9) +
  scale_fill_manual(values = pal) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(title = "Bathymetric corroboration (discrete depth zones)", x = NULL, y = "Ornamented proportion") +
  theme_geb +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        legend.position = "none")

# ---- MAIN Panel A (raw latitude °N) ----
grid_lat_deg <- seq(min(dt_refined$Lat_deg, na.rm = TRUE), max(dt_refined$Lat_deg, na.rm = TRUE), length.out = 200)
new_lat_raw <- CJ(Group = levels(dt_refined$Group), Lat_deg = grid_lat_deg)
new_lat_raw[, Lat_c := Lat_deg - lat_mean]
new_lat_raw[, Zone := factor(zone_ref, levels = zone_levels, ordered = FALSE)]
setDT(new_lat_raw)
new_lat_raw <- predict_with_ci(geo_model_raw, new_lat_raw)

p_lat_raw <- ggplot() +
  geom_ribbon(data = new_lat_raw[Group=="Brachiopoda"], aes(x = Lat_deg, ymin = lo, ymax = hi),
              fill = pal["Brachiopoda"], alpha = 0.18) +
  geom_ribbon(data = new_lat_raw[Group=="Bivalvia"], aes(x = Lat_deg, ymin = lo, ymax = hi),
              fill = pal["Bivalvia"], alpha = 0.18) +
  geom_line(data = new_lat_raw, aes(x = Lat_deg, y = pred, colour = Group, linetype = Group), linewidth = 1.1) +
  geom_point(data = lat_obs, aes(x = Lat_Bin, y = prop_orn, colour = Group, size = n_records), alpha = 0.6) +
  scale_colour_manual(values = pal) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "Latitudinal corroboration (temperature proxy)",
    subtitle = paste0("Lines: GLM predictions at Zone = ", zone_ref, " (±95% CI). Points: binned observations (≥50 records/bin)."),
    x = "Latitude (°N)", y = "Pr(ornamented)", size = "Records/bin"
  ) +
  theme_geb +
  theme(legend.position = "bottom")

fig6_main <- p_lat_raw + p_zone +
  plot_annotation(tag_levels = "A",
                  title = "Fig. 6. Geographic corroboration of environmental drivers",
                  subtitle = "Main text version uses raw latitude (°N) and discrete bathymetric zones.")

ggsave(file.path(OUT_DIR, "Figure6_GeoRobustness_rawLat.png"),
       fig6_main, width = 14, height = 6, dpi = 600, bg = "white")

# ---- SUPP Panel A (Latitude Z) ----
grid_lat_z <- seq(min(dt_refined$Lat_z, na.rm = TRUE), max(dt_refined$Lat_z, na.rm = TRUE), length.out = 200)
new_lat_z <- CJ(Group = levels(dt_refined$Group), Lat_z = grid_lat_z)
new_lat_z[, Zone := factor(zone_ref, levels = zone_levels, ordered = FALSE)]
setDT(new_lat_z)
new_lat_z <- predict_with_ci(geo_model_z, new_lat_z)

lat_obs_z <- copy(lat_obs)
lat_obs_z[, Lat_z := (Lat_Bin - lat_mean) / lat_sd]

p_lat_z <- ggplot() +
  geom_ribbon(data = new_lat_z[Group=="Brachiopoda"], aes(x = Lat_z, ymin = lo, ymax = hi),
              fill = pal["Brachiopoda"], alpha = 0.18) +
  geom_ribbon(data = new_lat_z[Group=="Bivalvia"], aes(x = Lat_z, ymin = lo, ymax = hi),
              fill = pal["Bivalvia"], alpha = 0.18) +
  geom_line(data = new_lat_z, aes(x = Lat_z, y = pred, colour = Group, linetype = Group), linewidth = 1.1) +
  geom_point(data = lat_obs_z, aes(x = Lat_z, y = prop_orn, colour = Group, size = n_records), alpha = 0.6) +
  scale_colour_manual(values = pal) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "Latitudinal corroboration (temperature proxy)",
    subtitle = paste0("Supplementary version (Latitude Z). Zone = ", zone_ref, " (±95% CI)."),
    x = "Latitude (Z)", y = "Pr(ornamented)", size = "Records/bin"
  ) +
  theme_geb +
  theme(legend.position = "bottom")

fig6_supp <- p_lat_z + p_zone +
  plot_annotation(tag_levels = "A",
                  title = "Supplementary Fig. Sx. Geographic corroboration (Latitude Z)",
                  subtitle = "Supplementary version uses z-scored latitude for comparability with standardized predictors.")

ggsave(file.path(OUT_DIR, "FigS_GeoRobustness_LatZ.png"),
       fig6_supp, width = 14, height = 6, dpi = 600, bg = "white")

# =============================================================================
# 2.7.4 Regional contingency analysis (Fig.7 + Table)
# =============================================================================

dt_refined[, Region := factor(Region, levels = c("Atl-North","Atl-South","Mediterranean"))]

if (!("Depth_log" %in% names(dt_refined))) dt_refined[, Depth_log := log10(Depth_num + 1)]
if (!("SFT_z" %in% names(dt_refined)))     dt_refined[, SFT_z := as.numeric(scale(as.numeric(SFT)))]
if (!("Depth_z" %in% names(dt_refined)))   dt_refined[, Depth_z := as.numeric(scale(as.numeric(Depth_log)))]

region_model <- glm(
  Orn_bin ~ SFT_z + Depth_z + Region * Group,
  data = dt_refined,
  family = binomial("logit")
)

fwrite(as.data.table(coef(summary(region_model)), keep.rownames = "Term"),
       file.path(OUT_DIR, "Table_Regional_Contingency_Model.csv"))

region_summary <- dt_refined[, .(prop_orn = mean(Orn_bin), n_records = .N), by = .(Region, Group)][order(Region, Group)]

new_reg <- CJ(Region = levels(dt_refined$Region), Group = levels(dt_refined$Group))
new_reg[, `:=`(SFT_z = 0, Depth_z = 0)]
setDT(new_reg)
new_reg <- predict_with_ci(region_model, new_reg)

p_raw_region <- ggplot(region_summary, aes(x = Region, y = prop_orn, fill = Group)) +
  geom_col(position = position_dodge(width = 0.8), color = "black", linewidth = 0.25, alpha = 0.9) +
  scale_fill_manual(values = pal) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0,1)) +
  labs(title = "Observed ornamentation by region", x = NULL, y = "Ornamented proportion") +
  theme_geb +
  theme(axis.text.x = element_text(angle = 15, hjust = 1),
        legend.position = "bottom")

p_pred_region <- ggplot(new_reg, aes(x = Region, y = pred, colour = Group, group = Group)) +
  geom_point(position = position_dodge(width = 0.5), size = 2.6) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.5), width = 0.15, linewidth = 0.7) +
  scale_colour_manual(values = pal) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0,1)) +
  labs(title = "Predicted ornamentation (controlling for SFT & depth)", x = NULL, y = "Predicted probability") +
  theme_geb +
  theme(axis.text.x = element_text(angle = 15, hjust = 1),
        legend.position = "none")

regional_fig <- p_raw_region + p_pred_region +
  plot_annotation(
    tag_levels = "A",
    title = "Biogeographic contingency in shell ornamentation",
    subtitle = "Three-region model (Mediterranean, Atl-North, Atl-South) controlling for benthic temperature and depth."
  )

ggsave(file.path(OUT_DIR, "Figure7_Regional_Analysis.png"),
       regional_fig, width = 14, height = 6, dpi = 600, bg = "white")
