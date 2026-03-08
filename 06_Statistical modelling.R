## =============================================================================
## 2.6 Statistical modelling and macroecological inference
## MAIN FIGURE 3 = raw predictors (SST, SFT, Depth[m] via log10(Depth+1))
## SUPPLEMENTARY = Z-scored predictors (SST_z, SFT_z, Depth_z)
##
## Outputs (OUT_DIR):
##  MAIN:
##   - Table3_Driver_Analysis_raw.csv
##   - Table4_Convergence_Analysis_raw.csv
##   - Figure3A_EffectSizes_raw.png
##   - Figure3BC_Responses_raw.png
##   - Pred_SFT_raw_withCI.csv
##   - Pred_Depth_raw_withCI.csv
##  SUPP:
##   - TableS3_Driver_Analysis_z.csv
##   - TableS4_Convergence_Analysis_z.csv
##   - FigS3A_EffectSizes_z.png
##   - FigS3BC_Responses_z.png
##   - Pred_SFT_z_withCI.csv
##   - Pred_Depth_z_withCI.csv
##
## Also saved for Section 2.7 scripts:
##   - dt_refined_for_validation.rds
##   - m_add.rds, m_int.rds   (MAIN models, for Moran's I etc.)
## =============================================================================

library(data.table)
library(ggplot2)
library(patchwork)

OUT_DIR <- "outputs"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 0) Load data
# -----------------------------------------------------------------------------
if (!exists("dt_final")) {
  dt_final <- readRDS("09_final_research_cleaned.rds")
}
setDT(dt_final)

TARGET_REGIONS <- c("Atl-North", "Atl-South", "Mediterranean")
TARGET_GROUPS  <- c("Bivalvia", "Brachiopoda")

dt <- dt_final[Region %chin% TARGET_REGIONS & Group %chin% TARGET_GROUPS]

# Clean Ornament
dt[, Ornament := trimws(as.character(Ornament))]
dt <- dt[!is.na(Ornament) & nzchar(Ornament)]

# Response: binary ornamentation (1=ornamented)
dt[, Orn_bin := fifelse(Ornament == "a", 0L, 1L)]

# Predictor presence
need_cols <- c("SST", "SFT", "Depth")
miss <- setdiff(need_cols, names(dt))
if (length(miss) > 0) stop("[STOP] Missing columns in dt_final: ", paste(miss, collapse = ", "))

dt <- dt[!is.na(SST) & !is.na(SFT) & !is.na(Depth)]

# Depth transform used in modelling (but plots can show Depth in meters)
dt[, Depth_m := as.numeric(Depth)]
dt[, Depth_log := log10(Depth_m + 1)]

# Group factor (Bivalvia reference)
dt[, Group := factor(Group, levels = c("Bivalvia", "Brachiopoda"))]

# Z-score versions for Supplementary
dt[, SST_z   := as.numeric(scale(as.numeric(SST)))]
dt[, SFT_z   := as.numeric(scale(as.numeric(SFT)))]
dt[, Depth_z := as.numeric(scale(as.numeric(Depth_log)))]

# Save modelling dataset for Section 2.7 scripts
dt_refined <- copy(dt)
saveRDS(dt_refined, file.path(OUT_DIR, "dt_refined_for_validation.rds"))

# -----------------------------------------------------------------------------
# Helper: export coefficient table
# -----------------------------------------------------------------------------
export_coef_table <- function(mod, out_csv) {
  tab <- as.data.table(coef(summary(mod)), keep.rownames = "Term")
  fwrite(tab, file.path(OUT_DIR, out_csv))
  invisible(tab)
}

# Helper: predictions with 95% CI (on link scale)
predict_with_ci <- function(mod, newdt) {
  pred_link <- predict(mod, newdata = newdt, type = "link", se.fit = TRUE)
  newdt[, `:=`(fit = as.numeric(pred_link$fit), se = as.numeric(pred_link$se.fit))]
  newdt[, `:=`(
    pred = plogis(fit),
    lo   = plogis(fit - 1.96 * se),
    hi   = plogis(fit + 1.96 * se)
  )]
  newdt
}

# Plot theme (consistent with your Fig.3 style)
pal <- c("Bivalvia" = "#0072B2", "Brachiopoda" = "#D55E00")

theme_geb <- theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.title = element_text(size = 11),
    axis.text  = element_text(size = 10),
    legend.title = element_blank(),
    legend.position = "bottom",
    panel.grid = element_blank()
  )

tag_style <- theme(
  plot.tag = element_text(face = "bold", size = 12),
  plot.tag.position = c(0.02, 0.98)
)

# =============================================================================
# MAIN MODELS (raw predictors)
# =============================================================================

# -----------------------------------------------------------------------------
# 1) MAIN Model 1: additive (raw)
# -----------------------------------------------------------------------------
m_add_raw <- glm(
  Orn_bin ~ SST + SFT + Depth_log + Group,
  data = dt,
  family = binomial(link = "logit")
)

# -----------------------------------------------------------------------------
# 2) MAIN Model 2: interaction / convergence (raw)
# -----------------------------------------------------------------------------
m_int_raw <- glm(
  Orn_bin ~ (SFT + Depth_log) * Group,
  data = dt,
  family = binomial(link = "logit")
)

# Export MAIN tables
export_coef_table(m_add_raw, "Table3_Driver_Analysis_raw.csv")
export_coef_table(m_int_raw, "Table4_Convergence_Analysis_raw.csv")

# Save MAIN models for Section 2.7 (Moran’s I etc.)
m_add <- m_add_raw
m_int <- m_int_raw
saveRDS(m_add, file.path(OUT_DIR, "m_add.rds"))
saveRDS(m_int, file.path(OUT_DIR, "m_int.rds"))

# -----------------------------------------------------------------------------
# 3) MAIN Figure 3A: unstandardized coefficients (raw units)
# -----------------------------------------------------------------------------
coef_dt <- as.data.table(coef(summary(m_add_raw)), keep.rownames = "Term")
keep_terms <- c("SST", "SFT", "Depth_log", "GroupBrachiopoda")
coef_dt <- coef_dt[Term %chin% keep_terms]

coef_dt[, `:=`(
  lo = Estimate - 1.96 * `Std. Error`,
  hi = Estimate + 1.96 * `Std. Error`
)]

label_map <- c(
  "SST" = "SST (°C)",
  "SFT" = "SFT (°C)",
  "Depth_log" = "log10(Depth+1)",
  "GroupBrachiopoda" = "Brachiopoda\n(vs Bivalvia)"
)
coef_dt[, Term_lab := factor(label_map[Term], levels = rev(label_map[keep_terms]))]

xr <- range(c(coef_dt$lo, coef_dt$hi), na.rm = TRUE)
pad <- 0.03 * diff(xr)

pA_raw <- ggplot(coef_dt, aes(x = Estimate, y = Term_lab)) +
  geom_vline(xintercept = 0, linewidth = 0.5, color = "grey35") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.18, linewidth = 0.75, color = "grey30") +
  geom_point(size = 2.4, color = "black") +
  coord_cartesian(xlim = c(xr[1] - pad, xr[2] + pad)) +
  labs(
    title = "Effect sizes (additive model; raw units)",
    x = "Log-odds coefficient (β, 95% CI)",
    y = NULL,
    tag = "A"
  ) +
  theme_geb + tag_style +
  theme(legend.position = "none")

ggsave(file.path(OUT_DIR, "Figure3A_EffectSizes_raw.png"),
       pA_raw, width = 6, height = 4, dpi = 600, bg = "white")

# -----------------------------------------------------------------------------
# 4) MAIN Figure 3B/C: response curves with 95% CI (raw axes)
#    - B: SFT (°C), hold Depth_log at mean
#    - C: Depth (m) on x-axis, model uses Depth_log, hold SFT at mean
# -----------------------------------------------------------------------------
# SFT grid (°C)
grid_sft <- seq(min(dt$SFT), max(dt$SFT), length.out = 250)
new_sft <- CJ(Group = levels(dt$Group), SFT = grid_sft)
new_sft[, Depth_log := mean(dt$Depth_log)]
new_sft <- predict_with_ci(m_int_raw, new_sft)
pred_sft_raw <- new_sft[, .(Group, SFT, pred, lo, hi)]
fwrite(pred_sft_raw, file.path(OUT_DIR, "Pred_SFT_raw_withCI.csv"))

# Depth grid (meters) -> compute Depth_log for prediction
# Use a log-spaced grid for nicer resolution across depth
dmin <- max(min(dt$Depth_m), 0)
dmax <- max(dt$Depth_m)
grid_depth_log <- seq(log10(dmin + 1), log10(dmax + 1), length.out = 250)
grid_depth_m <- 10^grid_depth_log - 1

new_dep <- CJ(Group = levels(dt$Group), Depth_m = grid_depth_m)
new_dep[, Depth_log := log10(Depth_m + 1)]
new_dep[, SFT := mean(dt$SFT)]
new_dep <- predict_with_ci(m_int_raw, new_dep)
pred_dep_raw <- new_dep[, .(Group, Depth_m, pred, lo, hi)]
fwrite(pred_dep_raw, file.path(OUT_DIR, "Pred_Depth_raw_withCI.csv"))

plot_resp_raw <- function(dat, xvar, xlabel, title_txt, tag_txt, x_log10 = FALSE) {
  dat <- copy(dat)
  dat[, Group := factor(Group, levels = c("Brachiopoda", "Bivalvia"))]  # draw Brachiopoda first
  
  p <- ggplot() +
    geom_ribbon(data = dat[Group=="Brachiopoda"],
                aes(x = get(xvar), ymin = lo, ymax = hi),
                fill = pal["Brachiopoda"], alpha = 0.22) +
    geom_ribbon(data = dat[Group=="Bivalvia"],
                aes(x = get(xvar), ymin = lo, ymax = hi),
                fill = pal["Bivalvia"], alpha = 0.22) +
    # CI boundary lines (ensure visibility)
    geom_line(data = dat, aes(x = get(xvar), y = lo, colour = Group), linewidth = 0.5, alpha = 0.9) +
    geom_line(data = dat, aes(x = get(xvar), y = hi, colour = Group), linewidth = 0.5, alpha = 0.9) +
    geom_line(data = dat, aes(x = get(xvar), y = pred, colour = Group, linetype = Group), linewidth = 1.2) +
    scale_colour_manual(values = pal) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(title = title_txt, x = xlabel, y = "Pr(ornamented)", tag = tag_txt) +
    theme_geb + tag_style
  
  if (x_log10) {
    p <- p + scale_x_continuous(trans = "log10", labels = scales::comma)
  }
  p
}

pB_raw <- plot_resp_raw(pred_sft_raw, "SFT", "SFT (°C)",
                        "Response to benthic temperature (SFT)", "B") +
  theme(legend.position = "none")

pC_raw <- plot_resp_raw(pred_dep_raw, "Depth_m", "Depth (m)",
                        "Response to depth", "C", x_log10 = TRUE)

fig_BC_raw <- (pB_raw | pC_raw) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(OUT_DIR, "Figure3BC_Responses_raw.png"),
       fig_BC_raw, width = 12, height = 4.8, dpi = 600, bg = "white")

# =============================================================================
# SUPPLEMENTARY MODELS (Z-scored predictors)
# =============================================================================

m_add_z <- glm(
  Orn_bin ~ SST_z + SFT_z + Depth_z + Group,
  data = dt,
  family = binomial(link = "logit")
)

m_int_z <- glm(
  Orn_bin ~ (SFT_z + Depth_z) * Group,
  data = dt,
  family = binomial(link = "logit")
)

export_coef_table(m_add_z, "TableS3_Driver_Analysis_z.csv")
export_coef_table(m_int_z, "TableS4_Convergence_Analysis_z.csv")

# FigS3A: standardized effect sizes (Z)
coef_z <- as.data.table(coef(summary(m_add_z)), keep.rownames = "Term")
keep_terms_z <- c("SST_z", "SFT_z", "Depth_z", "GroupBrachiopoda")
coef_z <- coef_z[Term %chin% keep_terms_z]
coef_z[, `:=`(lo = Estimate - 1.96 * `Std. Error`, hi = Estimate + 1.96 * `Std. Error`)]

label_map_z <- c(
  "SST_z" = "SST (Z)",
  "SFT_z" = "SFT (Z)",
  "Depth_z" = "log10(Depth+1) (Z)",
  "GroupBrachiopoda" = "Brachiopoda\n(vs Bivalvia)"
)
coef_z[, Term_lab := factor(label_map_z[Term], levels = rev(label_map_z[keep_terms_z]))]

xr2 <- range(c(coef_z$lo, coef_z$hi), na.rm = TRUE)
pad2 <- 0.03 * diff(xr2)

pA_z <- ggplot(coef_z, aes(x = Estimate, y = Term_lab)) +
  geom_vline(xintercept = 0, linewidth = 0.5, color = "grey35") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.18, linewidth = 0.75, color = "grey30") +
  geom_point(size = 2.4, color = "black") +
  coord_cartesian(xlim = c(xr2[1] - pad2, xr2[2] + pad2)) +
  labs(
    title = "Standardized effect sizes (additive model; Z)",
    x = "Log-odds coefficient (β, 95% CI)",
    y = NULL,
    tag = "A"
  ) +
  theme_geb + tag_style +
  theme(legend.position = "none")

ggsave(file.path(OUT_DIR, "FigS3A_EffectSizes_z.png"),
       pA_z, width = 6, height = 4, dpi = 600, bg = "white")

# FigS3BC: predictions on Z scale
grid_sft_z <- seq(min(dt$SFT_z), max(dt$SFT_z), length.out = 250)
new_sft_z <- CJ(Group = levels(dt$Group), SFT_z = grid_sft_z)
new_sft_z[, Depth_z := 0]
new_sft_z <- predict_with_ci(m_int_z, new_sft_z)
pred_sft_z <- new_sft_z[, .(Group, SFT_z, pred, lo, hi)]
fwrite(pred_sft_z, file.path(OUT_DIR, "Pred_SFT_z_withCI.csv"))

grid_dep_z <- seq(min(dt$Depth_z), max(dt$Depth_z), length.out = 250)
new_dep_z <- CJ(Group = levels(dt$Group), Depth_z = grid_dep_z)
new_dep_z[, SFT_z := 0]
new_dep_z <- predict_with_ci(m_int_z, new_dep_z)
pred_dep_z <- new_dep_z[, .(Group, Depth_z, pred, lo, hi)]
fwrite(pred_dep_z, file.path(OUT_DIR, "Pred_Depth_z_withCI.csv"))

pB_z <- plot_resp_raw(pred_sft_z, "SFT_z", "SFT (Z)",
                      "Response to benthic temperature (SFT)", "B") +
  theme(legend.position = "none")

pC_z <- plot_resp_raw(pred_dep_z, "Depth_z", "log10(Depth+1) (Z)",
                      "Response to depth", "C")

fig_BC_z <- (pB_z | pC_z) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(OUT_DIR, "FigS3BC_Responses_z.png"),
       fig_BC_z, width = 12, height = 4.8, dpi = 600, bg = "white")
