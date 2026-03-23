################################################################################
## Macroecological Determinants of Benthic Ornamentation
## Focus: Bivalvia vs. Brachiopoda (Depth <= 1000m)
## Updates: Added Beta Regression (Table S3) & Gaussian TCI GAM (Table S4)
################################################################################

# Load required packages
library(data.table)
library(dplyr)
library(ggplot2)
library(corrplot)  # For correlation matrix visualization
library(car)       # For Variance Inflation Factor (VIF)
library(mgcv)      # For Generalized Additive Models (GAMs)
library(MuMIn)     # For AICc model selection
library(ggeffects) # For extracting marginal effects
library(scales)    # For plot axis formatting
library(patchwork) # For combining multiple plots

# Ensure output directory exists
OUT_DIR <- "outputs"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

################################################################################
## SECTION 1: Data Preparation & Exploration (Correlation & VIF)
################################################################################
cat("\n--- SECTION 1: Data Preparation & Diagnostics ---\n")

# Assuming dt_final is already loaded in your environment
# dt_final <- readRDS("09_final_research_cleaned.rds")

# Format data and filter for depth <= 1000m
dt_analysis <- dt_final %>% 
  filter(Depth <= 1000) %>%
  mutate(
    Group = as.factor(Group),
    Region = as.factor(Region),
    species = as.factor(species),
    Orn_binary = ifelse(Binary_Trait > 0, 1, 0),
    TCI = as.numeric(TCI) 
  )

# 1a. Correlation Matrix
env_vars <- dt_analysis %>% select(SFT, SST, Depth) %>% na.omit()
cor_matrix <- cor(env_vars, method = "spearman")

png(file.path(OUT_DIR, "Fig_S6_Correlation_Matrix.png"), width = 5, height = 5, units = "in", res = 600)
corrplot(cor_matrix, method = "color", type = "upper", 
         addCoef.col = "black", tl.col = "black", diag = FALSE,
         title = "Spearman Correlation of Environmental Gradients", 
         mar = c(0, 0, 2, 0)) 
dev.off() 

# 1b. VIF Diagnostics (For Supplementary Table S6 / S1)
# Fit the full binomial GLM including environmental and categorical predictors
vif_model_full <- glm(Orn_binary ~ SFT + SST + Depth + Group + Region, 
                      family = binomial, data = dt_analysis)

cat("\n=== FOR TABLE S1: Adjusted VIF (Comparable to standard VIF) ===\n")
vif_raw <- vif(vif_model_full)
adjusted_vif <- vif_raw[, "GVIF^(1/(2*Df))"]^2
names(adjusted_vif) <- rownames(vif_raw)
print(round(adjusted_vif, 2))

################################################################################
## SECTION 2: Aggregating to Assemblage Level
################################################################################
cat("\n--- SECTION 2: Aggregating to Assemblage Level ---\n")

dt_agg <- dt_analysis %>%
  mutate(
    SFT_bin = floor(SFT),
    SST_bin = floor(SST),
    Depth_bin = floor(Depth / 50) * 50
  ) %>%
  group_by(Region, Group, SFT_bin, SST_bin, Depth_bin, species) %>%
  summarise(
    Orn_sp = max(Orn_binary, na.rm = TRUE),
    TCI_sp = max(TCI, na.rm = TRUE), 
    .groups = "drop"
  ) %>%
  group_by(Region, Group, SFT_bin, SST_bin, Depth_bin) %>%
  summarise(
    total_sp = n(),
    orn_sp = sum(Orn_sp),
    unorn_sp = total_sp - orn_sp,
    mean_tci = mean(TCI_sp, na.rm = TRUE), 
    SFT_mid = mean(SFT_bin) + 0.5,
    SST_mid = mean(SST_bin) + 0.5,
    Depth_mid = mean(Depth_bin) + 25,
    .groups = "drop"
  ) %>%
  filter(total_sp >= 5) 

################################################################################
## SECTION 3: Model Fitting & Selection (Binomial, Beta, and Gaussian TCI)
################################################################################
cat("\n--- SECTION 3: GAM Fitting & Output Generation ---\n")

# A. Binomial GAMs (Presence - Both Clades)
m_null_bin  <- gam(cbind(orn_sp, unorn_sp) ~ Group + Region, family = binomial, data = dt_agg, method = "ML")
m_sft_bin   <- gam(cbind(orn_sp, unorn_sp) ~ Group + Region + s(SFT_mid, by = Group, k = 4), family = binomial, data = dt_agg, method = "ML")
m_sst_bin   <- gam(cbind(orn_sp, unorn_sp) ~ Group + Region + s(SST_mid, by = Group, k = 4), family = binomial, data = dt_agg, method = "ML")
m_depth_bin <- gam(cbind(orn_sp, unorn_sp) ~ Group + Region + s(Depth_mid, by = Group, k = 4), family = binomial, data = dt_agg, method = "ML")
m_full_bin_ml <- gam(cbind(orn_sp, unorn_sp) ~ Group + Region + s(SFT_mid, by = Group, k = 4) + s(SST_mid, by = Group, k = 4) + s(Depth_mid, by = Group, k = 4), family = binomial, data = dt_agg, method = "ML")

cat("\n=== FOR TABLE 2: AICc Model Selection (Binomial) ===\n")
model_results <- model.sel(m_null_bin, m_sft_bin, m_sst_bin, m_depth_bin, m_full_bin_ml)
print(model_results)

# Refit Best Binomial using REML
m_full_bin <- gam(formula(m_full_bin_ml), family = binomial, data = dt_agg, method = "REML")

# B. Beta Regression GAM (Robustness check for Proportions)
N_bins <- nrow(dt_agg)
dt_beta <- dt_agg %>%
  mutate(
    raw_prop = orn_sp / total_sp,
    # Squeeze transformation strictly into (0, 1) bounds
    prop_beta = (raw_prop * (N_bins - 1) + 0.5) / N_bins
  )

m_beta_full <- gam(prop_beta ~ Group + Region + 
                     s(SFT_mid, by = Group, k = 4) + 
                     s(SST_mid, by = Group, k = 4) + 
                     s(Depth_mid, by = Group, k = 4), 
                   family = betar(link = "logit"), 
                   weights = total_sp, 
                   data = dt_beta, 
                   method = "REML")

# C. Gaussian GAM for TCI (BIVALVIA ONLY)
dt_agg_biv <- dt_agg %>% filter(Group == "Bivalvia")
m_full_tci <- gam(mean_tci ~ Region + s(SFT_mid, k = 4) + s(SST_mid, k = 4) + s(Depth_mid, k = 4), 
                  family = gaussian, weights = total_sp, data = dt_agg_biv, method = "REML")

# --- CONSOLE OUTPUTS FOR TABLES ---
cat("\n=== FOR TABLE 3: Summary of Best Binomial GAM (Presence) ===\n")
print(summary(m_full_bin))

cat("\n=== FOR TABLE S3: Summary of Beta Regression GAM (Robustness) ===\n")
print(summary(m_beta_full))

cat("\n=== FOR TABLE S4: Summary of Gaussian GAM (Complexity - TCI for Bivalvia) ===\n")
print(summary(m_full_tci))

################################################################################
## SECTION 4: Relative Importance (Deviance Partitioning - Binomial Model)
################################################################################
full_dev <- summary(m_full_bin)$dev.expl

drop_group  <- update(m_full_bin, . ~ . - Group - s(SFT_mid, by=Group, k=4) - s(SST_mid, by=Group, k=4) - s(Depth_mid, by=Group, k=4) + s(SFT_mid, k=4) + s(SST_mid, k=4) + s(Depth_mid, k=4))
drop_region <- update(m_full_bin, . ~ . - Region)
drop_sft    <- update(m_full_bin, . ~ . - s(SFT_mid, by=Group, k=4))
drop_sst    <- update(m_full_bin, . ~ . - s(SST_mid, by=Group, k=4))
drop_depth  <- update(m_full_bin, . ~ . - s(Depth_mid, by=Group, k=4))

importance_df <- data.frame(
  Predictor_Raw = c("Clade", "Region", "SFT", "SST", "Depth"),
  Deviance_Lost_Pct = c(
    (full_dev - summary(drop_group)$dev.expl) * 100,
    (full_dev - summary(drop_region)$dev.expl) * 100,
    (full_dev - summary(drop_sft)$dev.expl) * 100,
    (full_dev - summary(drop_sst)$dev.expl) * 100,
    (full_dev - summary(drop_depth)$dev.expl) * 100
  )
) %>% arrange(desc(Deviance_Lost_Pct)) %>% mutate(Predictor = factor(Predictor_Raw, levels = rev(Predictor_Raw)))

p_importance <- ggplot(importance_df, aes(x = Predictor, y = Deviance_Lost_Pct)) +
  geom_segment(aes(x = Predictor, xend = Predictor, y = 0, yend = Deviance_Lost_Pct), color = "grey70", linewidth = 1.5) +
  geom_point(size = 5, color = "#0072B2") +
  geom_text(aes(label = sprintf("%.2f%%", Deviance_Lost_Pct)), hjust = -0.4, size = 4.5) +
  coord_flip() +
  scale_y_continuous(limits = c(0, max(importance_df$Deviance_Lost_Pct) * 1.2), expand = expansion(mult = c(0, 0.1))) +
  labs(x = "", y = "Unique Deviance Explained (%)\n(Loss in explanatory power when term is removed)", title = "Relative Importance of Macroecological Drivers") +
  theme_classic(base_size = 14) +
  theme(axis.text.y = element_text(face = "bold", color = "black"), axis.text.x = element_text(color = "black"), plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "Fig2_Relative_Importance.png"), p_importance, width = 8, height = 4.5, dpi = 600)

################################################################################
## SECTION 5: High-Quality Visualizations
################################################################################
cat("\n--- SECTION 5: Generating Visualizations ---\n")
group_cols <- c("Bivalvia" = "#0072B2", "Brachiopoda" = "#D55E00")
dt_plot_agg <- dt_agg %>% mutate(raw_prop = orn_sp / total_sp)
dt_plot_biv <- dt_plot_agg %>% filter(Group == "Bivalvia") # For TCI plots

# ------------------------------------------------------------------------------
# 5A. FIGURE 2: 2x2 Marginal Effects (Presence & Complexity) 
# ------------------------------------------------------------------------------
pred_sft_bin   <- as.data.frame(ggpredict(m_full_bin, terms = c("SFT_mid [all]", "Group"))) %>% rename(Group = group) 
pred_depth_bin <- as.data.frame(ggpredict(m_full_bin, terms = c("Depth_mid [all]", "Group"))) %>% rename(Group = group) 
pred_sft_tci   <- as.data.frame(ggpredict(m_full_tci, terms = c("SFT_mid [all]"))) %>% mutate(Group = "Bivalvia") 
pred_depth_tci <- as.data.frame(ggpredict(m_full_tci, terms = c("Depth_mid [all]"))) %>% mutate(Group = "Bivalvia") 

custom_size_guide <- guide_legend(override.aes = list(shape = 21, fill = "grey70", color = "black", stroke = 0.5, alpha = 1))

# Row 1: Binary (Presence)
p_sft_bin <- ggplot() +
  geom_point(data = dt_plot_agg, aes(x = SFT_mid, y = raw_prop, size = total_sp, fill = Group), shape = 21, color = "white", alpha = 0.4, stroke = 0.5) +
  geom_ribbon(data = pred_sft_bin, aes(x = x, ymin = conf.low, ymax = conf.high, fill = Group), alpha = 0.25) +
  geom_line(data = pred_sft_bin, aes(x = x, y = predicted, color = Group), linewidth = 1.2) +
  scale_color_manual(values = group_cols) + scale_fill_manual(values = group_cols) +
  scale_size_continuous(name = "Species Richness", range = c(1, 6), guide = custom_size_guide) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0), labels = percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Predicted Proportion", title = "a) Marginal Effect of SFT (Presence)") + 
  theme_classic(base_size = 14) + theme(plot.title = element_text(face = "bold"))

p_depth_bin <- ggplot() +
  geom_point(data = dt_plot_agg, aes(x = Depth_mid, y = raw_prop, size = total_sp, fill = Group), shape = 21, color = "white", alpha = 0.4, stroke = 0.5) +
  geom_ribbon(data = pred_depth_bin, aes(x = x, ymin = conf.low, ymax = conf.high, fill = Group), alpha = 0.25) +
  geom_line(data = pred_depth_bin, aes(x = x, y = predicted, color = Group), linewidth = 1.2) +
  scale_color_manual(values = group_cols) + scale_fill_manual(values = group_cols) +
  scale_size_continuous(name = "Species Richness", range = c(1, 6), guide = custom_size_guide) +
  scale_x_continuous(limits = c(0, 1000), breaks = seq(0, 1000, 200)) + 
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0), labels = percent_format(accuracy = 1)) +
  labs(x = NULL, y = NULL, title = "b) Marginal Effect of Depth (Presence)") + 
  theme_classic(base_size = 14) + theme(plot.title = element_text(face = "bold"))

# Row 2: TCI (Complexity - BIVALVIA ONLY)
p_sft_tci <- ggplot() +
  geom_point(data = dt_plot_biv, aes(x = SFT_mid, y = mean_tci, size = total_sp, fill = Group), shape = 21, color = "white", alpha = 0.4, stroke = 0.5) +
  geom_ribbon(data = pred_sft_tci, aes(x = x, ymin = conf.low, ymax = conf.high, fill = Group), alpha = 0.25) +
  geom_line(data = pred_sft_tci, aes(x = x, y = predicted, color = Group), linewidth = 1.2) +
  scale_color_manual(values = group_cols, drop=FALSE) + scale_fill_manual(values = group_cols, drop=FALSE) +
  scale_size_continuous(name = "Species Richness", range = c(1, 6), guide = custom_size_guide) +
  scale_y_continuous(limits = c(0, max(dt_plot_biv$mean_tci)*1.1), expand = c(0, 0)) +
  labs(x = "Sea Floor Temperature (SFT, °C)", y = "Mean Trait Complexity (TCI)", title = "c) Marginal Effect of SFT (Bivalvia Complexity)") + 
  theme_classic(base_size = 14) + theme(plot.title = element_text(face = "bold"))

p_depth_tci <- ggplot() +
  geom_point(data = dt_plot_biv, aes(x = Depth_mid, y = mean_tci, size = total_sp, fill = Group), shape = 21, color = "white", alpha = 0.4, stroke = 0.5) +
  geom_ribbon(data = pred_depth_tci, aes(x = x, ymin = conf.low, ymax = conf.high, fill = Group), alpha = 0.25) +
  geom_line(data = pred_depth_tci, aes(x = x, y = predicted, color = Group), linewidth = 1.2) +
  scale_color_manual(values = group_cols, drop=FALSE) + scale_fill_manual(values = group_cols, drop=FALSE) +
  scale_size_continuous(name = "Species Richness", range = c(1, 6), guide = custom_size_guide) +
  scale_x_continuous(limits = c(0, 1000), breaks = seq(0, 1000, 200)) + 
  scale_y_continuous(limits = c(0, max(dt_plot_biv$mean_tci)*1.1), expand = c(0, 0)) +
  labs(x = "Depth (m)", y = NULL, title = "d) Marginal Effect of Depth (Bivalvia Complexity)") + 
  theme_classic(base_size = 14) + theme(plot.title = element_text(face = "bold"))

fig3_combined <- (p_sft_bin | p_depth_bin) / (p_sft_tci | p_depth_tci) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.title = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "Figure_3_Combined_Marginal_Effects.png"), fig3_combined, width = 12, height = 10, dpi = 600)


# ------------------------------------------------------------------------------
# 5B. Regional Marginal Breakdown - SPLIT INTO TWO SEPARATE FIGURES
# ------------------------------------------------------------------------------
# 1. Thermal Predictions (SFT)
pred_sft_bin_reg <- as.data.frame(ggpredict(m_full_bin, terms = c("SFT_mid [all]", "Group", "Region"))) %>% 
  rename(Group = group, Region = facet)
pred_sft_tci_reg <- as.data.frame(ggpredict(m_full_tci, terms = c("SFT_mid [all]", "Region"))) %>% 
  mutate(Group = "Bivalvia") %>% rename(Region = group) 

# 2. Bathymetric Predictions (Depth)
pred_depth_bin_reg <- as.data.frame(ggpredict(m_full_bin, terms = c("Depth_mid [all]", "Group", "Region"))) %>% 
  rename(Group = group, Region = facet)
pred_depth_tci_reg <- as.data.frame(ggpredict(m_full_tci, terms = c("Depth_mid [all]", "Region"))) %>% 
  mutate(Group = "Bivalvia") %>% rename(Region = group)

max_tci <- max(dt_plot_biv$mean_tci, na.rm = TRUE) * 1.1

# FIG S8: THERMAL
p_sft_bin_reg <- ggplot() +
  geom_point(data = dt_plot_agg, aes(x = SFT_mid, y = raw_prop, size = total_sp, fill = Group), shape = 21, color = "white", alpha = 0.4, stroke = 0.5) +
  geom_ribbon(data = pred_sft_bin_reg, aes(x = x, ymin = conf.low, ymax = conf.high, fill = Group), alpha = 0.25) +
  geom_line(data = pred_sft_bin_reg, aes(x = x, y = predicted, color = Group), linewidth = 1.2) +
  scale_color_manual(values = group_cols) + scale_fill_manual(values = group_cols) +
  scale_size_continuous(range = c(1, 6), guide = "none") + 
  facet_wrap(~ Region, ncol = 3) + 
  coord_cartesian(ylim = c(0, 1)) + 
  labs(x = NULL, y = "Predicted Proportion", title = "a) Regional Thermal Response (Presence)") +
  theme_bw(base_size = 14) + theme(plot.title = element_text(face="bold"), strip.background = element_rect(fill="grey90"))

p_sft_tci_reg <- ggplot() +
  geom_point(data = dt_plot_biv, aes(x = SFT_mid, y = mean_tci, size = total_sp, fill = Group), shape = 21, color = "white", alpha = 0.4, stroke = 0.5) +
  geom_ribbon(data = pred_sft_tci_reg, aes(x = x, ymin = conf.low, ymax = conf.high, fill = Group), alpha = 0.25) +
  geom_line(data = pred_sft_tci_reg, aes(x = x, y = predicted, color = Group), linewidth = 1.2) +
  scale_color_manual(values = group_cols, drop=FALSE) + scale_fill_manual(values = group_cols, drop=FALSE) +
  scale_size_continuous(name="Species Richness", range = c(1, 6), guide = custom_size_guide) +
  facet_wrap(~ Region, ncol = 3) + 
  coord_cartesian(ylim = c(0, max_tci)) +
  labs(x = "Sea Floor Temperature (SFT, °C)", y = "Mean TCI", title = "b) Regional Thermal Response (Bivalvia Complexity)") +
  theme_bw(base_size = 14) + theme(plot.title = element_text(face="bold"), strip.background = element_blank(), strip.text = element_blank())

fig_s3_thermal <- (p_sft_bin_reg / p_sft_tci_reg) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave(file.path(OUT_DIR, "Fig_S3_Regional_Thermal_Response.png"), fig_s3_thermal, width = 11, height = 8.5, dpi = 600)

# FIG S9: DEPTH
p_depth_bin_reg <- ggplot() +
  geom_point(data = dt_plot_agg, aes(x = Depth_mid, y = raw_prop, size = total_sp, fill = Group), shape = 21, color = "white", alpha = 0.4, stroke = 0.5) +
  geom_ribbon(data = pred_depth_bin_reg, aes(x = x, ymin = conf.low, ymax = conf.high, fill = Group), alpha = 0.25) +
  geom_line(data = pred_depth_bin_reg, aes(x = x, y = predicted, color = Group), linewidth = 1.2) +
  scale_color_manual(values = group_cols) + scale_fill_manual(values = group_cols) +
  scale_size_continuous(range = c(1, 6), guide = "none") + 
  facet_wrap(~ Region, ncol = 3) + 
  scale_x_continuous(breaks = seq(0, 1000, 250)) +
  coord_cartesian(ylim = c(0, 1)) + 
  labs(x = NULL, y = "Predicted Proportion", title = "a) Regional Bathymetric Response (Presence)") +
  theme_bw(base_size = 14) + theme(plot.title = element_text(face="bold"), strip.background = element_rect(fill="grey90"))

p_depth_tci_reg <- ggplot() +
  geom_point(data = dt_plot_biv, aes(x = Depth_mid, y = mean_tci, size = total_sp, fill = Group), shape = 21, color = "white", alpha = 0.4, stroke = 0.5) +
  geom_ribbon(data = pred_depth_tci_reg, aes(x = x, ymin = conf.low, ymax = conf.high, fill = Group), alpha = 0.25) +
  geom_line(data = pred_depth_tci_reg, aes(x = x, y = predicted, color = Group), linewidth = 1.2) +
  scale_color_manual(values = group_cols, drop=FALSE) + scale_fill_manual(values = group_cols, drop=FALSE) +
  scale_size_continuous(name="Species Richness", range = c(1, 6), guide = custom_size_guide) +
  facet_wrap(~ Region, ncol = 3) + 
  scale_x_continuous(breaks = seq(0, 1000, 250)) +
  coord_cartesian(ylim = c(0, max_tci)) +
  labs(x = "Depth (m)", y = "Mean TCI", title = "b) Regional Bathymetric Response (Bivalvia Complexity)") +
  theme_bw(base_size = 14) + theme(plot.title = element_text(face="bold"), strip.background = element_blank(), strip.text = element_blank())

fig_s4_bathymetric <- (p_depth_bin_reg / p_depth_tci_reg) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave(file.path(OUT_DIR, "Fig_S4_Regional_Bathymetric_Response.png"), fig_s4_bathymetric, width = 11, height = 8.5, dpi = 600)

cat("\n--- PIPELINE COMPLETE ---\n")
