################################################################################
## Macroecological Determinants of Benthic Ornamentation
## Focus: Bivalvia vs. Brachiopoda (Depth <= 1000m)
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

################################################################################
## SECTION 1: Data Preparation & Exploration (Correlation & VIF)
################################################################################
cat("\n--- SECTION 1: Data Preparation & Correlation Analysis ---\n")

# Assuming dt_final is already loaded in your environment
# dt_final <- readRDS("09_final_research_cleaned.rds")

# Format data and filter for <= 1000m
dt_analysis <- dt_final %>% 
  filter(Depth <= 1000) %>%
  mutate(
    Group = as.factor(Group),
    Region = as.factor(Region),
    species = as.factor(species),
    Orn_binary = ifelse(Orn_bin > 0, 1, 0)
  )

# 1a. Correlation Matrix among continuous predictors
env_vars <- dt_analysis %>% select(SFT, SST, Depth) %>% na.omit()
cor_matrix <- cor(env_vars, method = "spearman")

# Save correlation plot for Supplementary Material (Figure S3)
# Set to 600 DPI and 5x5 inches for publication-quality rendering
png("Fig_S3_Correlation_Matrix.png", width = 5, height = 5, units = "in", res = 600)

corrplot(cor_matrix, method = "color", type = "upper", 
         addCoef.col = "black", tl.col = "black", diag = FALSE,
         title = "Spearman Correlation of Environmental Gradients", 
         mar = c(0, 0, 2, 0)) # Increased top margin slightly to prevent title cutoff

dev.off() # Close the graphical device to actually save the file

# 1b. Variance Inflation Factor (VIF) Check
vif_model <- glm(Orn_binary ~ SFT + SST + Depth + Group + Region, 
                 family = binomial, data = dt_analysis)
cat("\nVariance Inflation Factors (VIF):\n")
print(vif(vif_model))

################################################################################
## SECTION 2: Aggregating to Assemblage Level (Mitigating bias)
################################################################################
cat("\n--- SECTION 2: Aggregating to Assemblage Level ---\n")

# Bin data (1°C for temps, 50m for depth)
dt_agg <- dt_analysis %>%
  mutate(
    SFT_bin = floor(SFT),
    SST_bin = floor(SST),
    Depth_bin = floor(Depth / 50) * 50
  ) %>%
  # Step 1: Species-level representation within each environmental bin
  group_by(Region, Group, SFT_bin, SST_bin, Depth_bin, species) %>%
  summarise(Orn_sp = max(Orn_binary, na.rm = TRUE), .groups = "drop") %>%
  # Step 2: Assemblage-level proportion
  group_by(Region, Group, SFT_bin, SST_bin, Depth_bin) %>%
  summarise(
    total_sp = n(),
    orn_sp = sum(Orn_sp),
    unorn_sp = total_sp - orn_sp,
    SFT_mid = mean(SFT_bin) + 0.5,
    SST_mid = mean(SST_bin) + 0.5,
    Depth_mid = mean(Depth_bin) + 25,
    .groups = "drop"
  ) %>%
  filter(total_sp >= 5) # Remove sparse bins to reduce statistical noise

################################################################################
## SECTION 3: Main Binomial GAM & AICc Model Selection
################################################################################
cat("\n--- SECTION 3: Model Selection (AICc) ---\n")

# Note: method = "ML" is required for valid AIC comparison of fixed effects
m_null  <- gam(cbind(orn_sp, unorn_sp) ~ Group + Region, 
               family = binomial, data = dt_agg, method = "ML")

m_sft   <- gam(cbind(orn_sp, unorn_sp) ~ Group + Region + s(SFT_mid, by = Group, k = 4), 
               family = binomial, data = dt_agg, method = "ML")

m_sst   <- gam(cbind(orn_sp, unorn_sp) ~ Group + Region + s(SST_mid, by = Group, k = 4), 
               family = binomial, data = dt_agg, method = "ML")

m_depth <- gam(cbind(orn_sp, unorn_sp) ~ Group + Region + s(Depth_mid, by = Group, k = 4), 
               family = binomial, data = dt_agg, method = "ML")

m_full  <- gam(cbind(orn_sp, unorn_sp) ~ Group + Region + 
                 s(SFT_mid, by = Group, k = 4) + 
                 s(SST_mid, by = Group, k = 4) + 
                 s(Depth_mid, by = Group, k = 4), 
               family = binomial, data = dt_agg, method = "ML")

model_results <- model.sel(m_null, m_sft, m_sst, m_depth, m_full)
print(model_results)

# Refit the full model using REML for accurate parameter estimation and Deviance partitioning
m_best_reml <- gam(formula(m_full), family = binomial, data = dt_agg, method = "REML")

################################################################################
## SECTION 4: Relative Importance (Deviance Partitioning)
################################################################################
cat("\n--- SECTION 4: Relative Importance of Drivers (Deviance Partitioning) ---\n")

full_dev <- summary(m_best_reml)$dev.expl
cat("Full Model Deviance Explained:", round(full_dev * 100, 2), "%\n\n")

# Calculate deviance lost when dropping each specific term
drop_group  <- update(m_best_reml, . ~ . - Group - s(SFT_mid, by=Group, k=4) - s(SST_mid, by=Group, k=4) - s(Depth_mid, by=Group, k=4) + s(SFT_mid, k=4) + s(SST_mid, k=4) + s(Depth_mid, k=4))
drop_region <- update(m_best_reml, . ~ . - Region)
drop_sft    <- update(m_best_reml, . ~ . - s(SFT_mid, by=Group, k=4))
drop_sst    <- update(m_best_reml, . ~ . - s(SST_mid, by=Group, k=4))
drop_depth  <- update(m_best_reml, . ~ . - s(Depth_mid, by=Group, k=4))

# SAFE CREATION of importance dataframe without sorting yet
importance_raw <- data.frame(
  Predictor_Raw = c("Clade", "Region", "SFT", "SST", "Depth"),
  Deviance_Lost_Pct = c(
    (full_dev - summary(drop_group)$dev.expl) * 100,
    (full_dev - summary(drop_region)$dev.expl) * 100,
    (full_dev - summary(drop_sft)$dev.expl) * 100,
    (full_dev - summary(drop_sst)$dev.expl) * 100,
    (full_dev - summary(drop_depth)$dev.expl) * 100
  )
)

# Now safely arrange and format for plotting
importance_df <- importance_raw %>% 
  arrange(desc(Deviance_Lost_Pct)) %>%
  mutate(Predictor = factor(Predictor_Raw, levels = rev(Predictor_Raw)))

print(importance_df)

################################################################################
## SECTION 5: Robustness Check - Beta Regression GAM
################################################################################
cat("\n--- SECTION 5: Robustness Check (Beta Regression) ---\n")

N_bins <- nrow(dt_agg)

dt_beta <- dt_agg %>%
  mutate(
    raw_prop = orn_sp / total_sp,
    # Squeeze transformation (Smithson & Verkuilen, 2006)
    prop_beta = (raw_prop * (N_bins - 1) + 0.5) / N_bins
  )

m_beta_full <- gam(prop_beta ~ Group + Region + 
                     s(SFT_mid, by = Group, k = 4) + 
                     s(SST_mid, by = Group, k = 4) + 
                     s(Depth_mid, by = Group, k = 4), 
                   family = betar(link = "logit"), 
                   weights = total_sp, # Weighting by species richness
                   data = dt_beta, 
                   method = "REML")

cat("\nFull Beta Model Deviance Explained:", round(summary(m_beta_full)$dev.expl * 100, 2), "%\n")

# Diagnostics (Uncomment to view in interactive session)
# par(mfrow=c(2,2))
# gam.check(m_beta_full)

################################################################################
## SECTION 6: High-Quality Visualizations
################################################################################
cat("\n--- SECTION 6: Generating Visualizations ---\n")

group_cols <- c("Bivalvia" = "#0072B2", "Brachiopoda" = "#D55E00")

# ==============================================================================
# FIGURE 3: Relative Importance Plot
# ==============================================================================
p_importance <- ggplot(importance_df, aes(x = Predictor, y = Deviance_Lost_Pct)) +
  geom_segment(aes(x = Predictor, xend = Predictor, y = 0, yend = Deviance_Lost_Pct), color = "grey70", linewidth = 1.5) +
  geom_point(size = 5, color = "#0072B2") +
  geom_text(aes(label = sprintf("%.2f%%", Deviance_Lost_Pct)), hjust = -0.4, size = 4.5) +
  coord_flip() +
  scale_y_continuous(limits = c(0, max(importance_df$Deviance_Lost_Pct) * 1.2), expand = expansion(mult = c(0, 0.1))) +
  labs(
    x = "",
    y = "Unique Deviance Explained (%)\n(Loss in explanatory power when term is removed)",
    title = "Relative Importance of Macroecological Drivers"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.y = element_text(face = "bold", color = "black"),
    axis.text.x = element_text(color = "black"),
    plot.title = element_text(face = "bold")
  )

# ==============================================================================
# FIGURE 2: Marginal Effects (Combined SFT & Depth)
# ==============================================================================
dt_plot_agg <- dt_agg %>% mutate(raw_prop = orn_sp / total_sp)

# --- 1. SFT Panel ---
pred_sft <- ggpredict(m_best_reml, terms = c("SFT_mid [all]", "Group"))
df_pred_sft <- as.data.frame(pred_sft) %>% rename(Group = group) 

p_sft <- ggplot() +
  geom_point(data = dt_plot_agg, aes(x = SFT_mid, y = raw_prop, size = total_sp, fill = Group), 
             shape = 21, color = "white", alpha = 0.4, stroke = 0.5) +
  geom_ribbon(data = df_pred_sft, aes(x = x, ymin = conf.low, ymax = conf.high, fill = Group), 
              alpha = 0.25, color = NA) +
  geom_line(data = df_pred_sft, aes(x = x, y = predicted, color = Group), linewidth = 1.2) +
  scale_color_manual(name = "Clade", values = group_cols) + 
  scale_fill_manual(name = "Clade", values = group_cols) +
  # Force proper guide rendering for size
  scale_size_continuous(
    name = "Species Richness", 
    range = c(1, 6),
    guide = guide_legend(override.aes = list(fill = "grey60", color = "grey30", alpha = 0.8))
  ) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0), labels = percent_format(accuracy = 1)) +
  labs(x = "Sea Floor Temperature (SFT, °C)", 
       y = "Predicted Proportion of Ornamented Species",
       title = "Marginal Effect of SFT",
       tag = "a") + 
  theme_classic(base_size = 14) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# --- 2. Depth Panel ---
pred_depth <- ggpredict(m_best_reml, terms = c("Depth_mid [all]", "Group"))
df_pred_depth <- as.data.frame(pred_depth) %>% rename(Group = group) 

p_depth <- ggplot() +
  geom_point(data = dt_plot_agg, aes(x = Depth_mid, y = raw_prop, size = total_sp, fill = Group), 
             shape = 21, color = "white", alpha = 0.4, stroke = 0.5) +
  geom_ribbon(data = df_pred_depth, aes(x = x, ymin = conf.low, ymax = conf.high, fill = Group), 
              alpha = 0.25, color = NA) +
  geom_line(data = df_pred_depth, aes(x = x, y = predicted, color = Group), linewidth = 1.2) +
  scale_color_manual(name = "Clade", values = group_cols) + 
  scale_fill_manual(name = "Clade", values = group_cols) +
  # Identical guide to allow patchwork collection
  scale_size_continuous(
    name = "Species Richness", 
    range = c(1, 6),
    guide = guide_legend(override.aes = list(fill = "grey60", color = "grey30", alpha = 0.8))
  ) +
  scale_x_continuous(limits = c(0, 1000), breaks = seq(0, 1000, 200)) + 
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(x = "Depth (m)", 
       y = NULL, 
       title = "Marginal Effect of Depth",
       tag = "b") + 
  theme_classic(base_size = 14) + 
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_blank(),  
    axis.ticks.y = element_blank(), 
    axis.line.y = element_blank()   
  )

# --- 3. Combine SFT and Depth using patchwork ---
p_marginal_combined <- (p_sft | p_depth) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",          
    legend.title = element_text(face = "bold"),
    legend.box = "horizontal",
    legend.margin = margin(t = 10, b = 10)
  )

# ==============================================================================
# Print and Save
# ==============================================================================
print(p_importance)
print(p_marginal_combined)

# Adjusted to GEB journal standard width (approx. 11-12 inches max for double-column legibility)
ggsave("Fig_Relative_Importance.png", p_importance, width = 8, height = 4.5, dpi = 600)
ggsave("Fig_Marginal_Combined.png", p_marginal_combined, width = 11, height = 5.5, dpi = 600)

# ==============================================================================
# SUPPLEMENTARY FIGURE: Marginal Effects by Region (Facet Wrap)
# ==============================================================================
cat("Generating Regional Breakdown Plots...\n")

# --- 1. SFT by Region ---
# Add "Region" as the third term to extract predictions for each region separately
pred_sft_reg <- ggpredict(m_best_reml, terms = c("SFT_mid [all]", "Group", "Region"))
df_pred_sft_reg <- as.data.frame(pred_sft_reg) %>% rename(Group = group, Region = facet) 

p_sft_reg <- ggplot() +
  geom_point(data = dt_plot_agg, aes(x = SFT_mid, y = raw_prop, size = total_sp, fill = Group), 
             shape = 21, color = "white", alpha = 0.4, stroke = 0.5) +
  geom_ribbon(data = df_pred_sft_reg, aes(x = x, ymin = conf.low, ymax = conf.high, fill = Group), 
              alpha = 0.25, color = NA) +
  geom_line(data = df_pred_sft_reg, aes(x = x, y = predicted, color = Group), linewidth = 1.2) +
  scale_color_manual(name = "Clade", values = group_cols) + 
  scale_fill_manual(name = "Clade", values = group_cols) +
  scale_size_continuous(name = "Species Richness", range = c(1, 6), 
                        guide = guide_legend(override.aes = list(fill = "grey60", color = "grey30", alpha = 0.8))) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0), labels = percent_format(accuracy = 1)) +
  facet_wrap(~ Region, ncol = 3) + # Split into 3 columns by Region
  labs(x = "Sea Floor Temperature (SFT, °C)", 
       y = "Predicted Proportion",
       title = "a) Marginal Effect of SFT by Region") + 
  theme_bw(base_size = 14) + # theme_bw looks better for faceted plots
  theme(plot.title = element_text(face = "bold"),
        strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(face = "bold", size = 12))

# --- 2. Depth by Region ---
pred_depth_reg <- ggpredict(m_best_reml, terms = c("Depth_mid [all]", "Group", "Region"))
df_pred_depth_reg <- as.data.frame(pred_depth_reg) %>% rename(Group = group, Region = facet) 

p_depth_reg <- ggplot() +
  geom_point(data = dt_plot_agg, aes(x = Depth_mid, y = raw_prop, size = total_sp, fill = Group), 
             shape = 21, color = "white", alpha = 0.4, stroke = 0.5) +
  geom_ribbon(data = df_pred_depth_reg, aes(x = x, ymin = conf.low, ymax = conf.high, fill = Group), 
              alpha = 0.25, color = NA) +
  geom_line(data = df_pred_depth_reg, aes(x = x, y = predicted, color = Group), linewidth = 1.2) +
  scale_color_manual(name = "Clade", values = group_cols) + 
  scale_fill_manual(name = "Clade", values = group_cols) +
  scale_size_continuous(name = "Species Richness", range = c(1, 6)) +
  scale_x_continuous(limits = c(0, 1000), breaks = seq(0, 1000, 200)) + 
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0), labels = percent_format(accuracy = 1)) +
  facet_wrap(~ Region, ncol = 3) + 
  labs(x = "Depth (m)", 
       y = "Predicted Proportion", 
       title = "b) Marginal Effect of Depth by Region") + 
  theme_bw(base_size = 14) + 
  theme(plot.title = element_text(face = "bold"),
        strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(face = "bold", size = 12))

# --- 3. Combine and Save ---
# Stack them vertically (SFT on top, Depth on bottom)
p_marginal_reg_combined <- (p_sft_reg / p_depth_reg) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",          
        legend.title = element_text(face = "bold"))

print(p_marginal_reg_combined)

# Save as a tall supplementary figure
ggsave("Fig_S8_Marginal_by_Region.png", p_marginal_reg_combined, width = 12, height = 9, dpi = 600)

cat("\n--- PIPELINE COMPLETE ---\n")
