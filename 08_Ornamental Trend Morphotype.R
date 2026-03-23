################################################################################
## Supplementary Analysis: Raw Macroecological Ornamental Trends 
## (Species- vs. Occurrence-level comparison)
## Output: Fig_S4_Raw_Combined_Trends.png
################################################################################

library(data.table)
library(dplyr)
library(ggplot2)
library(mgcv)
library(patchwork)

# -----------------------------------------------------------------------------
# 1. Setup and Data Loading
# -----------------------------------------------------------------------------
OUT_DIR <- "outputs"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

if(!exists("dt_final")){
  dt_final <- readRDS("09_final_research_cleaned.rds")
}
dt_final <- as.data.table(dt_final)

# Parameters
group_cols <- c("Bivalvia" = "#0072B2", "Brachiopoda" = "#D55E00")
MIN_BIN_N  <- 5

# -----------------------------------------------------------------------------
# 2. Ornamented Rate Calculation Function (UPDATED to Binary_Trait)
# -----------------------------------------------------------------------------
calc_rate <- function(data, var, bin_width, type){
  dt <- copy(data)
  
  if(var=="SFT"){ dt[, bin := floor(SFT/bin_width)*bin_width] }
  if(var=="Depth"){ dt[, bin := floor(Depth/bin_width)*bin_width] }
  
  if(type=="occurrence"){
    out <- dt %>%
      group_by(Group, bin) %>%
      # FIX: Changed Orn_bin to Binary_Trait
      summarise(rate = mean(Binary_Trait, na.rm=TRUE), n = n(), .groups="drop")
  }
  
  if(type=="species"){
    dt_sp <- dt %>%
      group_by(Group, species, bin) %>%
      # FIX: Changed Orn_bin to Binary_Trait
      summarise(Binary_Trait = max(Binary_Trait, na.rm=TRUE), .groups="drop")
    
    out <- dt_sp %>%
      group_by(Group, bin) %>%
      summarise(rate = mean(Binary_Trait), n = n(), .groups="drop")
  }
  
  out$type <- type
  out <- out %>% filter(n >= MIN_BIN_N)
  return(out)
}

# -----------------------------------------------------------------------------
# 3. Plotting Function (Fixed Axis Limits and GAM Sparsity Issues)
# -----------------------------------------------------------------------------
make_plot <- function(data, var, bin_width, region_name="All", show_y_label = TRUE){
  dt <- data
  if(region_name!="All"){ dt <- dt %>% filter(Region==region_name) }
  
  occ <- calc_rate(dt, var, bin_width, "occurrence")
  sp  <- calc_rate(dt, var, bin_width, "species")
  df  <- bind_rows(occ, sp)
  
  # Determine axis limits and labels dynamically
  if(var == "SFT"){ 
    x_limits <- c(0, 25)
    xlab <- "SFT (°C)" 
  }
  if(var == "Depth"){ 
    x_limits <- c(0, 5000)
    xlab <- "Depth (m)" 
  }
  
  p <- ggplot(df, aes(x = bin, y = rate, colour = Group, linetype = type)) +
    geom_point(aes(fill = Group), alpha = 0.35, size = 1.7, shape = 21) +
    # FIX 2: Reduced k=5 to k=3 to prevent GAM failures on sparse regional data (e.g., Mediterranean Brachiopods)
    geom_smooth(aes(fill = Group), method = "gam", formula = y ~ s(x, k = 3),
                method.args = list(optimizer = "efs"), se = TRUE, linewidth = 1.1) +
    scale_colour_manual(values = group_cols, name = "Group") +
    scale_fill_manual(values = group_cols, guide = "none") +
    scale_linetype_manual(
      name = "Data level",
      values = c(species = "solid", occurrence = "21"),
      labels = c(species = "Species mean", occurrence = "Occurrence mean")
    ) +
    # FIX 1: Use coord_cartesian instead of scale_..._continuous(limits=...) to prevent data dropping
    coord_cartesian(xlim = x_limits, ylim = c(0, 1)) +
    labs(x = xlab, y = ifelse(show_y_label, "Proportion Ornamented", "")) +
    theme_bw() +
    theme(legend.position = "bottom", legend.direction = "horizontal")
  
  if(!show_y_label) {
    p <- p + theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  }
  return(p)
}

# -----------------------------------------------------------------------------
# 4. Assembly and Save
# -----------------------------------------------------------------------------
cat("\nGenerating Figure S4 (Raw Combined Trends)...\n")
plot1 <- make_plot(dt_final, "SFT", 1, "All", TRUE) + ggtitle("All") + theme(plot.title = element_text(hjust = 0.5))
plot2 <- make_plot(dt_final, "SFT", 1, "Atl-North", FALSE) + ggtitle("Atl-North") + theme(plot.title = element_text(hjust = 0.5))
plot3 <- make_plot(dt_final, "SFT", 1, "Atl-South", FALSE) + ggtitle("Atl-South") + theme(plot.title = element_text(hjust = 0.5))
plot4 <- make_plot(dt_final, "SFT", 1, "Mediterranean", FALSE) + ggtitle("Mediterranean") + theme(plot.title = element_text(hjust = 0.5))

plot5 <- make_plot(dt_final, "Depth", 100, "All", TRUE) + theme(plot.title = element_text(hjust = 0.5))
plot6 <- make_plot(dt_final, "Depth", 100, "Atl-North", FALSE) + theme(plot.title = element_text(hjust = 0.5))
plot7 <- make_plot(dt_final, "Depth", 100, "Atl-South", FALSE) + theme(plot.title = element_text(hjust = 0.5))
plot8 <- make_plot(dt_final, "Depth", 100, "Mediterranean", FALSE) + theme(plot.title = element_text(hjust = 0.5))

row_SFT   <- plot1 + plot2 + plot3 + plot4 + plot_layout(ncol = 4)
row_Depth <- plot5 + plot6 + plot7 + plot8 + plot_layout(ncol = 4)

# Keep the & theme to force the collected legend to the bottom
combined_layout <- (row_SFT / row_Depth) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.margin = margin(t = 10, b = 10)
  )

ggsave(file.path(OUT_DIR, "Fig_S5_Raw_Combined_Trends.png"), combined_layout, width = 18, height = 8.5, dpi = 600)


################################################################################
## Supplementary Analysis: Bivalvia Morphotype Proportions (a-e)
## Focus: Integrating Trait Complexity Index (TCI) visually
## Output: 
## - Fig_S5_Bivalvia_Species.png
## - Fig_S6_Bivalvia_Occurrence.png
## - Table_S_Ornament_Type_Proportions.csv
################################################################################

cat("\nProcessing Bivalvia Morphotypes (TCI Integration)...\n")

# ------------------------------------------------------------------
# 1. Setup and Data Preparation (Mapping a-e to Descriptive TCI)
# ------------------------------------------------------------------
dt_biv <- dt_final[Group == "Bivalvia"]

# Clean raw Ornament strings just to be safe
dt_biv[, Ornament_clean := tolower(trimws(as.character(Ornament)))]

# Define valid morphotypes and map them to descriptive TCI labels
orn_map <- c(
  "a" = "Smooth/Commarginal (TCI 0)",
  "b" = "Radial (TCI 1)",
  "c" = "Concentric (TCI 1)",
  "e" = "Cancellate/Reticulate (TCI 2)",
  "d" = "Spines/Nodes (TCI 3)"
)

# Filter valid codes and apply new labels
dt_biv <- dt_biv[Ornament_clean %in% names(orn_map)]
dt_biv[, Morphotype := factor(orn_map[Ornament_clean], levels = orn_map)]

# Colour palette graded by complexity (Grey -> Blues/Greens -> Purple -> Red)
biv_colors <- c(
  "Smooth/Commarginal (TCI 0)" = "#999999",
  "Radial (TCI 1)" = "#377EB8",
  "Concentric (TCI 1)" = "#4DAF4A",
  "Cancellate/Reticulate (TCI 2)" = "#984EA3",
  "Spines/Nodes (TCI 3)" = "#E41A1C"
)

# ------------------------------------------------------------------
# 2. Calculation Functions
# ------------------------------------------------------------------
add_bins <- function(dt, var, bin_width) {
  dt <- copy(dt)
  if (var == "SFT")   dt[, bin := floor(SFT / bin_width) * bin_width]
  if (var == "Depth") dt[, bin := floor(Depth / bin_width) * bin_width]
  dt
}

calc_occ_prop <- function(dt, var, bin_width, region_name = "All") {
  dt <- copy(dt)
  if (region_name != "All") dt <- dt[Region == region_name]
  dt <- add_bins(dt, var, bin_width)
  
  occ_counts <- dt[, .(n_occ = .N), by = .(bin, Morphotype)]
  total_occ  <- dt[, .(total = .N), by = .(bin)]
  result <- merge(occ_counts, total_occ, by = "bin")
  result[, `:=`(prop = n_occ / total, Region = region_name, Variable = var, Level = "occurrence")]
  result <- result[total >= MIN_BIN_N]
  result[order(bin, Morphotype)]
}

calc_sp_prop <- function(dt, var, bin_width, region_name = "All") {
  dt <- copy(dt)
  if (region_name != "All") dt <- dt[Region == region_name]
  dt <- add_bins(dt, var, bin_width)
  
  # Select the most dominant morphotype per species per bin
  sp_ornament <- dt[, .N, by = .(species, Region, bin, Morphotype)
  ][, .SD[which.max(N)], by = .(species, Region, bin)
  ][, .(species, Region, bin, Morphotype)]
  
  sp_counts <- sp_ornament[, .(n_sp = .N), by = .(bin, Morphotype)]
  total_sp   <- sp_ornament[, .(total = .N), by = .(bin)]
  result <- merge(sp_counts, total_sp, by = "bin")
  result[, `:=`(prop = n_sp / total, Region = region_name, Variable = var, Level = "species")]
  result <- result[total >= MIN_BIN_N]
  result[order(bin, Morphotype)]
}

# ------------------------------------------------------------------
# 3. Compute Data Tables
# ------------------------------------------------------------------
regions   <- c("All", "Atl-North", "Atl-South", "Mediterranean")
variables <- c("SFT", "Depth")

prop_tables <- list()
for (reg in regions) {
  for (var in variables) {
    bw <- ifelse(var == "SFT", 1, 100)
    prop_tables[[paste(reg, var, "occ", sep = "_")]] <- calc_occ_prop(dt_biv, var, bw, reg)
    prop_tables[[paste(reg, var, "sp",  sep = "_")]] <- calc_sp_prop(dt_biv, var, bw, reg)
  }
}
all_proportions <- rbindlist(prop_tables, use.names = TRUE, fill = TRUE)
fwrite(all_proportions, file.path(OUT_DIR, "Table_S_Ornament_Type_Proportions.csv"))

# ------------------------------------------------------------------
# 4. Plotting Engine 
# ------------------------------------------------------------------
create_bivalvia_figure <- function(data, level_name) {
  panels <- list()
  for (v in c("SFT", "Depth")) {
    for (r in c("All", "Atl-North", "Atl-South", "Mediterranean")) {
      sub <- data[Variable == v & Region == r]
      
      xlim <- if(v == "SFT") c(0, 25) else c(0, 5000)
      xbreaks <- if(v == "SFT") seq(0, 25, by = 5) else seq(0, 5000, by = 1000)
      xlab <- if(v == "SFT") "SFT (°C)" else "Depth (m)"
      
      p <- ggplot(sub, aes(x = bin, y = prop, colour = Morphotype)) +
        geom_line(linewidth = 0.9) + 
        geom_point(size = 1.2, alpha = 0.7) +
        scale_colour_manual(values = biv_colors, name = "Bivalvia Morphotype & TCI") +
        scale_x_continuous(limits = xlim, breaks = xbreaks) +
        scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
        labs(
          title = r,
          x = if (r == "All") xlab else "",
          y = if (v == "SFT" & r == "All") "Proportion" else ""
        ) +
        theme_bw() +
        theme(
          plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8),
          axis.title = element_text(size = 9),
          panel.grid.minor = element_blank(),
          legend.title = element_text(face = "bold")
        )
      panels[[paste(v, r, sep = "_")]] <- p
    }
  }
  
  panel_order <- c("SFT_All", "SFT_Atl-North", "SFT_Atl-South", "SFT_Mediterranean",
                   "Depth_All", "Depth_Atl-North", "Depth_Atl-South", "Depth_Mediterranean")
  
  main_title <- paste0("Bivalvia Architectural Complexity (", tools::toTitleCase(level_name), "-weighted)")
  
  combined_plot <- wrap_plots(panels[panel_order], ncol = 4) +
    plot_layout(guides = "collect") + 
    plot_annotation(
      title = main_title, 
      tag_levels = "a",   
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0),
        plot.tag = element_text(size = 12, face = "bold")
      )
    ) & 
    theme(legend.position = "bottom") 
  
  return(combined_plot)
}

# ------------------------------------------------------------------
# 5. Build & Save 
# ------------------------------------------------------------------
cat("Generating Figure S5 & S6 (Bivalvia Morphotypes)...\n")
final_species    <- create_bivalvia_figure(all_proportions[Level == "species"], "species")
final_occurrence <- create_bivalvia_figure(all_proportions[Level == "occurrence"], "occurrence")

ggsave(file.path(OUT_DIR, "Fig_S5_Bivalvia_Species.png"), final_species, width = 11, height = 8, dpi = 600, bg = "white")
ggsave(file.path(OUT_DIR, "Fig_S6_Bivalvia_Occurrence.png"), final_occurrence, width = 11, height = 8, dpi = 600, bg = "white")

cat("\n--- Bivalvia Morphotype Pipeline Complete ---\n")
