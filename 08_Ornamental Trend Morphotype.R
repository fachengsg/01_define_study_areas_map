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
# 2. Ornamented Rate Calculation Function
# -----------------------------------------------------------------------------
calc_rate <- function(data, var, bin_width, type){
  dt <- copy(data)
  
  if(var=="SFT"){ dt[, bin := floor(SFT/bin_width)*bin_width] }
  if(var=="Depth"){ dt[, bin := floor(Depth/bin_width)*bin_width] }
  
  if(type=="occurrence"){
    out <- dt %>%
      group_by(Group, bin) %>%
      summarise(rate = mean(Orn_bin, na.rm=TRUE), n = n(), .groups="drop")
  }
  
  if(type=="species"){
    dt_sp <- dt %>%
      group_by(Group, species, bin) %>%
      summarise(Orn_bin = mean(Orn_bin), .groups="drop")
    
    out <- dt_sp %>%
      group_by(Group, bin) %>%
      summarise(rate = mean(Orn_bin), n = n(), .groups="drop")
  }
  
  out$type <- type
  out <- out %>% filter(n >= MIN_BIN_N)
  return(out)
}

# -----------------------------------------------------------------------------
# 3. Plotting Function (Preserved Original Layout)
# -----------------------------------------------------------------------------
make_plot <- function(data, var, bin_width, region_name="All", show_y_label = TRUE){
  dt <- data
  if(region_name!="All"){ dt <- dt %>% filter(Region==region_name) }
  
  occ <- calc_rate(dt, var, bin_width, "occurrence")
  sp  <- calc_rate(dt, var, bin_width, "species")
  df  <- bind_rows(occ, sp)
  
  xlabel <- ifelse(var=="SFT", "SFT (°C)", "Depth (m)")
  if(var=="SFT"){ xscale <- scale_x_continuous(limits=c(0,25)) }
  if(var=="Depth"){ xscale <- scale_x_continuous(limits=c(0,5000)) }
  
  p <- ggplot(df, aes(bin, rate, colour=Group, linetype=type)) +
    geom_point(aes(fill=Group), alpha=0.35, size=1.7, shape=21) +
    geom_smooth(aes(fill = Group), method="gam", formula=y~s(x,k=5),
                method.args = list(optimizer = "efs"), se=TRUE, size=1.1) +
    scale_colour_manual(values=group_cols, name="Group") +
    scale_fill_manual(values=group_cols, guide="none") +
    scale_linetype_manual(
      name="Data level",
      values=c(species="solid", occurrence="21"),
      labels=c(species="Species mean", occurrence="Occurrence mean")
    ) +
    xscale +
    scale_y_continuous(limits=c(0,1), expand=c(0,0)) +
    labs(x=xlabel, y=ifelse(show_y_label, "Ornamented rate", "")) +
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
plot1 <- make_plot(dt_final, "SFT", 1, "All", show_y_label = TRUE) + ggtitle("All") + theme(plot.title = element_text(hjust = 0.5))
plot2 <- make_plot(dt_final, "SFT", 1, "Atl-North", show_y_label = FALSE) + ggtitle("Atl-North") + theme(plot.title = element_text(hjust = 0.5))
plot3 <- make_plot(dt_final, "SFT", 1, "Atl-South", show_y_label = FALSE) + ggtitle("Atl-South") + theme(plot.title = element_text(hjust = 0.5))
plot4 <- make_plot(dt_final, "SFT", 1, "Mediterranean", show_y_label = FALSE) + ggtitle("Mediterranean") + theme(plot.title = element_text(hjust = 0.5))

plot5 <- make_plot(dt_final, "Depth", 100, "All", show_y_label = TRUE) + ggtitle("All") + theme(plot.title = element_text(hjust = 0.5))
plot6 <- make_plot(dt_final, "Depth", 100, "Atl-North", show_y_label = FALSE) + ggtitle("Atl-North") + theme(plot.title = element_text(hjust = 0.5))
plot7 <- make_plot(dt_final, "Depth", 100, "Atl-South", show_y_label = FALSE) + ggtitle("Atl-South") + theme(plot.title = element_text(hjust = 0.5))
plot8 <- make_plot(dt_final, "Depth", 100, "Mediterranean", show_y_label = FALSE) + ggtitle("Mediterranean") + theme(plot.title = element_text(hjust = 0.5))

row_SFT   <- plot1 + plot2 + plot3 + plot4 + plot_layout(ncol = 4)
row_Depth <- plot5 + plot6 + plot7 + plot8 + plot_layout(ncol = 4)

combined_layout <- (row_SFT / row_Depth) +
  plot_layout(guides = "collect") +
  theme(legend.position = "bottom", legend.box = "horizontal") & 
  theme(legend.position = "bottom")

ggsave(file.path(OUT_DIR, "Fig_S4_Raw_Combined_Trends.png"), combined_layout, width = 18, height = 8, dpi = 600)



################################################################################
## Supplementary Analysis: Bivalvia Morphotype Proportions (a-e)
## Target Journal: Global Ecology and Biogeography
## Output: 
## - Fig_S5_Bivalvia_Species.png
## - Fig_S6_Bivalvia_Occurrence.png
## - Table_S_Ornament_Type_Proportions.csv
################################################################################

library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)
library(cowplot)

# ------------------------------------------------------------------
# 1. Setup and Data Loading
# ------------------------------------------------------------------
OUT_DIR <- "outputs"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

if (!exists("dt_final")) {
  dt_final <- readRDS("09_final_research_cleaned.rds")
}
dt_final <- as.data.table(dt_final)

BIN_WIDTH_SFT   <- 1
BIN_WIDTH_DEPTH <- 100
MIN_BIN_N       <- 5   

# Colour palette for Bivalvia ornaments (a–e)
biv_colors <- c("a" = "#1f77b4", "b" = "#ff7f0e", "c" = "#2ca02c",
                "d" = "#d62728", "e" = "#9467bd")

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
  
  occ_counts <- dt[, .(n_occ = .N), by = .(Group, bin, Ornament)]
  total_occ  <- dt[, .(total = .N), by = .(Group, bin)]
  result <- merge(occ_counts, total_occ, by = c("Group", "bin"))
  result[, `:=`(prop = n_occ / total, Region = region_name, Variable = var, Level = "occurrence")]
  result <- result[total >= MIN_BIN_N]
  result[order(Group, bin, Ornament)]
}

calc_sp_prop <- function(dt, var, bin_width, region_name = "All") {
  dt <- copy(dt)
  if (region_name != "All") dt <- dt[Region == region_name]
  dt <- add_bins(dt, var, bin_width)
  
  sp_ornament <- dt[, .N, by = .(species, Group, Region, bin, Ornament)
  ][, .SD[which.max(N)], by = .(species, Group, Region, bin)
  ][, .(species, Group, Region, bin, Ornament)]
  
  sp_counts <- sp_ornament[, .(n_sp = .N), by = .(Group, bin, Ornament)]
  total_sp   <- sp_ornament[, .(total = .N), by = .(Group, bin)]
  result <- merge(sp_counts, total_sp, by = c("Group", "bin"))
  result[, `:=`(prop = n_sp / total, Region = region_name, Variable = var, Level = "species")]
  result <- result[total >= MIN_BIN_N]
  result[order(Group, bin, Ornament)]
}

# ------------------------------------------------------------------
# 3. Compute Data Tables
# ------------------------------------------------------------------
regions   <- c("All", "Atl-North", "Atl-South", "Mediterranean")
variables <- c("SFT", "Depth")

prop_tables <- list()
for (reg in regions) {
  for (var in variables) {
    bw <- ifelse(var == "SFT", BIN_WIDTH_SFT, BIN_WIDTH_DEPTH)
    prop_tables[[paste(reg, var, "occ", sep = "_")]] <- calc_occ_prop(dt_final, var, bw, reg)
    prop_tables[[paste(reg, var, "sp",  sep = "_")]] <- calc_sp_prop(dt_final, var, bw, reg)
  }
}
all_proportions <- rbindlist(prop_tables, use.names = TRUE, fill = TRUE)
fwrite(all_proportions, file.path(OUT_DIR, "Table_S_Ornament_Type_Proportions.csv"))

biv_data <- all_proportions[Group == "Bivalvia"]

SFT_LIM   <- c(0, 25)
DEPTH_LIM <- c(0, 5000)

# ------------------------------------------------------------------
# 4. Plotting Engine 
# ------------------------------------------------------------------
create_bivalvia_figure <- function(data, level_name) {
  panels <- list()
  for (v in c("SFT", "Depth")) {
    for (r in c("All", "Atl-North", "Atl-South", "Mediterranean")) {
      sub <- data[Variable == v & Region == r]
      
      xlim <- if(v == "SFT") SFT_LIM else DEPTH_LIM
      xbreaks <- if(v == "SFT") seq(0, 25, by = 5) else seq(0, 5000, by = 1000)
      xlab <- if(v == "SFT") "SFT (°C)" else "Depth (m)"
      
      # Minimal subplot title
      title_text <- r 
      
      p <- ggplot(sub, aes(x = bin, y = prop, colour = Ornament)) +
        geom_line(linewidth = 0.9) + 
        geom_point(size = 1.2, alpha = 0.7) +
        # Keep legend in original layers, patchwork will collect them automatically
        scale_colour_manual(values = biv_colors, name = "Bivalvia morphotype") +
        scale_x_continuous(limits = xlim, breaks = xbreaks) +
        scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
        labs(
          title = title_text,
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
  
  # Define the order of the 8 panels
  panel_order <- c("SFT_All", "SFT_Atl-North", "SFT_Atl-South", "SFT_Mediterranean",
                   "Depth_All", "Depth_Atl-North", "Depth_Atl-South", "Depth_Mediterranean")
  
  # Assemble the main title at top-left
  main_title <- paste0("Bivalvia (", level_name, ")")
  
  # Use native patchwork to combine plots, collect legends, and add main title
  combined_plot <- wrap_plots(panels[panel_order], ncol = 4) +
    plot_layout(guides = "collect") + # Automatically collect common legends from all subplots
    plot_annotation(
      title = main_title, # Add main title
      tag_levels = "a",   # Add a-h labels
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0),
        plot.tag = element_text(size = 12, face = "bold")
      )
    ) & 
    theme(legend.position = "bottom") # Force collected legend to the bottom
  
  return(combined_plot)
}

# ------------------------------------------------------------------
# 5. Build & Save (Super Clean!)
# ------------------------------------------------------------------
# Directly call the function to generate final plots, no need to extract legends separately
final_species    <- create_bivalvia_figure(biv_data[Level == "species"], "species")
final_occurrence <- create_bivalvia_figure(biv_data[Level == "occurrence"], "occurrence")

print(final_species)
print(final_occurrence)

ggsave(file.path(OUT_DIR, "Fig_S5_Bivalvia_Species.png"), final_species, width = 10, height = 8, dpi = 600, bg = "white")
ggsave(file.path(OUT_DIR, "Fig_S6_Bivalvia_Occurrence.png"), final_occurrence, width = 10, height = 8, dpi = 600, bg = "white")

cat("\n--- Bivalvia Morphotype Pipeline Complete ---\n")