################################################################################
# PLOT SI FIGURE S11 (AGGREGATED VACCINE EFFECTIVENESS)
################################################################################

# --- Some preparations before plotting ----

# Set working directory properly
setwd(getSrcDirectory(function(){})[1])

# Source R libraries and auxiliary functions
source("dependencies.R")

# Clear console
clc()

# Display FYI messages to user
message("[ THE IMPACT OF VACCINES ON DRUG-RESISTANT P. FALCIPARUM ]")
message("[ Plotting SI Figure S11 ]")

# ---- Set up data import ----

# Initialize lists for outputs and options
o = output_list = list()
o$pth = list()

# Load mapping table for facet labels (transmission, vaccination target group)
fig2_mapping = fread("./input/scenarios/figure2_mapping.csv")
list_scenarios = fig2_mapping$scenario

# Load mapping table for resistance scenarios
resist_scenarios = fread("./input/scenarios/resistance_scenarios_lim.txt")
id_order = resist_scenarios$id
id_order_label = resist_scenarios$name

# Construct paths to manuscript data and parameter tables
list_manuscript_pth = paste0(getwd(), "/output/", list_scenarios, 
                             "/8_manuscript/plot_data.rds")
list_parameters_pth = paste0(getwd(), "/output/", list_scenarios, 
                             "/1_parameters/sample_0/all_scenarios.rds")

# Path to folder with aggregated figures (more than one scenario)
fig_output_pth = file.path(getwd(), "output", "publication_figures")

# Create folder if it does not exist
if(!dir.exists(fig_output_pth))
  dir.create(fig_output_pth)

# Path to aggregated plotting data file
plot_dt_pth = file.path(fig_output_pth, "manuscript_ve.rds")

message("  * Load simulation data")

# ---- Collate all outputs ----

if(!file.exists(plot_dt_pth)) {
  # ... if aggregated plotting data file does not yet exist
  
  # Loop through all scenarios
  for(i in seq_along(list_manuscript_pth)) {
    
    # Load data if file exists
    if(file.exists(list_manuscript_pth[i]))
      output_dt = readRDS(list_manuscript_pth[i])$results
    else
      stop("! File does not exist: ", list_manuscript_pth[i])
    
    # Store scenario parameters
    parameters_dt = readRDS(list_parameters_pth[i])
    
    # Store scenario name
    set(output_dt, j = "scenario", value = list_scenarios[i])
    
    output_list[[i]] = output_dt
    
  }
  
  # Concatenate results data
  plot_dt = rbindlist(output_list)
  
  # Save aggregated plotting data file
  saveRDS(plot_dt, plot_dt_pth)
  
} else {
  
  # ... read file if it already exists
  plot_dt = readRDS(plot_dt_pth)
  
}

message("  * Prepare data for plotting")

# ---- Actual plotting from here on ----

# Add facet labels to data.table
plot_dt = fig2_mapping[plot_dt, on = "scenario"]

# Convert to factor variables in correct order
plot_dt[, id := factor(id, levels = id_order,
                       labels = id_order_label)]

# Calculate means and quantiles
plot_agg_dt = plot_dt[, .(mean = mean(change),
            lower = quantile(change, 0.025),
            upper = quantile(change, 0.975)),
        by = list(transmission, target_group, vaccine, indicator, access, pfpr, coverage, id)]

# Convert to factors for prettier plotting
plot_agg_dt[, transmission := factor(transmission, level = c("seasonal", "perennial"),
                                 label = c("Seasonal transmission", "Perennial transmission"))]

plot_agg_dt[, target_group := factor(target_group, levels = c("children", "all"),
                                 labels = c("Children only", "Mass vaccination"))]

plot_agg_dt[, vaccine := factor(vaccine, levels = c("rtssas01", "rtssas01like","longerlasting"),
                            labels = c("RTS,S/AS01", "RTS,S/AS01-like", "Longer-lasting\nvaccine"))]

# Only use low and moderate access to treatment
plot_agg_dt = plot_agg_dt[access %in% c(0.0356, 0.2412)]

plot_agg_dt[, access := ifelse(access == 0.0356, "Low access to treatment (10%)", "Moderate access to treatment (50%)")]

# Format parasite prevalence
pfpr_vals   = unique(plot_agg_dt$pfpr) %>% unique() %>% sort()

pfpr_labels = paste0("italic(Pf)PR[2-10]: ", format(round(100*pfpr_vals)),"%")

pfpr_labels = sprintf("paste(italic(Pf)*PR[2-10], ': %s%%')", round(100 * pfpr_vals))

names(pfpr_labels) = pfpr_vals

coverage_levels = plot_dt$coverage %>% unique() %>% sort()
coverage_labels = scales::percent(coverage_levels)

plot_agg_dt[, coverage := factor(coverage,
                                 levels = coverage_levels,
                                 labels = coverage_labels)]

# Color palette for plotting
pal = c("#468AB2", "#BF3227", "#EDCD64")

message("  * Build plot")

# Build actual plot with ggplot2
plot_ve =  
   plot_agg_dt[coverage != "100%", ] %>%
  ggplot(aes(x=id, y=mean, fill = id, group = factor(coverage), alpha = factor(coverage))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  scale_y_continuous(name = "Reduction in cumulative uncomplicated cases in children under five years of age",
                     labels = scales::percent,
                     expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(name = NULL)+
    scale_fill_manual(values = pal, name = "Resistance scenario") +
  geom_errorbar(aes(ymin = lower, ymax=upper), width = .3,
                position = position_dodge(width = 0.9), alpha = 1) +
  guides(fill = guide_legend(reverse = FALSE, ncol = 2, byrow = TRUE),
         alpha = guide_legend(ncol = 2, byrow = TRUE)) +
  facet_nested(transmission + target_group + vaccine ~ access + pfpr, 
             switch = "y", scales = "free_y",
             nest_line = element_line(colour = "black", linewidth = 0.5),
             labeller = labeller(pfpr = as_labeller(pfpr_labels, label_parsed))) +
     scale_alpha_manual(values = c(0.5, 0.75, 1.0),
                        name = "Vaccination coverage") +
  theme_amr() +
  theme(strip.placement = "outer", 
        strip.background = element_rect(fill = "white",
                                        color = NA),
        strip.text = element_text(size = 16),
        panel.background = element_rect(color = "black", linewidth = 0.8),
        legend.position = "bottom",
        # axis.text.x = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 16, color = "black"),
        axis.title.x = element_text(margin = margin(t = 15, b = 0)),
        legend.title = element_text(size = 16, color = "black"),
        legend.text = element_text(size = 16, color = "black"),
        legend.margin = margin(t = 0),
        panel.spacing.x = unit(15, units = "pt"),
        panel.spacing.y = unit(15, units = "pt"),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
        plot.margin = margin(t = 5, r = 20, b = 5, l = 5))

# Save plot to aggregated output folder
ggsave(file.path(fig_output_pth, "SI_Figure-S11_Vaccine-effectiveness.png"), plot_ve, width = 16, height = 18)  
  
message("[ Finished ]")