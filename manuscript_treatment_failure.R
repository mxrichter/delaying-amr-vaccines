################################################################################
# PLOT SI FIGURE S12 (AGGREGATED TREATMENT FAILURE)
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
message("[ Plotting Figure 3 and SI Figure S12 ]")

# ---- Set up data import ----

# Initialize lists for outputs and options
o = output_list = list()
o$pth = list()

o$pth$output = file.path(getwd(), "output")

# Load mapping table for facet labels (transmission, vaccination target group)
fig3_mapping = fread("./input/scenarios/figure3_mapping.csv")
list_scenarios = fig3_mapping$scenario

# Construct paths to manuscript data and parameter tables
list_manuscript_pth = paste0(getwd(), "/output/", list_scenarios, 
                             "/8_manuscript/plot_data.rds")
list_parameters_pth = paste0(getwd(), "/output/", list_scenarios, 
                             "/1_parameters/sample_0/all_scenarios.rds")

# Path to folder with aggregated figures (more than one scenario)
fig_output_pth = file.path(o$pth$output, "publication_figures")

# Create folder if it does not exist
if(!dir.exists(fig_output_pth))
  dir.create(fig_output_pth)

# Path to aggregated plotting data file
plot_dt_pth = file.path(fig_output_pth, "manuscript_treatment_failure.rds")

message("  * Load simulation data")

# ---- Collate all outputs ----

if(!file.exists(plot_dt_pth)) {
  # ... if aggregated plotting data file does not yet exist
  
  # Loop through all scenarios
  for(i in seq_along(list_manuscript_pth)) {
    
    # Load data if file exists
    if(file.exists(list_manuscript_pth[i]))
      output_dt = readRDS(list_manuscript_pth[i])$relative
    else
      stop("! File does not exist: ", list_manuscript_pth[i])
    
    # Exclude scenarios with potential elimination for mass vaccination scenarios
    if(grepl("altvax_massvax", list_scenarios[i])) {
      # In mass vaccination scenario with longer-lasting vaccine, PfPR2-10 40%
      output_dt = output_dt[pfpr >= 0.4, ]
    } else if(grepl("massvax", list_scenarios[i])) {
      # In RTS,S/AS01-like mass vaccination scenario, PfPR2-10 = {20%, 40%}
      output_dt = output_dt[pfpr %in% c(0.2, 0.4)]
    }
    
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

# Restrict analysis to vaccine-eligible children
plot_dt = plot_dt[age_group == 1, ]

# ---- Actual plotting from here on ----

# Color palette for plotting
pal = c("#094061", "#085685", "#2778a8", "#539dc9", "#a9c7d9") %>% rev()

# Calculate means and quantiles
plot_agg_dt = plot_dt[, .(mean = 100*mean(change),
            lower = 100*quantile(change, 0.025),
            upper = 100*quantile(change, 0.975)),
        by = list(scenario, access, id, prevalence_resistant, pfpr, coverage)]

# Add facet labels to data.table
plot_agg_dt = fig3_mapping[plot_agg_dt, on = "scenario"]

# Convert to factors for prettier plotting
plot_agg_dt[, transmission := factor(transmission, level = c("seasonal", "perennial"),
                                 label = c("Seasonal transmission", "Perennial transmission"))]

plot_agg_dt[, target_group := factor(target_group, levels = c("children", "all"),
                                 labels = c("Children only", "Mass vaccination"))]

plot_agg_dt[, vaccine := factor(vaccine, levels = c("rtssas01", "rtssas01like","longerlasting"),
                            labels = c("RTS,S/AS01", "RTS,S/AS01-like", "Longer-lasting\nvaccine"))]

# Format parasite prevalence
pfpr_vals   = unique(plot_agg_dt$pfpr) %>% unique() %>% sort()

pfpr_labels = paste0("italic(Pf)PR[2-10]: ", format(round(100*pfpr_vals)),"%")

pfpr_labels = sprintf("paste(italic(Pf)*PR[2-10], ': %s%%')", round(100 * pfpr_vals))

names(pfpr_labels) = pfpr_vals

freq_vals = unique(plot_agg_dt$prevalence_resistant) %>% sort()
freq_labels = scales::percent(freq_vals)

plot_agg_dt[, prevalence_resistant := factor(prevalence_resistant,
                                             levels = freq_vals,
                                             labels = freq_labels)]

coverage_vals = unique(plot_agg_dt$coverage) %>% sort()
coverage_labels = scales::percent(coverage_vals)

plot_agg_dt[, coverage := factor(coverage,
                                             levels = coverage_vals,
                                             labels = coverage_labels)]

# Basic plotting function (main text and SI plot differ by line break in y-axis label)
plot_fn = function(dt, si = FALSE) {
  gg_output = ggplot(dt, aes(x=prevalence_resistant, y=mean, group = coverage, fill = coverage)) +
    geom_bar(stat = "identity", width = 0.8, position = position_dodge(0.8, preserve = "single")) +
    geom_errorbar(aes(ymin=lower, ymax=upper), width = 0.4, 
                  position = position_dodge(0.8, preserve = "single")) +
    scale_fill_manual(name = "Vaccination coverage", values = pal) +
    scale_x_discrete(name = "Frequency of drug-resistant genotype") +
    theme_amr() +
    theme(strip.placement = "outer",
          axis.text.x = element_text(color = "black", size = 16,
                                     margin = margin(t = 10, b = 5)),
          axis.title.x = element_text(margin = margin(t = 10, b = 0)),
          axis.title.y = element_text(margin = margin(r = 10, l = 5)),
          legend.box.margin = margin(t = -5, b = 5),
          panel.spacing.x = unit(15, units = "pt"),
          panel.spacing.y = unit(15, units = "pt"),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
          plot.margin = margin(t = 5, r = 20, b = 5, l = 5))
  
  # If this is for SI, use certain order of line breaks
  if(si == FALSE) {
    gg_output = gg_output +
      scale_y_continuous(name = "Reduction in 30-day\ntreatment failure rate\n(in percentage points)")
  } else {
    gg_output = gg_output +
      scale_y_continuous(name = "Reduction in 30-day treatment failure rate (in percentage points)")
  }
}

# ---- Plot for main text ---

# Seasonal transmission, moderate access to treatment, very high resistance, 
# RTS,S/AS01 deployed to children
plot_dt_fig3 = plot_agg_dt %>%
  filter(access == 0.2412, id == "very_high",
         transmission == "Seasonal transmission", target_group == "Children only",
         vaccine == "RTS,S/AS01") 

# Assemble plot using auxiliary function
plot_fig3 = plot_fn(plot_dt_fig3) +
  facet_nested(. ~ pfpr,
               nest_line = element_line(color = "black", linewidth = 0.5),
               labeller = labeller(pfpr = as_labeller(pfpr_labels, label_parsed)),
               switch = "y")

# Save plot into publication figures folder
ggsave(file.path(fig_output_pth, "MANUSCRIPT_Figure-3_Treatment-failure-rates.png"), plot_fig3, width = 16, height = 5)  

# Extract treatment failure rate values mentioned in main text

# Values: largest reduction
# plot_dt_fig3[mean == max(mean)]

# Values: largest increase
# plot_dt_fig3[mean == min(mean)]

# ---- Plot for SI ----

# All transmission settings, moderate access to treatment, very high resistance,
# all vaccine types, all vaccine target groups
plot_dt_si = plot_agg_dt %>%
  filter(access == 0.2412, id == "very_high")

# Assemble plot using auxiliary function
plot_si = plot_fn(plot_dt_si, si = TRUE) +
  facet_nested(transmission + target_group + vaccine ~ pfpr,
               nest_line = element_line(color = "black", linewidth = 0.5),
               scales = "free_y", 
               labeller = labeller(pfpr = as_labeller(pfpr_labels, label_parsed)),
               switch = "y",
               render_empty = FALSE)

# Save plot into publication figures folder
ggsave(file.path(fig_output_pth, "SI_Figure-S12_Treatment-failure-rates.png"), plot_si, width = 16, height = 18)  
