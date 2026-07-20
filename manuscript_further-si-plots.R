################################################################################
# ADDITIONAL SI PLOTS (EMULATOR PERFORMANCE, VACCINE EFFICACY)
################################################################################

# Set working directory
setwd(getSrcDirectory(function(){})[1])

# ---- Some preparations for data import and plotting ----

# Source necessary R packages and scripts
source("dependencies.R")

# Path to folder with aggregated figures (more than one scenario)
fig_output_pth = file.path(getwd(), "output", "publication_figures")

if(!dir.exists(fig_output_pth))
  dir.create(fig_output_pth)

# Define color palette
pal = c("#468AB2", "#BF3227", "#EDCD64")

# ---- Produce emulator diagnostic plot (Figure S1) ----

# Mapping of scenario IDs to descriptions (plot facets)
scenarios = qc("spread_perennial" = "Perennial transmission\nRTS,S/AS01\nChildren", 
               "spread_seasonal" = "Seasonal transmission\nRTS,S/AS01\nChildren",
              "spread_perennial_massvax" = "Perennial transmission\nRTS,S/AS01-like\nMass vaccination", 
              "spread_seasonal_massvax" = "Seasonal transmission\nRTS,S/AS01-like\nMass vaccination",
              "spread_perennial_altvax" = "Perennial transmission\nAlternative vaccine\nChildren", 
              "spread_seasonal_altvax" = "Seasonal transmission\nAlternative vaccine\nChildren",
              "spread_perennial_altvax_massvax" = "Perennial transmission\nAlternative vaccine\nMass vaccination", 
              "spread_seasonal_altvax_massvax" = "Seasonal transmission\nAlternative vaccine\nMass vaccination")

# Build paths to respective model emulators
rds_paths = paste0("./output/", names(scenarios), 
                   "/7_gp_models/sample_0/f0001.selection_coefficient.rds")

# Auxiliary function to extract data from test and train set
format_emulator = function(x) {
  
  df = readRDS(x)
  
  test_set = data.table(actual = df$Z_test, predicted = df$Z_test_predicted, set = "Testing")
  train_set = data.table(actual = df$Z_train, predicted = df$Z_train_predicted, set = "Training")
  
  rbind(train_set, test_set)
  
}

# Extract data from test and train set for every scenario
df = lapply(rds_paths, FUN = format_emulator)

# Assign scenario ID as name of each list element
names(df) = names(scenarios)

# Concatenate results
df = rbindlist(df, idcol = "scenario")

# ---- Formatting for prettier plotting ----

# Parse training/testing as factor
df[, set := factor(set, levels = c("Training", "Testing"))]

# Define transmission setting label
df[, transmission := "Perennial transmission"]
df[grepl("seasonal", scenario), transmission := "Seasonal transmission"]

df[, transmission := factor(transmission, 
                            levels = c("Seasonal transmission", 
                                       "Perennial transmission"))]

# Define vaccine label
df[, vaccine := "RTS,S/AS01"]
df[grepl("mass", scenario), vaccine := "RTS,S/AS01-like"]
df[grepl("altvax", scenario), vaccine := "Longer-lasting vaccine"]
df[, vaccine := factor(vaccine, levels = c("RTS,S/AS01", "RTS,S/AS01-like", "Longer-lasting vaccine"))]

# Define target group label
df[, deployment := "Children only"]
df[grepl("mass", scenario), deployment := "Mass vaccination"]

df[, scenario := factor(scenario, levels = names(scenarios), labels = unname(scenarios))]

# Ensure that plot is quadratic
axis_limits = c(min(df$actual, min(df$predicted)), 
                max(df$actual, max(df$predicted)))



# Compile diagnostic emulator plot
diagnostic_plot = 
  df %>%
  ggplot(aes(x=actual, y=predicted, color=set)) +
  geom_point(alpha = .6, size = 1.5) +
  geom_abline(slope = 1, linewidth = 0.8) +
  scale_color_manual(values = pal, name = "Set") +
  scale_y_continuous(limits = axis_limits, name = "Predicted selection coefficient") +
  scale_x_continuous(limits = axis_limits, name = "Estimated selection coefficient") +
  facet_nested(transmission ~ deployment + vaccine,
               nest_line = element_line(color = "black", linewidth = 0.5)) +
  coord_fixed() +
  theme_amr(base_size = 10) +
    theme(axis.title.x = element_text(margin = margin(10, 0, -15, 0)),
          axis.title.y = element_text(margin = margin(r = 5, l = 5)),
          #axis.title.y = element_text(margin = margin(0, 5, 0, 0)),
          #legend.title = element_text(margin = margin(0, 0, -10, 0)),
          #legend.text = element_text(margin = margin(0, 0, -10, 0)),
          legend.box.margin = margin(0, 0, -10, 0),
          plot.margin = margin(t = 0, b = 0),
          strip.placement = "outside")

# Save emulator plot
ggsave(file.path(fig_output_pth, "SI_Figure-S3_Emulator-performance.png"), diagnostic_plot, width = 8, height = 5)

# ---- Produce vaccine efficacy plot ----

# Specify vaccine characteristics
vaccines = list(
  rtss = list(name = "RTS,S/AS01 vaccine",
              initial_efficacy = c(0.911, 0.7745),
              timing   = c(0, 365),
              decay = list(L = 0.61,
                           k = 0.69)),
  ll = list(name = "Longer-lasting vaccine",
            initial_efficacy = c(0.9, 0.89),
            timing = c(0, 365),
            decay = list(L = 2.56,
                         k = 0.75)))

# Auxiliary function for Weibull decay
fn_weibull = function(t, L = 1.0, k = 1.0) {
  
  # Formula for decay (see OpenMalaria wiki)
  y = exp(-(t/L)^k * log(2))
  
  return(y)
  
}

# Define timesteps (x-axis values)
timestep = 0:1095

vaccine_ids = names(vaccines)

vaccine_list = list()

# For every vaccine, ...
for(i in seq_along(vaccine_ids)) {
  
  # Pre-allocate data.table for temporal efficacy profile
  vaccine_df = data.table(id   = vaccine_ids[i],
                          time = timestep,
                          initial_efficacy = as.numeric(NA),
                          dose_interval = as.numeric(NA))
  
  # Set efficacy to initial efficacy every time dose is administered
  vaccine_df[time %in% vaccines[[i]]$timing, initial_efficacy := vaccines[[i]]$initial_efficacy]
  
  # Carry initial efficacy forward
  vaccine_df[, initial_efficacy := nafill(initial_efficacy, type = "locf")]
  
  # Number of doses
  vaccine_df[time %in% vaccines[[i]]$timing, dose_interval := seq(1, length(vaccines[[i]]$timing))]
  
  vaccine_df[, dose_interval := nafill(dose_interval, type = "locf")]
  
  # Timesteps since dose was administered
  vaccine_df[, timestep_decay := seq(0, .N - 1), by = list(id, dose_interval)]
  
  # Calculate protective efficacy after decay
  vaccine_df[, efficacy_decay := initial_efficacy*fn_weibull(timestep_decay/365, 
                                                             L = vaccines[[i]]$decay$L,
                                                             k = vaccines[[i]]$decay$k)]
  
  vaccine_list[[i]] = vaccine_df
  
}

# Extract "name" element from list of vaccines
dict_id = lapply(vaccines, function(x) { x$name }) %>%
  unlist()

# Concatenate vaccine data
full_vaccine_df = rbindlist(vaccine_list)

# Format vaccine_id for prettier plotting
full_vaccine_df[, id := factor(id, levels = names(dict_id),
                               labels = unname(dict_id))]

# Build color palette
pal = c("#22333B", "#4C8577")

pal = c("#154d6e", "#468AB2", "#b8d4e6")

# Compile vaccine efficacy plot
efficacy_plot = 
  full_vaccine_df %>%
  ggplot(aes(x=time, y=efficacy_decay, fill=id, group = id, color=id)) +
  geom_area(position = "identity", alpha = 0.6) +
  scale_y_continuous(name = "Vaccine efficacy (%)", labels = scales::percent,
                     limits = c(0, 1), expand = expansion(mult = c(0.0, 0.01))) +
  scale_x_continuous(name = "Time since third vaccine dose (in days)",
                     limits = range(timestep), expand = expansion(mult = c(0.0, 0.0)),
                     breaks = seq(0, 365*3, by = 180)) +
  scale_color_manual(name = "Vaccine", labels = dict_id, values = pal) +
  scale_fill_manual(name = "Vaccine", labels =dict_id, values = pal) +
  theme_bw() +
  theme(
    legend.title = element_text(size = 12, color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    axis.title.x = element_text(margin = margin(t = 15, b = 0)),
    axis.text = element_text(size = 12, color = "black"),
    plot.margin = margin(t = 5, r = 10, b = 5, l = 5),
    legend.position = "bottom"
  )

# Save vaccine efficacy plot
ggsave(file.path(fig_output_pth, "SI_Figure-S1_Vaccine-efficacy.png"), efficacy_plot, width = 8, height = 4)

