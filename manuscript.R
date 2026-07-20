################################################################################
# MANUSCRIPT
#                                                                             
# Create majority of manuscript figures
#
################################################################################

#-------------------------------------------------------------------------------
# Main function for creating manuscript figures
#-------------------------------------------------------------------------------
prepare_manuscript = function(o) {
  
  message("* Preparing manuscript figures")
  
  # Plotting options
  plot_seasonality = FALSE
  
  # Load simulation parameters
  p = load_parameters(o)
  
  # Shortcut to parameter grid
  params = p$parameters$fixed_grid
  
  # Load reduced set of scenarios
  resistance_scenarios = fread(paste0(o$pth$input, 
                                      "/scenarios/resistance_scenarios_lim.txt"))
  
  # If specified, plot seasonality profile
  if(plot_seasonality)
    plot_seasonality(o, p)
  
  # ---- Go through different scenario types and choose appropriate plotting function ----
  
  # Choose EIR to match certain PfPR2-10 level
  if(grepl("eir", o$analysis_name) & !grepl("match", o$analysis_name) & 
     !grepl("spread", o$analysis_name)) {
    align_eir(o, p)
    return(NULL)
  }
  
  # Treatment failure rates
  if(grepl("treat", o$analysis_name) & !grepl("lhc", o$analysis_name)) { 
    plot_treatment_failures(o, p, resistance_scenarios)
    return(NULL) 
  }
  
  # Sensitivity analysis for selection coefficient
  if(grepl("spread", o$analysis_name) | grepl("lhc", o$analysis_name)) { 
    plot_sensitivity_results(o, p)
    return(NULL)
  }
  
  # Vaccine effectiveness
  if(grepl("morbidity", o$analysis_name)) {
    plot_vaccine_effectiveness(o, p, scenarios = resistance_scenarios)
    return(NULL)
  }
  
  stop(paste0("! Plotting function for scenario '", o$analysis_name, 
              "' not defined"))
  
}

# ------------------------------------------------------------------------------
# Return grid of scenarios
# ------------------------------------------------------------------------------
list_scenarios = function(o, p, sample_id = 0, return_grid = FALSE) {
  
  all_sims = readRDS(paste0(o$pth$parameters, "/sample_", sample_id, "/all_scenarios.rds"))
  
  set(all_sims, j = c("cid", "seed", "sid", "scenario"), value = NULL)
  
  all_sims = all_sims %>% unique()
  
  if(return_grid) {
    
    output = all_sims[order(fid)]
    
  } else {
    
    output = all_sims$fid %>% sort()
    
  }
  
  return(output)
  
}

# ------------------------------------------------------------------------------
# Aggregate all relevant OpenMalaria measures
# ------------------------------------------------------------------------------
aggregate_output = function(o, p, t_start = 1, vars, summarize_age_groups = FALSE,
                            scenarios = NULL) { 
  
  parameters_df = list_scenarios(o, p, return_grid = TRUE)
  
  fids = parameters_df$fid %>% sort()
  
  # Initialize empty lists
  output_list = pfpr2_10 = eir = list()
  
  # Define age groups based on XML input
  age_groups = c(unlist(p$age_groups), 100)
  
  # Filter indices corresponding to age groups for 2-10 year-olds
  age_groups_2_10 = which(age_groups > 2.0 & age_groups <= 10.0)
  
  # Check if treatment scenario --> survey time step is 5 days
  if(grepl("treat", o$analysis_name)) {
    
    # Store time step before vaccination is introduced
    pfpr2_10_cutoff = ((p$burn_in - p$lead_time) * (365/5)) %>% round(digits = 0)
    
    # Include relevant age groups (cohort should be empty, but still include)
    age_groups_2_10 = c(age_groups_2_10, 1000 + age_groups_2_10)
    
    eir_mult = 1
    
  } else {
    
    pfpr2_10_cutoff = ((p$burn_in - p$lead_time) * p$timeline$surveys_per_year) %>% round(digits = 0)
    
    eir_mult = 6
    
  }
  
  pb = startProgressBar(length(fids))
  
  # Loop through settings
  for(i in seq_along(fids)) {
    
    # Read processed OpenMalaria output
    output_df = readRDS(paste0(o$pth$processed_out, "sample_0/", fids[i], ".rds"))
    
    # Exclude first survey
    output_df = output_df[survey > 1, ]
    output_df[, survey := survey - 1]
    
    if(nrow(output_df) == 0) {
      warning(paste0("! Results for scenario ", fids[i], " not found."))
      next
    }
    
    # ---- PfPR2-10 analysis for morbidity scenarios ----
    
    if(grepl("morbidity", o$analysis_name) || grepl("treat", o$analysis_name)) {
      
      # Calculate PfPR2-10 prevalence
      pfpr2_10[[i]] = extract_pfpr_2_10(o, p, output_df, age_groups_2_10, pfpr2_10_cutoff)
      
      # Parse seed as number
      pfpr2_10[[i]][, seed := as.numeric(lapply(strsplit(scenario, split = "_"), "[", 4))]
      
      # Reassign setting ID
      pfpr2_10[[i]][, scenario := fids[i]]
      
      eir[[i]] = output_df[indicator == "simulatedEIR" & survey <= pfpr2_10_cutoff]
      set(eir[[i]], j = "group", value = NULL)
      
      # Parse seed as number
      eir[[i]][, seed := as.numeric(lapply(strsplit(scenario, split = "_"), "[", 4))]
      
      # Reassign setting ID
      eir[[i]][, scenario := fids[i]]
      
    }
    
    # Exclude all data captured before start date
    output_df = output_df %>%
      filter(survey > t_start,
             indicator %in% vars) %>%
      as.data.table()
    
    # Collapse age groups to under five years of age (1) and five years and older (2)
    if(summarize_age_groups) {
      
      output_df = output_df %>% 
        mutate(group = case_when(
          group %in% c("1", "2") ~ "1",
          group %in% c("3", "4") ~ "2",
          
          group %in% c("1001", "1002") ~ "1001",
          group %in% c("1003", "1004") ~ "1002",
          
          
          group %in% c("0000001","0000002") ~ "0000001",
          group %in% c("0000003","0000004") ~ "0000002",
          
          group %in% c("0001001","0001002") ~ "0001001",
          group %in% c("0001003","0001004") ~ "0001002",
          
          group %in% c("1000001","1000002") ~ "1000001",
          group %in% c("1000003","1000004") ~ "1000002",
          
          group %in% c("1001001","1001002") ~ "1001001",
          group %in% c("1001003","1001004") ~ "1001002",
          
          TRUE ~ NA
          
        )) %>%
        group_by(scenario, survey, indicator, group) %>%
        summarize(value = sum(value)) %>%
        as.data.table()
      
    }
    
    # Pause if any NAs left in group column
    if(nrow(output_df[is.na(group),]) > 0)
      browser()
    
    # Parse seed as number
    output_df[, seed := as.numeric(lapply(strsplit(scenario, split = "_"), "[", 4))]
    
    # Reassign setting ID
    output_df[, scenario := fids[i]]
    
    # Calculate cumulative sum
    output_df[order(survey), cum_value := cumsum(value), by = list(scenario, indicator, group, seed)]
    
    # Store output in list format
    output_list[[i]] = output_df
    
    setTxtProgressBar(pb, i)
    
  }
  
  close(pb)
  
  # Concatenate all list output
  output_list = rbindlist(output_list)
  
  if(length(pfpr2_10) > 0) {
    
    # Breaks differ by scenario type
    if(grepl("treat", o$analysis_name)) {
      
      # ... in treatment failure scenarios, time step is 5 days
      surveys_per_year = 365/5
      
      year_breaks = seq(1, surveys_per_year*p$lead_time, by = surveys_per_year)
      
    } else if(grepl("morbidity", o$analysis_name)) {
      
      # ... in vaccine effectiveness scenarios, time step is one month
      survey_max = ((p$burn_in - p$lead_time) * p$timeline$surveys_per_year) %>% round(digits = 0)
      
      year_breaks = seq(from = 1, to = survey_max, by = round(p$timeline$surveys_per_year))
      
    } else {
      
      stop("PfPR cannot be calculated")
      
    }
    
    # Concatenate PfPR2-10 data
    pfpr2_10_df = rbindlist(pfpr2_10)
    
    # Assign year of observation
    pfpr2_10_df[, year := findInterval(survey, vec = year_breaks)]
    
    # Calculate annualized EIR
    pfpr2_10_annualized_df = 
      pfpr2_10_df[, .(pfpr2_10 = mean(pfpr2_10)), by = list(scenario, seed, year)]
    
    # Only consider period after burn-in
    pfpr2_10_annualized_df = pfpr2_10_annualized_df[year <= (p$burn_in - p$lead_time)]
    
    # Merge with scenario parameters
    pfpr2_10_annualized_df = 
      p$parameters$fixed_grid[pfpr2_10_annualized_df, on = c("fid" = "scenario")]
    
    # Aggregate PfPR2-10 values
    pfpr2_10_aggregated_df = pfpr2_10_annualized_df[, .(min = min(pfpr2_10),
                                                        lower = quantile(pfpr2_10, 0.025),
                                                        mean = mean(pfpr2_10),
                                                        median = median(pfpr2_10),
                                                        upper = quantile(pfpr2_10, 0.975),
                                                        max = max(pfpr2_10)),
                                                    by = list(eir, access)]
    
    # ---- EIR before vaccination ----
    
    # Concatenate EIR data
    eir_df = rbindlist(eir)
    
    # Assign year of observation
    eir_df[, year := findInterval(survey, vec = year_breaks)]
    
    # Calculate annualized EIR
    eir_annualized_df = eir_df[, .(sim_eir = sum(eir_mult*value)), by = list(scenario, seed, year, indicator)]
    
    # Only consider period after burn-in
    eir_annualized_df = eir_annualized_df[year <= (p$burn_in - p$lead_time)]
    
    # Merge with scenario parameters
    eir_annualized_df = 
      p$parameters$fixed_grid[eir_annualized_df, on = c("fid" = "scenario")]
    
    # Aggregate EIR values
    eir_aggregated_df = eir_annualized_df[, .(min = min(sim_eir),
                                              lower = quantile(sim_eir, 0.025),
                                              mean = mean(sim_eir),
                                              median = median(sim_eir),
                                              upper = quantile(sim_eir, 0.975),
                                              max = max(max(sim_eir))),
                                          by = list(eir, access)]
    
    # Store everying in list
    pfpr2_10_output = list(age_groups      = age_groups,
                           age_groups_2_10 = age_groups_2_10,
                           pfpr2_10_cutoff = pfpr2_10_cutoff,
                           year_breaks     = year_breaks,
                           output          = pfpr2_10_aggregated_df)
    
  } else {
    
    eir_output = NULL
    
    pfpr2_10_output = NULL
    
  }
  
  # Collate aggregated output
  aggregated_output = list(output     = output_list,
                           prevalence = pfpr2_10_output,
                           eir        = eir_aggregated_df)
  
  return(aggregated_output)
  
}

# ------------------------------------------------------------------------------
# Extract parasite prevalence in two-to-ten year-old children
# ------------------------------------------------------------------------------
extract_pfpr_2_10 = function(o, p, output_df, age_groups, t_end) {
  
  # Use nPatent to calculate PfPR2-10
  output_short_df = output_df[indicator %in% c("nPatent", "nHost") &
                                survey <= t_end &
                                group %in% age_groups, ]
  
  # Aggregate nHost and nPatent over age groups
  output_short_df = output_short_df[, .(value = sum(value)), 
                                    by = list(scenario, indicator, survey)]
  
  # Convert to long format to wide format
  output_short_df = output_short_df %>% pivot_wider(names_from = "indicator") %>%
    as.data.table()
  
  # Calculate parasite prevalence
  pfpr2_10_df = output_short_df[, "pfpr2_10" := nPatent/nHost]
  
  # Remove redundant columns
  pfpr2_10_df[, c("nPatent","nHost") := c(NULL, NULL)]
  
  return(pfpr2_10_df)
  
}

# ------------------------------------------------------------------------------
# Format values in percentage
# ------------------------------------------------------------------------------
format_pct = function(vals) {
  
  return(paste0(round(vals * 100, 0), "%"))
  
}

#-------------------------------------------------------------------------------
# Plotting function for results of sensitivity analysis
#-------------------------------------------------------------------------------
plot_sensitivity_results = function(o, p) {
  
  pal = c("#468AB2", "#BF3227", "#B5C751", "#EDCD64", "#868686", "#C5D3E5", 
          "#EBD7A3")
  
  # Load mapping table for plot file names
  mapping_titles = fread(file.path(o$pth$input, "scenarios", "plot_titles_spread.txt"))
  
  # ---- Preparation for plotting ----
  
  # Add line break to labels for Sobol indices and partner drug resistance
  p$labels$idx_order = qc(first = "First-\norder", total = "Total")
  p$labels$parameters["lf_resistance_level"] = "Resistance to\nlong-acting drug"
  
  factors = list(coverage = "Vaccination\ncoverage",
                 ar_mkr_resistant = "Resistance to\nartemisinin",
                 lf_resistance_level = "Resistance to\nlong-acting drug",
                 access = "Treatment\naccess",
                 eir = "EIR\n")
  
  factors_flat = list(coverage = "Vaccination coverage",
                      ar_mkr_resistant = "Resistance to artemisinin",
                      lf_resistance_level = "Resistance to long-acting drug",
                      access = "Treatment access",
                      eir = "EIR")
  
  # ---- Loading and preparing data ----
  
  # Path to output from sensitivity analysis
  sensitivity_path = paste0(o$pth$manuscript, "global_sensitivity_analysis.rds")
  
  # Load data or return NULL if file not available
  if(file.exists(sensitivity_path)) {
    results_df = readRDS(sensitivity_path)
  } else {
    stop("! Results for sensitivity analysis not available")
  }
  
  fields = list("median", "lower", "upper")
  
  # Shortcut to sobol indices
  sobol_df = results_df$sobol_indices
  
  # Shortcut to data from global sensitivity analysis
  sensitivity_df = results_df$sensitivity_data
  
  # List of LHC outcomes that GPs were fitted to
  outcomes = sobol_df$outcome %>% unique()
  
  # List of settings
  fids = sobol_df$fid %>% unique()
  
  # Loop through LHC outcomes and settings
  for(this_outcome in outcomes) {
    
    for(this_fid in fids) {
      
      # ---- Sensitivity analysis ----
      
      # Extract data from GSA for this scenario and outcome
      this_sensitivity_df = sensitivity_df[fid == this_fid & outcome == this_outcome]
      
      # List of values on y-axis 
      cols = unlist(fields)
      
      # Determine min and max value for values to be plotted on y-axis
      y_min = min(this_sensitivity_df[, ..cols])
      y_max = max(this_sensitivity_df[, ..cols])
      
      subplots = list()
      
      # Loop through every single parameter
      for(this_factor in seq_along(factors)) {
        
        subplot_df = this_sensitivity_df[variable == names(factors)[this_factor]] 
        
        subplots[[this_factor]] = subplot_df %>%
          ggplot(aes_string(x="x", y=fields[[1]])) +
          geom_ribbon(aes_string(ymin=fields[[2]], ymax=fields[[3]]), 
                      alpha = .3, size = 0, fill = pal[this_factor]) +
          geom_line(size = 1.5, color = pal[this_factor]) +
          scale_y_continuous(limits = c(y_min, y_max), 
                             expand = expansion(mult = c(0, 0.01)))
        
        # Based on parameter, create custom x-axis
        if(names(factors)[this_factor] == "coverage") {
          
          subplots[[this_factor]] = subplots[[this_factor]] +
            scale_x_continuous(breaks = c(0.0, 0.25, 0.5, 0.75, 1),
                               labels = scales::percent,
                               expand = expansion(mult = c(0,0)))
          
        } else if(names(factors)[this_factor] == "ar_mkr_resistant") {
          
          subplots[[this_factor]] = subplots[[this_factor]] +
            scale_x_reverse(breaks = c(max(subplot_df$x), min(subplot_df$x)), 
                            labels = c("Low", "Very high"),
                            expand = expansion(mult = c(0,0)))
          
        } else if(names(factors)[this_factor] == "lf_resistance_level") {
          
          subplots[[this_factor]] = subplots[[this_factor]] +
            scale_x_continuous(breaks = c(min(subplot_df$x), max(subplot_df$x)), 
                               labels = c("None", "Very high"),
                               expand = expansion(mult = c(0,0)))
          
        } else if(names(factors)[this_factor] == "access") {
          
          subplots[[this_factor]] = subplots[[this_factor]] +
            scale_x_continuous(breaks = c(min(subplot_df$x), 0.2412, max(subplot_df$x)), 
                               labels = c("10%", "50%", "80%"),
                               expand = expansion(mult = c(0,0)))
          
        } else if(names(factors)[this_factor] == "eir") {
          
          subplots[[this_factor]] = subplots[[this_factor]] +
            scale_x_continuous(breaks = c(5, 20, 40, 60, 80),
                               expand = expansion(mult = c(0,0)))
          
        }
        
        subplots[[this_factor]] = subplots[[this_factor]] +
          ggtitle(factors[[this_factor]]) +
          theme(axis.title = element_text(size = 16, color = "black"),
                axis.text = element_text(size = 16, color = "black"),
                strip.text = element_text(size = 16, color = "black"),
                plot.title = element_text(size = 16, hjust = 0.5),
                plot.margin = margin(10, 40, 0, 10),
                axis.title.y = element_blank(),
                axis.title.x = element_blank(),
                axis.text.x = element_text(margin = margin(10, 50, 5, 5)),
                panel.border = element_rect(color = "black", fill = NA))
        
      }
      
      # Concatenate all subplots
      subplots_all = ggarrange(plotlist = subplots)
      
      # Create blank y-axis
      synthetic_axis = ggplot(NULL) +
        geom_blank() +
        scale_y_continuous(name = "Selection coefficient") +
        theme(axis.title.y = element_text(size = 16, color="black", 
                                          margin = margin(0, 10, 0, 10)),
              panel.background = element_blank(),
              plot.background = element_blank(),
              panel.grid = element_blank())
      
      # Put together synthetic axis and subplots for each model parameter
      sensitivity_subplot = ggarrange(plotlist = list(synthetic_axis, subplots_all), 
                                      widths = c(1, 30))
      
      # ---- Sobol indices ----
      
      # Extract Sobol indices for this scenario and outcome
      this_sobol_df = sobol_df[fid == this_fid & outcome == this_outcome]
      
      # Format factor names for plotting      
      this_sobol_df[, parameter := factor(parameter, levels = names(factors_flat),
                                          labels = unname(unlist(factors_flat)))]
      
      # Put together plot for Sobol indices
      sobol_plot = this_sobol_df %>%
        ggplot(aes(x=fid, y=norm_value, fill = parameter)) +
        geom_bar(stat = "identity", position = "stack") +
        facet_wrap(~ idx_order, labeller = labeller(idx_order = p$labels$idx_order)) +
        scale_y_continuous(name = NULL, limits = c(0, NA), expand = expansion(mult = c(0.01, 0.01))) +
        scale_fill_manual(values = pal, name = "Parameter") +
        facet_wrap(~ idx_order, labeller = labeller(idx_order = p$labels$idx_order)) +
        theme_amr() +
        theme(
          axis.text = element_text(16),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          legend.margin = margin(5, 0, 0, 0),
          plot.background = element_rect(fill = "white", color = "white"),
          legend.text = element_text(margin = margin(l = 2, r = 30, unit = "pt")),
          legend.box.margin = margin(10, 0, 10, 0),
          legend.title = element_text(margin = margin(r = 20)))
      
      # Merge plots (global sensitivity analysis, Sobol indices)
      # Warnings are due to overlapping elements in plot
      manuscript_plot = suppressWarnings(ggarrange(sensitivity_subplot, sobol_plot,
                                                   labels = c("A", "B"),
                                                   common.legend = TRUE, legend = "bottom",
                                                   widths = c(4,1),
                                                   font.label = list(size = 18)) +
                                           theme(plot.background = element_rect(fill = "white", color = NA)))

      # Save plot
      if(length(fids) == 1 & this_outcome == "selection_coefficient") {
        ggsave(file.path(o$pth$publication_figures, mapping_titles[analysis_name == o$analysis_name, plot_title]), 
               manuscript_plot, width = 16, height = 8)
      } else {
        stop(paste0("! File name for this plot not specified: ", this_outcome))
      }
      
    }
    
  }
  
}

#-------------------------------------------------------------------------------
# Plotting function for vaccine effectiveness analysis
#-------------------------------------------------------------------------------
plot_vaccine_effectiveness = function(o, p, scenarios) {
  
  # Color palette
  pal = c("#468AB2", "#BF3227", "#EDCD64", "#868686", "#C5D3E5", 
          "#EBD7A3")
  
  # Path to consolidated data file
  manuscript_df_file = paste0(o$pth$manuscript, "manuscript_data.rds")
  
  if(file.exists(manuscript_df_file)) {
  
    # ... load file if it already exists  
    output_list = readRDS(manuscript_df_file)
    
  } else {
    
    # ... if not, start consolidating data
    
    # Filter for model parameters of predefined resistance scenarios
    params_reduced = p$parameters$fixed_grid %>%
      inner_join(resistance_scenarios, by = c("ar_mkr_resistant", "lf_resistance_level"))
    
    # Extract data from post-processed files
    output_list = aggregate_output(o, p, round(p$burn_in * p$timeline$surveys_per_year) + p$timeline$surveys_per_generation, 
                                   vars = c("nEPIVaccinations", "nMassVaccinations",
                                            "simulatedEIR", "nUncomp", "nSevere", "nHost", "nPatent"),
                                   scenarios = params_reduced)
    
    # Save extracted data for future use
    saveRDS(output_list, manuscript_df_file)
  } 
  
  # Parse 'group' column from OpenMalaria output
  output_df = extract_groups(output_list$output)
  
  # Selection of OpenMalaria output variables
  vars = c("nTreatments1", "nTreatments2", "nSevere", "nUncomp", "nEPIVaccinations", "nMassVaccinations", "simulatedEIR")
  
  # Filter simulation output
  output_df = output_df[indicator %in% vars, ]
  
  # Filter for first- and second-line treatment
  treatments_df = output_df[indicator %in% c("nTreatments1", "nTreatments2")]
  
  # Consolidate first- and second-line treatment into single column
  treatments_df[, indicator := "nTreatments"]
  
  treatments_df = treatments_df[, .(value = sum(value), cum_value = sum(cum_value)), by = list(scenario, survey, indicator, seed, genotype, cohort, age_group)]
  
  # Concatenate both data.tables
  output_df = rbindlist(list(output_df, treatments_df), use.names = TRUE)
  
  # Data.table with scenario parameters
  params_df = readRDS(paste0(o$pth$parameters, "sample_0/all_scenarios.rds"))
  
  # Load applicable PfPR2-10 table (based on seasonality parameter)
  if(grepl("seasonal", o$analysis_name))
    eir_aligned = readRDS(paste0(o$pth$output, "/morbidity_seasonal_eir/8_manuscript/manuscript_data_emulator.rds"))$eir_aligned[parameter == "pfpr_aggregated"]
  else
    eir_aligned = readRDS(paste0(o$pth$output, "/morbidity_perennial_eir/8_manuscript/manuscript_data_emulator.rds"))$eir_aligned[parameter == "pfpr_aggregated"]
  
  # Rename columns for easier merging
  eir_aligned = eir_aligned[, .(access, pfpr = target_value, eir = input_eir)]
  
  # Shortcut to PfPR2-10 data
  prevalence_df = output_list$prevalence$output
  
  # List of resistance scenarios that are used for plotting
  params_reduced = params_df %>%
    inner_join(scenarios, by = c("ar_mkr_resistant", "lf_resistance_level")) 
  
  # Parse EIR as strings to facilitate merging of data.tables
  params_reduced[, eir := as.character(eir)]
  eir_aligned[, eir := as.character(eir)]
  
  # Merge EIR and PfPR2-10 data
  params_reduced = params_reduced %>%
    right_join(eir_aligned, by = c("access", "eir")) 
  
  # Use uncomplicated malaria episodes only
  results_df = copy(output_list$output[indicator == "nUncomp"])
  
  # Vector with age groups
  age_groups = c(unlist(p$age_groups), 100)
  
  # Parse age group ID into actual ages (in years)
  results_df$age_group = age_groups[results_df$group %>% as.numeric()]
  
  # Only use data for children <= 5.0 years
  aggregated_df = results_df[age_group <= 5.0, .(value = sum(value), cum_value = sum(cum_value)), 
                             by = list(scenario, survey, indicator, seed)]
  
  # Only use last simulated data point (five year mark)
  aggregated_df = aggregated_df[survey == max(survey) & scenario %in% params_reduced$fid, 
                                .(scenario, indicator, seed, value = cum_value)]
  
  # Reference scenarios assume 0% vaccination coverage
  reference_scenarios = params_reduced[coverage == min(coverage), .(eir, access, pfpr, ar_mkr_resistant, lf_resistance_level, reference_scenario = fid, id)] %>%
    unique()
  
  # All other scenarios are considered impact scenarios
  impact_scenarios = params_reduced[!(fid %in% reference_scenarios$reference_scenario)]
  
  # Build reference scenario-impact scenario mapping table
  impact_scenarios = impact_scenarios %>%
    left_join(reference_scenarios, 
              by = c("eir", "access", "pfpr", "ar_mkr_resistant", "lf_resistance_level", "id")) %>%
    select(eir, access, pfpr, coverage, scenario = fid, reference_scenario, id) %>%
    unique()
  
  # Data for impact scenarios
  impact_df = aggregated_df[scenario %in% impact_scenarios$scenario]
  
  # Data for reference_scenarios
  reference_df = aggregated_df[scenario %in% impact_scenarios$reference_scenario,
                               .(reference_scenario = scenario, indicator, seed, reference_value = value)]
  
  # Merge data of reference and impact scenarios
  impact_df = impact_df %>%
    left_join(impact_scenarios, by = "scenario") %>%
    left_join(reference_df, by = c("seed", "indicator", "reference_scenario"))
  
  setDT(impact_df)
  
  # Calculate and summarize change compared to baseline
  impact_df[, change := 1 - value/reference_value]
  
  # Remove absolute values and keep relative changes only
  impact_df[, c("value", "reference_value") := .(NULL, NULL)]
  
  # ---- Prepare output list ----
  
  # Save simulated settings
  set_fids = params_reduced[, .N, by = list(eir, access, pfpr)]
  
  # Collate data for output file
  output_list = list(fids = set_fids,
                     results = impact_df)
  
  # Save plotting data for use in manuscript_vaccine_effectiveness.R
  saveRDS(output_list, file.path(o$pth$manuscript, "plot_data.rds"))
  
  # ---- Actual plotting ----
  
  # Only plotting that happens here is for morbidity_seasonal_ve scenario (Figure 2)
  
  if(o$analysis_name != "morbidity_seasonal_ve")
    stop(paste0("! Use manuscript_vaccine_effectiveness.R to plot results for scenario '", o$analysis_name, "'"))
  
  plot_df = copy(impact_df)
  
  # ---- Assign labels and sort values for plotting ----
  
  # Parasite prevalence labels
  pfpr_vals   = unique(set_fids$pfpr) %>% unique() %>% sort()
  pfpr_labels = sprintf("paste(italic(Pf)*PR[2-10], ' : %s%%')", round(100 * pfpr_vals))
  
  names(pfpr_labels) = pfpr_vals
  
  plot_df[, pfpr_label := factor(pfpr, levels = pfpr_vals)]
  
  plot_df[, id := factor(id, 
                         levels = scenarios$id,
                         labels = scenarios$name)]
  
  # Aggregate data for plotting (mean, quantiles)
  plot_df = plot_df[, .(lower = quantile(change, 0.025), upper = quantile(change, 0.975),
                        mean = mean(change)), 
                    by = list(scenario, indicator, eir, access, pfpr, coverage, reference_scenario, id, pfpr_label)]
  
  # Only use vaccination coverage < 100% and settings with low and moderate access to treatment
  plot_df = plot_df[coverage < 1.0 & access < 0.5, ]
  
  # Format coverage and access to treatment values
  coverage_vals = plot_df$coverage %>% unique() %>% sort()
  coverage_labels = scales::percent(coverage_vals)
  
  plot_df[, coverage := factor(coverage, levels = coverage_vals,
                               labels = coverage_labels)]
  
  access_vals = unique(set_fids$access) %>% unique() %>% sort()
  access_labels = paste0(c("Low", "Moderate", "High"), " access to treatment (", scales::percent(c(0.1, 0.5, 0.8)), ")")
  
  plot_df[, access := factor(access, levels = access_vals,
                             labels = access_labels)]
  
  # Build actual plot
  impact_plot = 
    plot_df %>%
    ggplot(aes(x=coverage, y=mean, fill=id)) +
    geom_bar(stat = "identity", width = 0.8, position= position_dodge(0.9, preserve = "single")) +
    geom_errorbar(aes(ymin=lower,ymax=upper), width = .2, position = position_dodge(0.9, preserve = "single")) +
    scale_fill_manual(name = "Resistance scenario", values = pal) +
    scale_x_discrete(name = "Vaccination coverage") +
    scale_y_continuous(name = "Reduction in cumulative uncomplicated cases\nin children under five years of age",
                       labels = scales::percent,
                       expand = expansion(mult = c(0, 0.05))) +
    facet_nested(. ~ access + pfpr,
                 nest_line = element_line(linewidth = 0.5, color = "black"),
                 labeller = labeller(pfpr = as_labeller(pfpr_labels, label_parsed))) +
    theme_amr(base_size = 14) +
    theme(axis.title = element_text(size = 14),
          axis.title.x = element_text(margin = margin(t = 10)),
          axis.title.y = element_text(margin = margin(r = 10, l = 5)),
          axis.text = element_text(size = 14),
          legend.box.margin = margin(t = -5),
          legend.title = element_text(margin = margin(r = 30)),
          legend.text = element_text(size = 14, margin = margin(r = 20, l = 5)),
          strip.text = element_text(size = 14),
          legend.position = "bottom",
          plot.background = element_rect(fill = "white", color = "white"),
          legend.box.background = element_rect(fill = "white", color = "white"),
          panel.border = element_rect(color = "black", linewidth = 0.8, fill = NA))
  
  # Save plot into publication figures folder
  ggsave(file.path(o$pth$publication_figures, "MANUSCRIPT_Figure-2_Vaccine-effectiveness.png"), plot = impact_plot, width = 12, height = 6)
  
}

#-------------------------------------------------------------------------------
# Parse OpenMalaria group column (cohorts, age groups)
#-------------------------------------------------------------------------------
extract_groups = function(output_df) {
  
  # Ensure each group is a seven-digit number (leading zeros)
  output_df[, group := str_pad(group, width = 7, side = "left", pad = "0")]
  
  # Auxiliary function to parse million, thousand, and single unit
  f = function(x) { list(substr(x, 1, 1) %>% as.numeric(), substr(x, 2, 4) %>% as.numeric(), substr(x, 5, 7) %>% as.numeric()) }  
  
  # Apply auxiliary function so it returns seperate columns
  output_df[ , c("genotype", "cohort", "age_group") := f(group), by = .I]
  
  # Remove source column
  output_df[, group := NULL]
  
  # Add row number 
  output_df[, idx := .I]
  
  # Store list of OpenMalaria output metrics
  metrics = output_df$indicator %>% unique()
  
  # Loop through all metrics
  for(m in metrics) {
    
    # Create subset of data
    output_df_subset = output_df[indicator == m, .(idx, indicator, genotype, cohort, age_group)]
    
    # For each metric, store unique values of genotype, cohorts, and age groups
    genotypes = output_df_subset$genotype %>% unique()
    
    cohorts = output_df_subset$cohort %>% unique()
    
    age_groups = output_df_subset$age_group %>% unique()
    
    # If only one value is available for a column, store NA as default value
    if(length(genotypes) == 1)
      output_df[output_df_subset$idx, genotype := NA]
    
    if(length(cohorts) == 1)
      output_df[output_df_subset$idx, cohort := NA]
    
    if(length(age_groups) == 1)
      output_df[output_df_subset$idx, age_group := NA]
    
  }
  
  # Remove index column
  output_df[, idx := NULL]
  
  # Return formatted output data.table
  return(output_df)
  
}

# ------------------------------------------------------------------------------
# Plot seasonality profile
# ------------------------------------------------------------------------------
plot_seasonality = function(o, p) {
  
  # Only do this when seasonality is given as vector
  if(length(p$seasonality$values) > 1) {
    
    # Parse seasonality vector as data.table and add x-axis values
    seasonality_df= data.table(t     = 1:length(p$seasonality$values),
                               value = p$seasonality$values) 
    
    # Construct plot
    seasonality_plot = seasonality_df %>%
      ggplot(aes(x=t, y=value)) +
      geom_line() +
      scale_x_continuous(breaks = 1:length(p$seasonality$values),
                         expand = expansion(mult = c(0.01,0.01)),
                         name = "Month") +
      scale_y_continuous(name = "Proportion of annual infectious\nmosquito bites per person") +
      theme_amr(base_size = 10) +
      theme(axis.title.y = element_text(margin = margin(0, 5, 0, 0)),
            axis.title.x = element_text(margin = margin(5, 0, 0, 0)))
    
    # Save plot in manuscript folder
    ggsave(file.path(o$pth$publication_figures, "SI_Figure-S2_Seasonality-profile.png"), 
           plot = seasonality_plot, width = 7, height = 3)
    
  }
  
}

#-------------------------------------------------------------------------------
# Function for plotting treatment failures
#-------------------------------------------------------------------------------
plot_treatment_failures = function(o, p, resistance_scenarios) {
  
  # Read PfPR2-10 details from separate file
  eir_aligned_df = readRDS(paste0(o$pth$output, 
                                  "/", p$eir_filter, "/8_manuscript/manuscript_data_emulator.rds"))$eir_aligned[parameter == "pfpr_aggregated"]
  
  # Extract list of scenarios
  params_df = list_scenarios(o, p, return_grid = TRUE)
  
  # Rename columns in PfPR2-10 data.table for easier joining
  eir_aligned_df = eir_aligned_df[, .(access, eir = input_eir, pfpr = target_value)]
  
  # Parse EIR as string in both data.tables
  eir_aligned_df[, eir := as.character(eir)]
  params_df[, eir := as.character(eir)]
  
  # Join parameters with PfPR2-10 data and scenarios
  params_df = params_df %>%
    left_join(eir_aligned_df, by = c("eir", "access")) %>%
    left_join(resistance_scenarios, by = c("ar_mkr_resistant", "lf_resistance_level")) %>%
    as.data.table()
  
  # List of scenarios that belong to reduced parameter set
  scenarios_reduced = params_df$fid %>% unique()
  
  # List of settings (transmission intensity x access)
  settings_df = params_df[, .(eir, access)] %>% unique()
  
  # Previously calibrated inputEIR and PfPR2-10 values
  settings_df = settings_df[eir_aligned_df, on = list(eir, access)]
  
  # List of variables to be aggregated from simulation outputs
  vars = c("nPatent", "nHost", "nInfect")  
  
  # Check whether aggregated data exists for plotting
  if(file.exists(paste0(o$pth$manuscript, "manuscript_data.rds"))) {
    
    # Load data if already aggregated earlier
    output_list = readRDS(paste0(o$pth$manuscript, "manuscript_data.rds"))
    
  } else {
    
    # Aggregate simulation output into single file
    output_list = aggregate_output(o, p, (p$burn_in*(365/5)), vars = vars,
                                   summarize_age_groups = TRUE, 
                                   scenarios = params_df)
    
    # Save data frame into separate file
    saveRDS(output_list, paste0(o$pth$manuscript, "manuscript_data.rds")) 
    
  }
  
  # Ensure only number of infected and host are in data
  treat_raw_df = output_list$output[indicator %in% vars & scenario %in% scenarios_reduced]
  
  # Disentangle group stored in OpenMalaria output file (age group, cohort, genotype)
  treat_df = extract_groups(treat_raw_df)
  
  # Choose scenario which we use to determine deployment of cohort
  deployment_fid = params_df %>%
    filter(id == "very_high", coverage == min(coverage)) %>%
    arrange(access, prevalence_resistant, pfpr) %>%
    tail(1) %>%
    select(fid) %>%
    pull()
  
  # Extract survey dates at end of observation period (after 30d)
  survey_dates = treat_df %>%
    select(-cum_value) %>%
    filter(cohort == 1, age_group == 1, seed == 1, indicator == "nHost", scenario == deployment_fid) %>%
    filter(value == 0) %>%
    pull(survey) - 1
  
  # Remove duplicates (due to multiple scenarios)
  survey_dates = unique(survey_dates) %>% sort()
  
  # ---- Calculate treatment failure rates ----
  
  treat_df = treat_df %>%
    select(-cum_value) %>%
    filter(cohort == 1, survey %in% survey_dates) %>%
    pivot_wider(names_from = "indicator") %>%
    group_by(scenario, seed, genotype, cohort, age_group) %>%
    summarize(nHost = sum(nHost), nPatent = sum(nPatent)) %>%
    mutate(failure_rate = nPatent/nHost) %>%
    as.data.table()
  
  set(treat_df, j = c("genotype", "cohort"), value = NULL)
  
  treat_df = treat_df[!is.nan(failure_rate)]
  
  # ---- Calculate reductions in treatment failure rates compared to 0% vaccination coverage ----
  
  # Use the following columns from the data.table
  cols = c("fid", "eir", "coverage", "coverage_reduced", "access", "ar_mkr_resistant", 
           "lf_resistance_level", "prevalence_resistant", "id", "pfpr")
  
  # Exclude column coverage in reference data.table (coverage = 0%)
  reference_cols = setdiff(cols, "coverage")
  
  # Use 0% coverage as reference 
  reference_params_df = params_df[coverage == 0, ..reference_cols]
  
  # Rename scenario column to reference_scenario
  setnames(reference_params_df, "fid", "reference_scenario")
  
  # List of reference scenarios
  reference_scenarios = reference_params_df[, reference_scenario]
  
  # Data.table with results from reference scenarios
  reference_df = treat_df[scenario %in% reference_scenarios, 
                          .(reference_scenario = scenario, seed, age_group, 
                            reference_failure_rate = failure_rate)]
  
  # Data.table with results from impact scenarios (coverage greater than 0%)
  impact_params_df = params_df[!(fid %in% reference_scenarios), ..cols]
  
  # Rename scenario column
  setnames(impact_params_df, "fid", "scenario")
  
  # List of impact scenarios
  impact_scenarios = impact_params_df$scenario
  
  # Data.table with results from impact_scenarios
  impact_df = treat_df[scenario %in% impact_scenarios, 
                       .(scenario, seed, age_group, failure_rate)]  
  
  # Columns that data.tables should be joined on
  join_cols = c("eir", "coverage_reduced", "access", "ar_mkr_resistant",
                "lf_resistance_level", "prevalence_resistant", "id", "pfpr")
  
  # Mapping table of impact scenarios and reference scenarios
  reference_mapping = reference_params_df[impact_params_df, 
                                          on = setNames(join_cols, join_cols)][, .(scenario, reference_scenario)]
  
  # Merge mapping table into impact scenarios
  impact_df = reference_mapping[impact_df, on = "scenario"]
  
  # Merge in values from reference scenarios
  treat_rel_df = reference_df[impact_df, on = c("reference_scenario", "seed", "age_group")]
  
  # Calculate reduction/increase in treatment failure rate
  treat_rel_df[, change := reference_failure_rate - failure_rate]
  
  # Merge scenario parameters into data.table
  treat_rel_df = impact_params_df[treat_rel_df, on = "scenario"]
  
  # ---- "Absolute" treatment failure rates (not used anymore) ----
  
  # Simply copy original data.table
  treat_abs_df = copy(treat_df)
  
  # Remove nHost and nPatent
  set(treat_abs_df, j = c("nHost", "nPatent"), value = NULL)
  
  # Join scenario parameters into results data.table
  treat_abs_df = treat_abs_df %>%
    left_join(params_df, by = c("scenario" = "fid"))
  
  # Store both relative and absolute values
  output_list = list(relative = treat_rel_df,
                     absolute = treat_abs_df)
  
  # Save output to manuscript folder
  saveRDS(output_list, file = file.path(o$pth$manuscript, "plot_data.rds"))
  
  return(NULL)

}

# ------------------------------------------------------------------------------
# Solve for EIR that results in pre-specified PfPR2-10 setting
# ------------------------------------------------------------------------------
align_eir = function(o, p) {
  
  set.seed(o$seed + 10000)
  
  params = p$parameters$fixed_grid
  
  scenarios_df = readRDS(paste0(o$pth$parameters, "sample_0/all_scenarios.rds"))
  
  fids = scenarios_df$fid %>% unique()
  
  survey_max = ((p$burn_in - p$lead_time) * p$timeline$surveys_per_year) %>% round(digits = 0)
  
  year_breaks = seq(from = 1, to = survey_max, by = round(p$timeline$surveys_per_year))
  
  if(!file.exists(paste0(o$pth$manuscript, "eir_mapping.rds"))) {
    
    pfpr_df_list = eir_df_list = list()
    
    vars = c("nPatent", "nHost", "simulatedEIR", "inputEIR")
    
    pb = startProgressBar(length(fids))
    
    age_groups = c(unlist(p$age_groups), 100)
    
    for(i in seq_along(fids)) {
      
      # ---- Prepare data for calculation ----
      
      sim_output_df = readRDS(paste0(o$pth$processed_out, "sample_0/", fids[i], ".rds"))
      
      sim_output_df = sim_output_df[survey > 1, ]
      sim_output_df[, survey := survey - 1]
      
      sim_output_df = sim_output_df[indicator %in% vars & survey <= survey_max]
      
      # ---- Parasite prevalence ----
      
      pfpr_raw_df = sim_output_df[indicator %in% c("nPatent", "nHost")]
      
      pfpr_raw_df[, age_group := age_groups[as.numeric(group)]]
      
      pfpr_raw_df = pfpr_raw_df[age_group > 2.0 & age_group <= 10.0, ]
      
      pfpr_raw_df = pfpr_raw_df[, .(value = sum(value)), by = list(scenario, survey, indicator)]
      
      pfpr_raw_df_wide = data.table::dcast(pfpr_raw_df, scenario + survey ~ indicator, value.var = "value")
      
      pfpr_raw_df_wide[, pfpr := nPatent/nHost]
      
      pfpr_df_list[[i]] = pfpr_raw_df_wide[, .(scenario, survey, pfpr)]
      
      # ---- EIR ----
      
      eir_df = sim_output_df[indicator %in% c("simulatedEIR", "inputEIR")]
      
      set(eir_df, j = "group", value = NULL)
      
      eir_df_list[[i]] = eir_df
      
      setTxtProgressBar(pb, i)
      
    }
    
    close(pb)
    
    # ---- Preparing EIR data ----
    
    eir_df = rbindlist(eir_df_list)
    
    eir_df[, year := findInterval(survey, year_breaks)]
    
    eir_df = eir_df[year <= (p$burn_in - p$lead_time)]
    
    eir_df[, value := (p$timestep/5)*value]
    
    eir_annualized_df = eir_df[, .(value = sum(value)), by = list(scenario, indicator, year)]
    
    eir_aggregated_df = eir_annualized_df[, .(value = mean(value)), by = list(scenario, indicator)]
    
    eir_aggregated_df[, c("sid", "fid", "cid", "seed") := tstrsplit(scenario, "_")]
    
    eir_aggregated_df[, seed := as.numeric(seed)]
    
    eir_annualized_df[, c("sid", "fid", "cid", "seed") := tstrsplit(scenario, "_")]
    
    eir_annualized_df[, seed := as.numeric(seed)]
    
    # ---- Aggregating parasite prevalence ----
    
    pfpr_df = rbindlist(pfpr_df_list)
    
    pfpr_df[, year := findInterval(survey, year_breaks)]
    
    pfpr_df = pfpr_df[year <= (p$burn_in - p$lead_time), ]
    
    pfpr_annualized_df = pfpr_df[, .(value = mean(pfpr)), by = list(scenario, year)]
    
    pfpr_aggregated_df = pfpr_annualized_df[, .(value = mean(value)), by = list(scenario)]
    
    pfpr_aggregated_df[, c("sid", "fid", "cid", "seed") := tstrsplit(scenario, "_")]
    
    pfpr_aggregated_df[, seed := as.numeric(seed)]
    
    pfpr_annualized_df[, c("sid", "fid", "cid", "seed") := tstrsplit(scenario, "_")]
    
    pfpr_annualized_df[, seed := as.numeric(seed)]
    
    # ---- Aggregating EIR ----
    
    output_list = list(eir_annualized  = eir_annualized_df,
                       pfpr_annualized = pfpr_annualized_df,
                       eir_aggregated = eir_aggregated_df,
                       pfpr_aggregated = pfpr_aggregated_df)
    
    saveRDS(output_list, paste0(o$pth$manuscript, "eir_mapping.rds"))
    
  } else {
    
    output_list = readRDS(paste0(o$pth$manuscript, "eir_mapping.rds"))
    
  }
  
  output_list$eir_aggregated = output_list$eir_aggregated[indicator == "simulatedEIR"]
  
  continuous_parameters = p$parameters$continuous$parameter
  
  cols = c("scenario", "sid_cid", continuous_parameters, "value")
  
  emulator = optim_results = list()
  
  target_values = list()
  target_values$eir_aggregated = c(5, 80)
  target_values$pfpr_aggregated = c(0.1, 0.2, 0.4)
  
  eir_input = seq(from = 0.1, to = 600, by = 0.1)
  
  eir_input_std = normalise_0to1(eir_input, x_min = p$parameters$continuous[parameter == "eir", lb],
                                 x_max = p$parameters$continuous[parameter == "eir", ub]) %>%
    as.matrix() %>% unname()
  
  var_names = c("eir_aggregated", "pfpr_aggregated")
  
  var_names = c("pfpr_aggregated")
  
  for(i in seq_along(fids)) { 
    
    for(name in var_names) {
      
      emu_df = output_list[[name]][fid == fids[i]] %>%
        left_join(scenarios_df, by = c("scenario", "sid", "fid", "cid", "seed")) %>%
        as.data.table()
      
      emu_df[, sid_cid := paste0(sid, "_", cid)]
      
      emu_df = emu_df[, ..cols]
      
      keys = unique(emu_df$sid_cid) %>% unique()
      
      n_test = round(o$gp_test_split * length(keys))
      
      # Store order of continuous settings
      order_keys = emu_df$sid_cid
      
      X_std = emu_df[, ..continuous_parameters] %>% 
        normalise_params(continuous_parameters, p$parameters$continuous[, .(lb, ub)])
      
      test_idx = sample.int(length(keys), n_test, replace = FALSE)
      test_keys = keys[test_idx]
      
      train_keys = keys[!(keys %in% test_keys)]
      
      emu_std_df = cbind(sid_cid = order_keys, X_std, value = emu_df$value)
      
      X_train = emu_std_df[sid_cid %in% train_keys, ..continuous_parameters] %>% 
        as.matrix() %>% unname()
      
      Z_train = emu_std_df[sid_cid %in% train_keys, value]
      
      train_data = find_reps(X_train, Z_train, rescale = FALSE, normalize = FALSE,
                             return.Zlist = FALSE)
      
      # ---- Test data ----
      
      test_data = emu_std_df[sid_cid %in% test_keys] %>% 
        group_by(sid_cid, across(all_of(continuous_parameters))) %>% 
        summarize(value = mean(value), .groups = "drop") %>% 
        as.data.table()
      
      X_test = test_data[, ..continuous_parameters] %>% as.matrix() %>% unname()
      Z_test = test_data[, value]
      
      gp_model = mleHetGP(X = list(X0 = train_data$X0, Z0 = train_data$Z0, mult = train_data$mult),
                          Z = train_data$Z, covtype = o$gp_kernel, maxit = o$gp_max_iter)
      
      gp_model = strip(gp_model)
      
      gp_model = rebuild(gp_model, robust = TRUE)
      
      Z_train_predicted = predict(gp_model, X_train)
      
      Z_test_predicted = predict(object = gp_model, X_test)
      
      results_train = data.table(set = "train", 
                                 actual = Z_train, 
                                 predicted = Z_train_predicted$mean)
      
      results_test = data.table(set = "test", 
                                actual = Z_test, 
                                predicted = Z_test_predicted$mean)
      
      results = rbindlist(list(results_train, results_test), use.names = FALSE)
      
      lims_x = min(results$actual, results$predicted)
      lims_y = max(results$actual, results$predicted)
      
      emulator[[fids[i]]][[name]] = list(gp_model    = gp_model,
                                         std_params  = X_std$std_params,
                                         X_train     = X_train,
                                         Z_train     = Z_train,
                                         Z_train_predicted = Z_train_predicted$mean,
                                         X_test      = X_test,
                                         Z_test      = Z_test,
                                         Z_test_predicted = Z_test_predicted$mean)
      
      
      results$set = factor(results$set, levels = c("train", "test"),
                           labels = c("Train", "Test"))
      
      performance_text = paste0("Correlation in test set: ", round(cor(Z_test, Z_test_predicted$mean), digits = 4))
      
      diagnostic_plot = ggplot(results, aes(x = actual, y = predicted, color = set)) +
        geom_point(size = 3, alpha = .8) + 
        scale_color_manual(name = "Dataset", values = o$pal[2:length(o$pal)]) +
        scale_x_continuous(limits = c(lims_x, lims_y), name = "Actual") +
        scale_y_continuous(limits = c(lims_x, lims_y), name = "Predicted") + 
        geom_abline(slope = 1) +
        annotate("text", label = performance_text, x = lims_y, y = lims_x, size = 4, hjust = 1) +
        theme_amr() +
        theme(legend.position = "bottom",
              axis.title.y = element_text(margin = margin(0, 15, 0, 5)),
              axis.title.x = element_text(margin = margin(15, 0, 0, 0)),
              plot.margin = margin(10, 15, 5, 5))
      
      ggsave(paste0(o$pth$manuscript, "EMULATOR_", fids[i], "_", name, ".png"), plot = diagnostic_plot, height = 10, width = 10)
      
      target_eirs = lapply(target_values[[name]], FUN = function(x, obj, input_vals) {  
        
        val_predicted = predict(obj$gp_model, input_vals)
        
        err = (x - val_predicted$mean)^2
        
        err_min = min(err)
        
        idx_err_min = which(err == err_min)
        
        list(target_value = x, 
             idx_err_min = idx_err_min, 
             error_min = err_min)
        
      }, emulator[[fids[i]]][[name]], eir_input_std)
      
      target_eirs = rbindlist(target_eirs)
      
      optim_results[[fids[i]]][[name]] = target_eirs
      
    }
    
    optim_results[[fids[i]]] = rbindlist(optim_results[[fids[i]]], idcol = "parameter")
    
  }
  
  optim_results = rbindlist(optim_results, idcol = "fid")
  
  optim_results[, input_eir := eir_input[idx_err_min]]
  
  optim_results_access = p$parameters$fixed_grid[, .(fid, access)][optim_results, on = "fid"]
  
  optim_eir = optim_results_access[, .(parameter, access, target_value, input_eir)][order(parameter, access, target_value)]
  
  emulator$eir_aligned = optim_eir
  
  saveRDS(emulator, paste0(o$pth$manuscript, "manuscript_data_emulator.rds"))
  
}
