################################################################################
# SUMMARIZE RESULTS
#
# Calculate selection coefficient and other outcome metrics
#
################################################################################

# ------------------------------------------------------------------------------
# Main function to summarize postprocessed OpenMalaria output
# ------------------------------------------------------------------------------
summarize_results = function(o, p, m) {
  
  # Skip process step if not specified in options list or not a spread scenario
  if(!(6 %in% o$do_step)) 
    return(NULL)
  
  if((grepl("treat", o$analysis_name) & 
     !grepl("lhc", o$analysis_name)) || grepl("morbidity", o$analysis_name)) 
    return(NULL)
  
  message("* Calculating summary metrics")
  
  # Load scenario definitions
  scenario_df = readRDS(file.path(h(m)$pth$parameters, "all_scenarios.rds"))
  
  # Shortcut to scenario grid
  fids = p$parameters$fixed_grid$fid
  
  # ---- Loop through processed results and calculate selection coefficient ----
  
  # Initialize progress bar
  pb = startProgressBar(length(fids))
  
  # Loop through scenarios
  for(s in seq_along(fids)) {
    
    # Read processed OpenMalaria output
    results_df = readRDS(file.path(h(m)$pth$processed_out, 
                                   paste0(fids[s], ".rds")))
    
    # List of scenarios
    scenarios = unique(results_df$scenario)
    
    # Pre-allocate list for results
    output_list = list()
    
    # Loop through list of scenarios
    for(i in scenarios) {
      
      # Shortcut to results of scenario s  
      df = results_df[scenario == i]

      # Calculate spread metrics (selection coefficient, T50)
      spread_list = calculate_spread(p, df)
      
      # Cumulative uncomplicated episodes
      cumulative_nUncomp = df[indicator == "nUncomp", sum(value)]
      
      # Redundant list of treatment failure results (not in emulator anymore)
      failure_list = NULL
      
      # Collate list of outputs for scenarios
      output_list[[i]] = rbindlist(list(spread_list,
                                        failure_list,
                                        list(scenario = i,
                                             outcome = "cumulative_nUncomp",
                                             value = cumulative_nUncomp)))
      
    }
    
    # Concatenate all results
    output_list = rbindlist(output_list)
    
    # Append scenario ID
    output_list[, fid := fids[s]]
    
    # Save summary metrics
    saveRDS(output_list, file.path(h(m)$pth$summarised, paste0(fids[s], ".rds")))
    
    # Update the progress bar
    setTxtProgressBar(pb, s)
    
  }
  
  # Close the progress bar
  close(pb)
  
}

# ------------------------------------------------------------------------------
# Estimate selection coefficient
# ------------------------------------------------------------------------------
calculate_spread = function(p, results_df) {
  
  # Store scenario ID
  scenario_id = results_df$scenario %>% unique()
  
  # Define time at which regression will start (when drug concentration of the previous drug is equal to zero + 1 parasite generations)
  t_start = results_df[indicator == "nHostDrugConcNonZero" & 
                         group == "2000001" & value == 0, survey]
  
  t_start = t_start[t_start >= (p$timeline$resistance_start + p$timeline$surveys_per_generation)][1]
  
  # Check whether perennial or seasonal setting
  if(length(p$seasonality$values) == 1) {
    
    # Define time at which regression end (2 years after the beginning of the regression)
    t_end = t_start + round(2 * p$timeline$surveys_per_year + p$timeline$surveys_per_generation)
    
    # Determine three last observations in regression horizon
    t_last = seq(t_end, t_end - 11)
    
    # Calculate cumulative cases in last three observations
    cumul_episodes = results_df[indicator == "nUncomp" & survey %in% t_last, value]
    
    # Prevalence check
    prevalence_check = results_df[indicator %in% c("nHost", "nPatent") & 
                                    survey %in% t_last, .(survey, indicator, value)]
    
    prevalence_check = prevalence_check %>%
      data.table::dcast(survey ~ indicator) 
    
    mean_pfpr = prevalence_check[, pfpr := nPatent/nHost][, pfpr] %>% mean()
    
    # If PfPR is less than 1%, throw a message to console
    if(mean_pfpr < 0.01) {
      message(paste0("\n~ Possible elimination in scenario ", scenario_id, ": Prevalence ", mean_pfpr))
    }
    
    # Only run regression if still cases/prevalent infections at end of regression horizon
    if(sum(cumul_episodes) > 0 & mean_pfpr > 0.01) {
      
      # Filter for survey data within interval [t_start, t_end]
      regression_df = results_df[survey > t_start & survey < t_end & indicator == "innoculationsByGenotype" & group == 1, .(y = value)]
      
      # Remove data for which frequency of resistant genotype is too high or too low
      regression_df = regression_df[y >= o$prevalence_limits[1] & y <= o$prevalence_limits[2]]
      
      # Add time variable for regression
      regression_df[, x := .I/(p$timeline$surveys_per_year)]
      
      # Ensure at least two observations and no NAs in vector
      if(nrow(regression_df) >= 2 & !is.na(sum(regression_df$y))) {
        
        # Linear regression of logit prevalence on time
        slope = lm(log(y/(1-y)) ~ x, regression_df)$coefficients[2] %>%
          unname()
        
        # Convert slope into selection coefficient by dividing through number of parasite generations
        selection_coefficient = slope/(365/p$parasite_generation)
        
        # Frequency of resistant genotype at the outset of the regression
        init_prevalence = regression_df[1, y]
        
        # Selection coefficient in time rather than parasite generations
        s_adjusted = selection_coefficient*(365/p$parasite_generation)
        
        # Time it takes until 50% frequency is reached
        t50 = (1/s_adjusted)*(log(0.5/(1-0.5)) - log(init_prevalence/(1-init_prevalence)))
        
        # Compile list of output
        output = list(scenario = scenario_id,
                      outcome = c("selection_coefficient", "t50"),
                      value = c(selection_coefficient, t50)) %>%
          as.data.table()
        
      } else {
        
        # Return empty list of output instead
        output = list(scenario = scenario_id, 
                      outcome = c("selection_coefficient", "t50"), 
                      value = NA) %>%
          as.data.table()
        
      }
      
      
      
    } else {
      
      # Return empty list of output instead
      output = list(scenario = scenario_id, 
                    outcome = c("selection_coefficient", "t50"), 
                    value = NA) %>%
        as.data.table()
      
    }
    
  } else {
    
    # Capture 2 years exactly to avoid seasonal bias
    
    # Define time at which regression end (2 years after the beginning of the regression) and 6 surveys on top for rolling average
    t_end = t_start + round(2 * p$timeline$surveys_per_year + p$timeline$surveys_per_generation) + 6
    
    # ---- Prevalence check ---- 
    
    # Calculate average prevalence in second year
    pfpr_df = results_df[survey %in% seq(t_end, t_end - 11) &
                 indicator %in% c("nPatent", "nHost")] %>%
      pivot_wider(names_from = "indicator")
    
    # Convert back to data.table (pivot_wider parses it as data.frame)
    setDT(pfpr_df)
    
    # Calculate parasite prevalence 
    pfpr_df[, pfpr := nPatent/nHost]
    
    # Average parasite prevalence in second year
    mean_pfpr = mean(pfpr_df$pfpr)
    
    # Filter for survey data within interval [t_start, t_end]
    raw_df = results_df[survey > t_start & survey < t_end & indicator == "innoculationsByGenotype" &
                          group == 1, .(y = value)]
    
    # Add time variable for regression
    raw_df[, x := .I/(p$timeline$surveys_per_year)]
    
    # Pre-allocate data.table for results
    ma_df = data.table(x = numeric(18), y = numeric(18))
    
    # Calculate 6-month rolling average
    for(i in (1 + 6):24) {
      
      ma_y = mean(raw_df[(i-6):(i+6), y])
      
      ma_df[i - 6, y := ma_y]
      
      ma_df[i - 6, x := raw_df[i, x]]
      
    }
    
    cropped_y = ma_df[1:12, y]
    
    # Only take values that fall within prevalence range
    cropped_y = cropped_y[cropped_y >= o$prevalence_limits[1] & cropped_y <= o$prevalence_limits[2]]
    
    # If not enough observations are available or prevalence too low, return empty output list
    if(length(cropped_y) == 0 || mean_pfpr < 0.01) {
      
      output = list(scenario = scenario_id, 
                    outcome = c("selection_coefficient", "t50"), 
                    value = NA) %>%
        as.data.table()
      
    } else {
      
      # Enough observations are available and parasite prevalence is sufficiently high, continue with regression
      regression_df = data.table(x = ma_df[1:length(cropped_y), x], y = cropped_y)
      
      # Ensure at least two observations and no NAs in vector
      if(nrow(regression_df) >= 2 & !is.na(sum(regression_df$y))) {
        
        # Linear regression of logit prevalence on time
        slope = lm(log(y/(1-y)) ~ x, regression_df)$coefficients[2] %>%
          unname()
        
        # Convert slope into selection coefficient by dividing through number of parasite generations
        selection_coefficient = slope/(365/p$parasite_generation)
        
        # Store initial prevalence of resistant genotype
        init_prevalence = regression_df[1, y]
        
        # Transform selection coefficient back into time rather than parasite generations
        s_adjusted = selection_coefficient*(365/p$parasite_generation)
        
        # Calculate time it takes until frequency reaches 50%
        t50 = (1/s_adjusted)*(log(0.5/(1-0.5)) - log(init_prevalence/(1-init_prevalence)))
        
        # Return list of outputs
        output = list(scenario = scenario_id,
                      outcome = c("selection_coefficient", "t50"),
                      value = c(selection_coefficient, t50)) %>%
          as.data.table()
        
      } else {
        
        # Either condition is violated, so return empty list of outputs instead
        output = list(scenario = scenario_id, 
                      outcome = c("selection_coefficient", "t50"), 
                      value = NA) %>%
          as.data.table()
        
      }
      
    }
    
  }
  
  return(output)
  
}

#-------------------------------------------------------------------------------
# Function for calculating treatment failures if treatment failure rate is 
# LHC outcome (discontinued)
#-------------------------------------------------------------------------------
calculate_treatment_failures = function(o, results_df) {
  
  # List of variables to be aggregated from simulation outputs
  vars = c("nInfect", "nHost", "nPatent")
  
  # Ensure only number of infected and host are in data
  treat_df = results_df[indicator %in% vars]
  
  # Disentangle group stored in OpenMalaria output file (age group, cohort, genotype)
  treat_df = extract_groups(treat_df)
  
  # Extract survey dates at end of observation period (after 30d) and exclude first observation
  survey_dates = treat_df %>%
    filter(cohort == 1, age_group == 1, indicator == "nHost") %>%
    filter(value == 0, survey != min(results_df$survey)) %>%
    pull(survey) - 1
  
  # Remove duplicates (due to multiple scenarios)
  survey_dates = unique(survey_dates)
  
  # ---- Calculate treatment failure rates ----
  
  treat_df = treat_df %>%
    filter(cohort == 1, survey %in% survey_dates) %>%
    pivot_wider(names_from = "indicator") %>%
    group_by(scenario, genotype, cohort, age_group) %>%
    summarize(nHost = sum(nHost), 
              nInfect = sum(nInfect),
              nPatent = sum(nPatent)) %>%
    mutate(failure_nInfect = nInfect/nHost,
           failure_nPatent = nPatent/nHost) %>%
    pivot_longer(cols = c("failure_nInfect", "failure_nPatent")) %>%
    as.data.table()
  
    treat_df = treat_df[!(is.na(value) | is.na(value))]
    
    output_df = treat_df %>%
      mutate(outcome = paste0(name, "_age_", age_group)) %>%
      select(scenario, outcome, value)
    
    return(output_df)
  
}