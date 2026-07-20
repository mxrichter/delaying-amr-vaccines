################################################################################
# CLUSTER JOBS
#
# Perform a job from a cluster array submission.
#
################################################################################

# ------------------------------------------------------------------------------
# Perform cluster job of choice (see job_type input)
# ------------------------------------------------------------------------------
run_cluster_job = function(o, job_type, task_id, sample_id = NULL) {
  
  # ---- Processing OpenMalaria simulations ----
  
  if ((grepl("postprocess", job_type))) {
    
    # Read in scenario parameters
    p = load_parameters(o)
    
    parameters_df = readRDS(paste0(o$pth$parameters, "sample_", sample_id, "/all_scenarios.rds"))
    
    fids = parameters_df$fid %>% unique() %>% sort()
    
    process_fid = fids[task_id]
    
    # if(!file.exists(paste0(o$pth$processed_out, "sample_", sample_id, "/", process_fid, ".rds"))) {
      
      # Only use simulations which belong to the scenario
      all_sims = parameters_df[fid == process_fid, scenario]
      
      # Split list of simulations into blocks (otherwise, risk of memory issues)
      scenario_array = split(all_sims, 1:ceil(length(all_sims)/20))
      
      # Preallocate list for full results table
      all_output_list = list()
      
      # Looping through blocks of scenarios
      for(a in 1:length(scenario_array)) {
        
        # Shortcut to simulation IDs
        scenarios = scenario_array[[a]]
        
        # Construct paths and file names to output files
        scenario_files = paste0(o$pth$raw_out, "sample_", sample_id, "/", scenarios, "_out.txt")
        
        # Preallocate list for results of scenario blocks
        output_list = list()
        
        # ---- Load and format simulation results ----
        
        # Iterate through individual simulation results
        for (s in seq_along(scenario_files)) {
          
          # Check if results file exists, otherwise throw warning message
          if(file.exists(scenario_files[s])) {
            output_df = fread(scenario_files[s]) 
          } else {
            message(paste0("! Scenario ", scenario_files[s], " not found"))
            next
          }
          
          # Assign column names to make output more readable
          colnames(output_df) = o$om_output_headers
          
          # Replace indicator IDs by names and format group to seven-character strings
          output_df = output_df %>%
            left_join(p$monitoring, by = c("indicator" = "id")) %>%
            select(survey, indicator = name, group, value) %>%
            mutate(group = as.character(group)) %>%
            mutate(group = ifelse(indicator %in% c("innoculationsPerAgeGroup", "nInfectByGenotype", "nPatentByGenotype"), str_pad(group, 7, pad = "0"), group)) %>%
            select(survey, indicator, group, value)
          
          # Remove first data point (captures everything which happened until first survey)
          output_df = output_df[survey > 1]
          output_df[, survey := survey - 1]
          
          if(grepl("treat_", o$analysis_name) & grepl("_eir", o$analysis_name)) {
            
            output_df = output_df[survey <= (p$burn_in - p$lead_time)*73]
            
            output_df = output_df[indicator %in% c("nHost", "nInfect", "nPatent",
                                                   "inputEIR", "simulatedEIR")]
            
          }
          
          # ---- Calculate prevalence of resistant and sensitive genotype ----
          
          # Deep copy of results table with relevant output variables
          pct_df = copy(output_df)[indicator %in% c("innoculationsPerAgeGroup", "nInfectByGenotype", "nPatentByGenotype")]
          
          # Extract genotype and age group based on single column "group"
          pct_df[, genotype := substr(group, 1, 1)]
          pct_df[, age_group := substring(group, 2)]
          
          # Build new indicator innoculationsByGenotype
          pct_df = pct_df %>%
            filter(indicator == "innoculationsPerAgeGroup") %>%
            group_by(survey, indicator, genotype) %>%
            summarize(value = sum(value), .groups = "drop") %>%
            mutate(indicator = "innoculationsByGenotype") %>%
            select(survey, indicator, group = genotype, value)
          
          # Calculate prevalence based on absolute values and set initial prevalence to 0
          pct_df = pct_df %>% 
            group_by(survey, indicator) %>%
            mutate(value = value/sum(value)) %>%
            mutate(value = ifelse(survey == 1, 0, value))
          
          # Concatenate simulation outcomes
          output_list[[scenarios[s]]] = rbindlist(list(output_df, pct_df))
          
        }
        
        message(paste0("    > Array job ", a, "/", length(scenario_array), " finished."))
        
        # Concatenate simulation outcomes per scenario block
        all_output_list[[a]] = rbindlist(output_list, idcol = "scenario")
        
      }
      
      # Concatenate output from all scenarios across blocks
      all_output_list = rbindlist(all_output_list)
      
      # Store in separate folder
      saveRDS(all_output_list, paste0(o$pth$processed_out, "sample_", sample_id, "/", process_fid, ".rds"))
      
    # } 
    
  }
  
}
