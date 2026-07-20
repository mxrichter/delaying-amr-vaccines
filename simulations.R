################################################################################
# SIMULATIONS
# 
# Core function for running OpenMalaria simulations and fitting Gaussian 
# Process (GP) emulator
#                                                                          
################################################################################

#-------------------------------------------------------------------------------
# Main function running OpenMalaria simulations and fitting GP emulator
#-------------------------------------------------------------------------------
run_simulations = function(o) {
  
  message(paste0("* Running analysis: ", o$analysis_name))
  
  # Which steps of the workflow should be run
  o$do_step = c(1:7)
  
  # Set seed
  set.seed(o$seed)
  
  # Initialize "running variables" necessary for adaptive sampling
  m = NULL; i_adaptive = 1
  
  # Load parameters that are not stochastic
  p = load_parameters(o)

  # Start whole workflow
  while(i_adaptive <= o$sampling_max_iter) {
    
    m = adjust_input(o, p, m)
    
    # Print iteration to console if using adaptive sampling
    if(o$adaptive_sampling)
      message(paste0("  > Adaptive sampling round: ", m$sample_num))
    
    # Iteration assumed to be completed if parameter set already exists
    if(!check_existing_simulations(m)) {
      
      # Step 1) Sample parameter space for simulation
      sample_parameters(o, p, m) # see parse_parameters.R
      
      # Step 2) Build sed patterns as preliminary step for XML files
      n_jobs_unique = build_sed_patterns(o, p, m) # see openmalaria_setup.R

      # Step 3) Generate XML files for OpenMalaria simulations
      generate_xml_files(o, p, n_jobs = n_jobs_unique, m) # see openmalaria_setup.R
      
    }
    
    # Step 4) Run OpenMalaria simulations
    run_openmalaria(o, p, m) # see openmalaria_run.R
    
    # Step 5) Post-processing 
    postprocess_results(o, p, m) # see postprocess.R
      
    # Step 6) Summarized the results
    summarize_results(o, p, m) # see summary_results.R
      
    # Step 7) Fit the Gaussian process
    run_gp(o, p, m) # see GP.R
      
    # Exit while loop if we do not want to use adaptive sampling
    if(!o$adaptive_sampling)
      break
      
    # 8) Start adaptive sampling
    m = adaptive_sampling(o, m) # see adaptive_sampling.R
    
    # Increment iteration of adaptive sampling
    i_adaptive = i_adaptive + 1
    
  }
  
}

#-------------------------------------------------------------------------------
# Load all parameters that are not stochastic or depend on adaptive sampling
#-------------------------------------------------------------------------------
load_parameters = function(o) {
  
  # Load and parse scenarios.yaml
  y = parse_scenarios_yaml(o)
  
  # Remove snippets element (contains mostly monitoring variables, which are loaded separately)
  y$snippets = NULL
  
  # Remove seasonality profile (irrelevant for now)
  y$seasonality_profiles = NULL
  
  # ---- Fixed parameters (i.e., settings) ----
  
  if(is.list(y$parameters$fixed))
    y$parameters$fixed$eir = unlist(y$parameters$fixed$eir)
  
  # Merge fixed parameters into one data.table
  param_fixed_df = do.call(CJ, y$parameters$fixed)
  
  # Append identifier for set of fixed parameters
  if(grepl("treat", o$analysis_name)) {
    
    param_fixed_df[, fid := paste0("f", str_pad(.I, 5, pad = "0"))]  
    
  } else {
    
    param_fixed_df[, fid := paste0("f", str_pad(.I, 4, pad = "0"))]
    
  }
  
  # Necessary for join-operations later in the script
  setkey(param_fixed_df, fid)
  
  # Store everything in parameter list
  y$parameters$fixed_grid = param_fixed_df

  return(y)
  
}

#-------------------------------------------------------------------------------
# Create additional directories and options table in adaptive sampling process
#-------------------------------------------------------------------------------
adjust_input = function(o, p, m = NULL, overwrite = NULL) {

  if(!is.null(overwrite)) 
    m = list(sample_num = overwrite)
  
  # Initialize at first round of adaptive sampling (0), else increment by +1
  if(is.null(m)) 
    m = list(sample_num = 0)
  else 
    m$sample_num = m$sample_num + 1
  
  # Concatenate name of sub-folder for current iteration of adaptive sampling
  folder_snippet = paste("sample_", m$sample_num)
  
  # List of directories in output folder which should be amended
  sample_path_ids = c("parameters", "sed_patterns", "sim_xml", "raw_out", 
                      "processed_out", "summarised", "gp_models")
  
  # Add updated paths to sample variable
  m[[paste0("sample_", m$sample_num)]]$pth = lapply(sample_path_ids, 
                       FUN = function(s, paths, sample_num) 
                         paste0(paths[[s]], paste0("sample_", sample_num)), 
                       o$pth, m$sample_num)
  
  # Assign previously used list names
  names(m[[paste0("sample_", m$sample_num)]]$pth) = sample_path_ids
  
  # lapply creates output which is not needed here; no further need for variable
  lapply_output = lapply(m[[paste0("sample_", m$sample_num)]]$pth, 
         FUN = function(path) if(!dir.exists(path)) dir.create(path)) %>%
    suppressMessages()

  return(m)
    
}

#-------------------------------------------------------------------------------
# Simple handle to parameters of adaptive sampling round i
#-------------------------------------------------------------------------------
h = function(m, i = NULL) {
  
  if(is.null(i))
    return(m[[paste0("sample_", m$sample_num)]])
  else
    return(m[[paste0("sample_", i)]])
  
} 

#-------------------------------------------------------------------------------
# Check whether iteration was already simulated 
#-------------------------------------------------------------------------------
check_existing_simulations = function(m) {
  
  return(file.exists(file.path(h(m)$pth$parameters, "all_scenarios.rds"))) 
  
} 

#-------------------------------------------------------------------------------
# Check whether iteration has gone through a process step already; return
# setting if not yet processed
#-------------------------------------------------------------------------------
check_existing = function(o, fids, folder = NULL, m = NULL, sample_id = NULL) {
  
  # If m is NULL, check_existing called from cluster job
  if(is.null(m)) {
    
    file_path = paste0(o$pth$processed_out, "sample_", sample_id, "/", fids, ".rds")
    
    if(file.exists(file_path))
      return(TRUE)
    else
      return(FALSE)
    
  } else {
    
    file_paths = paste0(h(m)$pth[[folder]], "/", fids, ".rds")
    
    flag_exists = file.exists(file_paths)
    
    return(fids[!flag_exists])
    
  }
  
}