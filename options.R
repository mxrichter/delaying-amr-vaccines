################################################################################
# OPTIONS
# 
# Set options for analyses and override with separate file if necessary
#                                                                          
################################################################################

# ------------------------------------------------------------------------------
# Set options for simulations
# ------------------------------------------------------------------------------
set_options = function(quiet = FALSE) {

  message("* Setting options and importing parameters")
  
  # Initialize list of options (o) and parameters (p)
  o = list()
  
  # ---- General options ----

  # Starting seed
  o$seed = 12345 # For sensitivity analysis, seed is set to 5*o$seed
  
  # Analysis name (must be defined as scenario in scenarios.yaml)
  o$analysis_name = "spread_seasonal"
  
  # Store user name (useful for scicore operations)
  o$user = Sys.info()[["user"]]

  # Number of different seeds to simulate
  o$n_seeds = 10
  
  # Force simulations even if already done
  o$force_simulate = FALSE
  
  # Force summarizing (= postprocessing, summarizing, fitting emulator)
  o$force_summarize = FALSE
  
  # Custom palette for plotting
  o$pal = c("#127475", "#553E4E", "#FF9B42", "#C0C999", "#EE6C4D", "#B8DBD9", "#A499B3")
  
  o$pal_seq = list(c("#C6B5CF", "#2B2031"), c("#8CCFC9", "#122B29"), c("#A9C7B9", "#18251F"))
  
  o$pal_div = c("#22333B", "#E08E45", "#904E55")
  
  o$quantiles = c(0.25, 0.75)
  
  # --- Latin hypercube (LHC) sampling ----
  
  # Number of LHC samples to generate initially (prior to any GP fitting)
  o$lhc_samples = 1000
  
  # Outcome metrics for emulator
  o$lhc_outcomes = c("selection_coefficient", "t50")
  
  # ---- Gaussian process options ----
  
  # Select GP kernel function
  o$gp_kernel = "Matern5_2" # see hetGP documentation for options
  
  # Proportion of data withheld for testing GP performance
  o$gp_test_split = 0.2
  
  # Maximum number of iterations of optimisation algorithm
  o$gp_max_iter = 1000
  
  # Use python (scikit-learn) or R (hetGP) for Gaussian processes?
  o$gp_engine = "R"
  
  # ---- Adaptive sampling ----
  
  # Run adaptive sampling to improve emulator accuracy
  o$adaptive_sampling = FALSE
  
  # Maximum number of adaptive sampling attempts
  o$sampling_max_iter = 3
  
  # Number of points to resimulate per fixed parameter set
  o$n_adaptive_samples = 50
  
  # ---- OpenMalaria options ----
  
  # Column headers for OM output
  o$om_output_headers = c("survey", "group", "indicator", "value")
  
  # ---- Cluster options ----
  
  # Choose cluster partition for parallel jobs
  o$cluster_partition = "scicore"
  
  o$job_time  = "01:30:00"  # Use "HH:MM:SS" or "D-HH:MM:SS" format
  o$job_queue = "6hours"  # Queue to use (check out scicore wiki if unfamiliar)
  
  o$batch_size = 32
  
  o$cpus_per_task = o$batch_size
  
  # Memory allocated to each single job
  o$job_memory = "15GB"
  
  # Maximum number of jobs that can be run at the same time
  o$job_limit = 500
  
  # Define names for cluster log and error files
  o$log_file = "scicore_log.txt"
  o$err_file = "scicore_error.txt"
  
  # ---- Global sensitivity analysis ----
  
  # Number of LHC samples for calculation of Sobol indices
  o$sobol_n_samples = 50000
  
  # Number of Sobol-Jansen bootstraps
  o$sobol_n_bootstraps = 5000
  
  # Plotting options (two sets of measures available)
  o$sensitivity_metrics = list(mean = c("mean", "min", "max"),
                               median = c("median", "lower", "upper"))
  
  # Default set of parameters for plotting
  o$sensitivity_statistic = "mean"
  
  # Number of bins for plotting sensitivity
  o$n_bins = 300
  
  # ---- Selection coefficient ----
  
  # Lower and upper limit for prevalence in calculation of selection coefficient
  o$prevalence_limits = c(0.3, 0.9)
  
  # ---- Other options ----
  
  # Which attributes are replaced in XML file apart from the setting variables?
  o$other_sed_vars = c("population_size", "init_date", "start_date", "end_date", 
                       "deployment_date", "timestep", "ar_mkr_sensitive", 
                       "lf_sensitive_level")
  
  # Override options if specified in my_options.yaml ----
  o = override_options(o)

  # ---- Set up folder structure ----
  
  # Take care of folder tree
  o = set_directories(o) # see directories.R
  
  # ---- Scenario data ----
  
  # Load and parse scenarios.yaml --> We need the name of the XML template
  y = parse_scenarios_yaml(o)
  
  # Path to XML template of scenario
  o$pth$xml_base = file.path(o$pth$input, paste0(y$xml_template, ".xml"))

  return(o)

}

# ------------------------------------------------------------------------------
# Override options with input from my_options.yaml
# ------------------------------------------------------------------------------
override_options = function(o) {

  # Load and parse yaml file with custom options
  o_override = read_yaml("./config/my_options.yaml")
  
  # Store labels of options set in my_options.yaml
  opts = names(o_override)
  
  # Replace standard attribute
  for(i in opts) {
    
    o[[i]] = o_override[[i]]
    
  }
  
  # Return full options list
  return(o)
  
}

# ------------------------------------------------------------------------------
# Parse scenario data from scenarios.yaml
# ------------------------------------------------------------------------------
parse_scenarios_yaml = function(o) {
  
  # Read and parse scenarios yaml file
  y = read_yaml(file.path(o$pth$config, "scenarios.yaml"))
  
  # ---- Labels for plotting ----
  
  y$labels$parameters = y$labels$parameters %>% unlist()
  
  # ---- Scenario-specific YAML ----
  
  scenarios = lapply(y$scenarios, "[[", "id") %>% unlist()
  
  names(y$scenarios) = scenarios
  
  if(!(o$analysis_name %in% scenarios))
    stop(paste0("! No parameterization found for analysis ", o$analysis_name))
  
  scenario_yaml = y$scenarios[[o$analysis_name]]
  
  y$scenarios = NULL

  # ---- Fix age group format ---
  # y$age_groups = y$age_groups %>% unlist()
  
  # scenario_yaml$age_groups = scenario_yaml$age_groups %>% unlist()
  
  # y = update_list(y, scenario_yaml)
  
  # Overwrite attributes defined in main text by scenario attributes
  
  for(n in names(scenario_yaml)) {
    
    y[n] = scenario_yaml[n]
    
  }
  
  # Remove separate scenario attributes
  y$scenarios = NULL
  
  # ---- Scenario parameters ----
  
  # Store type of parameters (usually fixed and continuous) 
  parameter_type = names(scenario_yaml$parameters)
  
  # Lookup continuous parameters and parse lower and upper bounds into data.table
  y$parameters$continuous = lapply(scenario_yaml$parameters$continuous, "[[", "array_lhc") %>%
    rbindlist(idcol = "parameter")
  
  # ---- Seasonality profiles (to be updated once seasonal setting is included) ----
  
  # Store copy of seasonality profile input (source to be replaced later)
  seasonality_profiles = y$seasonality_profiles
  
  # IDs provided for seasons
  seasonality_ids = lapply(seasonality_profiles, "[[", "id") %>% unlist()
  
  if(!(y$seasonality %in% seasonality_ids))
    stop(paste0("! Seasonality profile '", y$seasonality, "' not defined."))
  
  names(seasonality_profiles) = seasonality_ids
  
  seasonality = seasonality_profiles[[y$seasonality]]
  
  y$seasonality = seasonality
  
  # --- Vaccination deployment ----
  
  vaccination_strategies = y$vaccination_strategies
  
  names(vaccination_strategies) = lapply(vaccination_strategies, "[[", "id") %>% unlist()
  
  if(!(y$vaccine_deployment %in% names(vaccination_strategies)))
     stop(paste0("! Vaccination strategy '", y$vaccine_deployment, "' not defined."))
  
  vaccine_deployment = vaccination_strategies[[y$vaccine_deployment]]
  
  y$vaccine_deployment = vaccine_deployment
  
  y$vaccination_strategies = NULL
  
  # Store fixed parameters that define setting
  y$parameters$fixed = scenario_yaml$parameters$fixed
  
  # ---- Timing parameters ----
  
  # Number of surveys per parasite generation
  surveys_per_generation = y$parasite_generation / y$timestep
  
  # Number of surveys per year
  surveys_per_year = 365 / y$timestep
  
  # Time step at which drug resistance is introduced
  resistance_start = round(y$burn_in * surveys_per_year)
  
  y$timeline = list(surveys_per_generation = surveys_per_generation,
                    surveys_per_year = surveys_per_year,
                    resistance_start = resistance_start)
  
  # ---- Monitoring options ----
  
  # Concatenate list format of monitoring options into data.table
  y$monitoring = rbindlist(y$monitoring)

  # ---- Other options ----
  
  y$start_date = format(ymd(y$init_date) + years(y$burn_in), "%Y-%m-%d")
  
  y$end_date = format(ymd(y$init_date) + years(y$sim_duration), "%Y-%m-%d")
  
  y$deployment_date = format(ymd(y$start_date) - years(y$lead_time), "%Y-%m-%d")
  
  return(y)
  
}