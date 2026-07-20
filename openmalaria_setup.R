################################################################################
# OPENMALARIA SETUP
#
# Prepare text files as input for OpenMalaria input XML files
################################################################################

# ------------------------------------------------------------------------------
# Generate sed patterns as XML input
# ------------------------------------------------------------------------------
build_sed_patterns = function(o, p, m) {

  # Skip process step if not specified in list of options
  if(!(2 %in% o$do_step)) return(NULL)
  
  message("  - Generating sed patterns for XML files")
  
  # Load parameter data.table
  param_df = readRDS(paste0(h(m)$pth$parameters, "/", "all_scenarios.rds"))
  
  # List of parameter names that should be replaced
  param_ids = c(names(p$parameters$fixed), p$parameters$continuous$parameter, 
                "seed", "prevalence_sensitive", "access_complement")
  
  # Frequency of sensitive genotype and access complement implicitly defined
  param_df[, prevalence_sensitive := 1.0 - prevalence_resistant]
  param_df[, access_complement := 1.0 - round(access, 10)]
  
  # Initiate progress bar
  pb = startProgressBar(nrow(param_df))
  
  # Loop through all different simulations
  for(i in 1:nrow(param_df)) {
    
    # Convert all parameter columns to long format
    param_df_wide = param_df[i, ..param_ids] %>%
      pivot_longer(cols = param_ids)
    
    # Add additional patterns which are neither fixed nor continuous
    additional_params = list(name = o$other_sed_vars, 
                             value = p[o$other_sed_vars] %>% unlist() %>% unname())
    
    param_df_wide = rbind(param_df_wide, additional_params)
    
    # ---- Add parameter by parameter to pattern ----
    
    seasonality = list(name = "seasonality",
                       value = extract_seasonality(o, p))
    
    age_groups = list(name = "agegroup_groups",
                      value = extract_age_groups(o, p))
    
    cohort_list = extract_cohort(o, p)
    
    cohort_monitoring = list(name  = "cohort_monitoring",
                             value = cohort_list$monitoring)
    
    cohort_removal = list(name  = "cohort_removal",
                          value = cohort_list$removal)
    
    cohort_deployment = list(name  = "cohort_deployment",
                          value = cohort_list$deployment)
    
    cohort_definition = list(name  = "cohort_definition",
                             value = cohort_list$definition)
    
    vaccine_decay = list(name = "vaccine_decay",
                            value = extract_vaccine_decay(o, p))
    
    vaccine_efficacy_sensitive = list(name = "vaccine_efficacy_sensitive",
                            value = extract_vaccine_efficacy(o, p, "sensitive"))
    
    vaccine_efficacy_resistant = list(name = "vaccine_efficacy_resistant",
                            value = extract_vaccine_efficacy(o, p, "resistant"))
    
    vaccine_deployment = list(name = "vaccine_deployment",
                      value = extract_deployment(o, p))
    
    param_df_wide = rbind(cohort_monitoring, cohort_definition, vaccine_deployment, 
                          cohort_deployment, param_df_wide, seasonality, age_groups, 
                          vaccine_decay, vaccine_efficacy_sensitive, vaccine_efficacy_resistant, cohort_removal)
    
    # Paste name-value pairs in sed format
    sed_pattern = paste0("s$@", param_df_wide$name, "@$", param_df_wide$value, "$g;")
    
    # Store the sed pattern in file
    write.table(sed_pattern, file = paste0(h(m)$pth$sed_patterns, "/", param_df[i, scenario], ".txt"), 
      row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    setTxtProgressBar(pb, i)
    
  }
  
  # Close progress bar
  close(pb)
  
  return(nrow(param_df))
  
}

# ------------------------------------------------------------------------------
# Generate XML files based on sed patterns (run on the cluster)
# ------------------------------------------------------------------------------
generate_xml_files = function(o, p, n_jobs, m) {

  # Skip step if not specified in list of options
  if(!(3 %in% o$do_step)) return(NULL)
  
  message("  - Generating XML files for OpenMalaria")
  
  # Create a new log file for the cluster jobs
  log_file = file.path(o$pth$log_files, "scicore_log.txt")
  create_bash_log(log_file)

  # Add country sub directory to file paths
  sed_country = paste0(h(m)$pth$sed_patterns, "/")
  xml_country = paste0(h(m)$pth$sim_xml, "/")
  
  # Create sequence of XMLs that are run at the same time
  n_seq = min(100, n_jobs)
  n_submit = ceil(n_jobs / n_seq)

  # Construct slurm array command for running in parallel
  slurm_array = paste0("--array=1-", n_submit, "%", min(n_submit, 300))
  
  slurm_array = paste(slurm_array, paste0("--output=", o$pth$log_files, "/scicore_%A_%a.txt"))

  # Concatenate system command
  sys_command = paste(
    "sbatch", slurm_array, "bash_make_xmls.sh",
    n_seq, o$pth$xml_base, sed_country, xml_country, log_file)

  # Invoke this command
  system(sys_command)

  # Wait for all cluster jobs to complete
  wait_for_jobs(o, log_file, n_submit)
  
}

# ------------------------------------------------------------------------------
# Generate XML snippet for seasonality pattern
# ------------------------------------------------------------------------------
extract_seasonality = function(o, p) {
  
  # Load XML snippet for seasonality section
  xml_template = readr::read_file(paste0(o$pth$input, "/snippets/seasonality_", 
                                         p$seasonality$pattern, ".txt"))
  
  # Parse numbers as strings (required by str_replace_all)
  vals = as.character(p$seasonality$values)
  
  # If only one number is provided, profile is interpreted as perennial
  if(length(vals) == 1)
    vals = rep(vals, times = 12)
  
  # Build sed pattern
  names(vals) = paste0("@seasonality", 1:12, "@")  
  
  # Replace values in the XML snippet using sed patterns
  snippet = str_replace_all(xml_template, vals)
  
  # Fix line breaks so that replacing snippets on the cluster works well
  snippet = gsub("\\n", "\\\\n", snippet)
  
  # Return complete seasonality snippet
  return(snippet)
  
}

# ------------------------------------------------------------------------------
# Generate XML snippet for vaccine efficacy
# ------------------------------------------------------------------------------
extract_vaccine_efficacy = function(o, p, genotype = "sensitive") {
  
  # Load XML snippet for vaccine efficacy section
  xml_template = readr::read_file(paste0(o$pth$input, 
                                         "/snippets/vaccine_efficacy.txt"))
  
  if(genotype == "sensitive") {
    
    # Parse numbers as strings (required by str_replace_all)
    vals = (p$vaccine$efficacy/p$vaccine$resistance_scaler) %>% 
      as.character()
    
  } else {
    
    vals = p$vaccine$efficacy %>% as.character()
    
  }
  
  # Replace values in the XML snippet using sed patterns
  snippet = str_replace(xml_template, "@vaccine_efficacy@", vals) %>% 
    paste(collapse = "\n")
  
  # Fix line breaks so that replacing snippets on the cluster works well
  snippet = gsub("\\n", "\\\\n", snippet)
  
  # Return complete seasonality snippet
  return(snippet)
  
}

# ------------------------------------------------------------------------------
# Generate XML snippet for vaccine efficacy decay
# ------------------------------------------------------------------------------
extract_vaccine_decay = function(o, p) {
  
  # Load XML snippet for vaccine efficacy decay section
  xml_template = readr::read_file(paste0(o$pth$input, 
                                         "/snippets/vaccine_decay.txt"))
  
  # Named vector to replace multiple placeholders
  key_vals = c("@half_life@" = as.character(p$vaccine$half_life),
               "@k_decay@"   = as.character(p$vaccine$k_decay))
  
  # Replace values in the XML snippet using sed patterns
  snippet = str_replace_all(xml_template, key_vals) %>% 
    paste(collapse = "\n")
  
  # Fix line breaks so that replacing snippets on the cluster works well
  snippet = gsub("\\n", "\\\\n", snippet)
  
  # Return complete vaccine efficacy decay snippet
  return(snippet)
  
}

# ------------------------------------------------------------------------------
# Generate XML snippet for cohort pattern
# ------------------------------------------------------------------------------
extract_cohort = function(o, p) {
  
  snippet_list = list(monitoring = "",
                      definition = "",
                      deployment = "",
                      removal = "")
  
  if(p$cohorts$include) {
    
    # ---- Monitoring of cohort ----
    
    # Import XML snippet
    xml_template = readr::read_file(paste0(o$pth$input, "/snippets/cohort_monitoring.txt"))
    
    # Replace id of cohort for tracking
    xml_template = str_replace(xml_template, "@cohort_tracking@", p$cohorts$tracking)
    
    # Fix line breaks so that replacing snippets on the cluster works well
    snippet_list$monitoring = gsub("\\n", "\\\\n", xml_template)
    
    if(!(p$cohorts$tracking == "vaccine")) {
      
      xml_template = readr::read_file(paste0(o$pth$input, "/snippets/cohort_definition.txt"))
      
      xml_template = str_replace(xml_template, "@cohort_definition@", xml_template)
      
      snippet_list$definition = gsub("\\n", "\\\\n", xml_template)
      
      xml_template = readr::read_file(paste0(o$pth$input, "/snippets/cohort_deployment.txt"))
      
      # Replace id of cohort for tracking
      xml_template = str_replace(xml_template, "@recruitment_start@", as.character(p$cohorts$begin_recruitment))
      
      snippet_list$deployment = gsub("\\n", "\\\\n", xml_template)
      
    }
    
    # ---- Removal from cohort ----
    
    if(p$cohorts$max_time > 0) {
      
      # Import XML snippet
      xml_template = readr::read_file(paste0(o$pth$input, "/snippets/cohort_removal.txt"))
      
      # Replace placeholder after which individual is removed from cohort  
      snippet_list$removal = str_replace(xml_template, "@cohort_max_time@", as.character(p$cohorts$max_time))
      
    }
    
  }

  return(snippet_list)
  
}

# ------------------------------------------------------------------------------
# Generate XML snippet for age group monitoring pattern
# ------------------------------------------------------------------------------
extract_age_groups = function(o, p) {
  
  age_groups = p$age_groups
  
  snippet = ""
  
  if(length(age_groups) > 0) {
    
    # ---- Monitoring of age groups (other than 0-99) ----
    
    # Import XML snippet
    xml_template = readr::read_file(paste0(o$pth$input, 
                                           "/snippets/agegroup_group.txt"))
    
    snippet = str_replace_all(xml_template, "@age_upperbound@", 
                              as.character(age_groups)) %>% paste(collapse = "\n")
    
    # Fix line breaks so that replacing snippets on the cluster works well
    snippet = gsub("\\n", "\\\\n", snippet)
    
  }
  
  return(snippet)
  
}

# ------------------------------------------------------------------------------
# Generate XML snippet for deployment pattern
# ------------------------------------------------------------------------------
extract_deployment = function(o, p) {
  
  if(grepl("perennial", o$analysis_name) & grepl("massvax", o$analysis_name)) {
    deployment_schedule = "seasonal"
  } else {
    deployment_schedule = p$vaccine_deployment$id
  }
  
  # Import roll out-specific deployment XML snippet
  xml_template = readr::read_file(paste0(o$pth$input, "/snippets/deployment_", 
                                         deployment_schedule, ".txt"))
  
  # Introduce restriction of vaccine deployment
  xml_snippet = readr::read_file(paste0(o$pth$input, "/snippets/deployment_restriction.txt"))
  
  if(p$cohorts$include == TRUE & p$cohorts$tracking == "intervention_cohort") {
    
    xml_snippet = str_replace(xml_snippet, "@cohort_tracking@", p$cohorts$tracking)
    
  } else {
    
    xml_snippet = ""
    
  }
  
  xml_template = str_replace(xml_template, "@deployment_restriction@", xml_snippet)
  
  
  # First, let's build the snippet for seasonal vaccination roll out
  if(p$vaccine_deployment$id == "seasonal") {

    deployment_template = xml_template
    
    # ---- Sort out timing of deployment ----
    
    # Determine peak month based on seasonality profile
    peak_season = which(p$seasonality$values == max(p$seasonality$values))
    
    # How many months in advance of the peak should vaccine be administered
    lead_time = 3.5
    
    # Determine date at which peak occurs
    peak_date = ymd(p$deployment_date) + months(peak_season)
    
    # Subtract lead time from peak date (convoluted approach with modulo operator to extract decimal number)
    deployment_init = peak_date - months(ceil(lead_time) + 1) +
      days((lead_time %% 1)*30 - 1)
    
    # Build yearly deployment schedule based on simulation horizon
    deployment_schedule = deployment_init + years(0:(p$lead_time + p$sim_duration - p$burn_in))
    
    # Ensure that deployment ends with end of the simulation horizon
    deployment_schedule = deployment_schedule[deployment_schedule <= p$end_date] %>% as.character()
    
    # ---- Deployment third dose ----
    
    # Import snippet for third dose
    xml_dose3 = readr::read_file(paste0(o$pth$input, "/snippets/deployment_seasonal_dose3.txt"))
    
    # Parse numbers as strings (required by str_replace_all)
    dose3_vals = p$vaccine_deployment$doses$initial$age %>% as.character()
    
    # Build sed pattern
    names(dose3_vals) = paste0("@", c("min_age", "max_age"), "@")
    
    # Replace values in the XML snippet using sed patterns
    dose3_snippet = str_replace_all(xml_dose3, dose3_vals)
    
    dose3_snippet = str_replace_all(dose3_snippet, "@deployment_dose3_individual@", 
                                    deployment_schedule) %>%
      paste(collapse = "\n")
    
    dose3_deployment = str_replace(deployment_template, "@deployment_schedule@", dose3_snippet)
    
    dose3_deployment = str_replace(dose3_deployment, "@deployment_name@", "initial")
    
    n_booster = p$vaccine_deployment$doses$booster$age %>% length() %>% as.character()
    
    dose3_deployment = str_replace_all(dose3_deployment, "@n_booster@", n_booster)
    
    if(p$cohorts$include == TRUE & p$cohorts$tracking == "intervention_cohort") {
      
      dose3_deployment = str_replace_all(dose3_deployment, "@coverage@", "1.0")
    
    }
    
    # ---- Deployment booster dose ----
    
    xml_booster = readr::read_file(paste0(o$pth$input, 
                                          "/snippets/deployment_seasonal_booster.txt"))
    
    # Extract schedule for booster doses
    booster_schedule = p$vaccine_deployment$doses$booster$age
    
    # Determine number of booster doses (set to 1 if only one age-pair provided)
    n_boosters = ifelse(is.list(booster_schedule), length(booster_schedule), 1)
    
    # Convert booster schedule to matrix format
    booster_matrix = matrix(unlist(booster_schedule), 
                            nrow = n_boosters, byrow = TRUE)

    # Create auxiliary vector for setting correct order of doses
    if(p$vaccine_deployment$booster_conditional) {
      
      # Only provide consecutive boosters if previous shots had been received    
      previous_doses = 1:n_boosters
      max_doses = previous_doses + 1 
      
    } else {
      
      # Only initial vaccination required
      previous_doses = rep(1, n_boosters)
      max_doses = rep(n_boosters + 1, n_boosters)
      
    }
    
    # Combine all of these properties into matrix form
    booster_df = cbind(booster_matrix, previous_doses, max_doses) %>%
      apply(MARGIN = 2, FUN = as.character) %>%
      matrix(ncol = 4)
    
    # Initialize list for booster xml snippets
    booster_xmls = c()
    
    for(i in 1:n_boosters) {
      
      # Shortcut to booster series i
      booster_vals = booster_df[i, ]
      
      # Build sed pattern
      names(booster_vals) = paste0("@", c("min_age", "max_age", "previous_doses", 
                                          "max_doses"), "@")
      
      # Replace values in the XML snippet using sed patterns
      this_booster_snippet = str_replace_all(xml_booster, booster_vals)
      
      # Expand deployment to deployment schedule
      this_booster_snippet = str_replace_all(this_booster_snippet, 
                                             "@deployment_booster_individual@", 
                                             deployment_schedule) %>%
        paste(collapse = "\n")
      
      # Assign name based on number of booster for readability
      this_booster_xml = str_replace(deployment_template, "@deployment_name@", 
                                     paste0("booster", i))  
      
      # Merge deployment schedule into deployment snippet
      this_booster_xml = str_replace(this_booster_xml, "@deployment_schedule@", 
                                     this_booster_snippet)  
      
      # Append to list of booster doses
      booster_xmls = c(booster_xmls, this_booster_xml)
      
    }
    
    # Concatenate booster snippets
    booster_deployment = booster_xmls %>% paste(collapse = "\n")
    
    # Concatenate third dose and booster doses
    xml_template = c(dose3_deployment, booster_deployment) %>%
      paste(collapse = "\n")
  
  } else if(p$vaccine_deployment$id == "age") { # Age-based deployment
    
    # ---- Deployment third dose ----
    
    # Import snippet for third dose
    
    if(grepl("massvax", o$analysis_name)) {
      
      
      deployment_template = xml_template
      
      duration_vaccination = year(p$end_date) - year(p$deployment_date)
      
      deployment_dates = (ymd(p$deployment_date) + years(0:duration_vaccination)) %>% 
        as.character()
      
      xml_dose3 = readr::read_file(paste0(o$pth$input, "/snippets/deployment_", 
                                          deployment_schedule, "_dose3.txt"))
      
      min_age  = seq(from = p$vaccine_deployment$doses$initial$age[1],
                     to   = p$vaccine_deployment$doses$initial$age[2],
                     by   = 1)
      
      age_rng  = cbind(min_age = min_age, max_age = min_age + 1) %>%
        as.data.table()
      
      n_booster = duration_vaccination
    
      dose3_vals = c(p$vaccine_deployment$doses$initial$age[1], 
                     p$vaccine_deployment$doses$initial$age[2], 1) %>%
        as.character()
      names(dose3_vals) = c("@min_age@", "@max_age@", "@n_booster@")
        
      dose3_snippet = str_replace_all(xml_dose3, dose3_vals)
      
      dose3_snippet = str_replace_all(dose3_snippet, "@deployment_dose3_individual@", 
                                      deployment_dates) %>%
        paste(collapse = "\n")
      
      dose3_deployment = deployment_template %>%
        str_replace("@deployment_name@", "initial") %>%
        str_replace("@deployment_schedule@", dose3_snippet)

      # ---> BOOSTER =--
      
      single_xmls = c()
      
      xml_booster = readr::read_file(paste0(o$pth$input, "/snippets/deployment_", 
                                            deployment_schedule, "_booster.txt"))
      
      for(i in 1:nrow(age_rng)) {
        
        booster_vals = c(age_rng[i, min_age], age_rng[i, max_age], "1", nrow(age_rng))
        names(booster_vals) = c("@min_age@", "@max_age@", "@previous_doses@", "@max_doses@")
      
        this_booster = str_replace_all(xml_booster, booster_vals) %>%
          str_replace_all("@deployment_booster_individual@", deployment_dates) %>%
          paste(collapse = "\n")
        
        booster_vals_2 = c(paste0("booster_", i), this_booster)
        names(booster_vals_2) = c("@deployment_name@", "@deployment_schedule@")
        
        this_deployment = str_replace_all(deployment_template, booster_vals_2)
        
        single_xmls = c(single_xmls, this_deployment)
          
      }
      
      booster_deployment = single_xmls %>% paste(collapse = "\n")
      
      xml_template = c(dose3_deployment, booster_deployment) %>%
        paste(collapse = "\n")
      
      
    } else {
      
      xml_dose3 = readr::read_file(paste0(o$pth$input, 
                                          "/snippets/deployment_age_dose3.txt"))
      
      # Parse numbers as strings (required by str_replace_all)
      dose3_vals = p$vaccine_deployment$doses$initial$age %>% as.character()
      
      dose3_vals = c(dose3_vals, "0", "1")
      
      # Build sed pattern
      names(dose3_vals) = c("@age@", "@min_dose@", "@max_dose@")
      
      # Replace values in the XML snippet using sed patterns
      dose3_snippet = str_replace_all(xml_dose3, dose3_vals)
      
      if(p$cohorts$include == TRUE & p$cohorts$tracking == "intervention_cohort") {
        
        dose3_snippet = str_replace(dose3_snippet, "@coverage@", as.character("1.0"))
        
      }
      
      xml_template = str_replace(xml_template, "@deployment_dose3@", dose3_snippet)
      
      # ---- Deployment booster dose ----
      
      # Import snippet for booster dose
      xml_booster = readr::read_file(paste0(o$pth$input, 
                                            "/snippets/deployment_age_booster.txt"))
      
      # Extract vaccination schedule
      booster_schedule = p$vaccine_deployment$doses$booster$age
      
      # Determine number of booster shots
      n_boosters = length(booster_schedule)
      
      # Create auxiliary vector for setting correct order of doses
      min_doses = 1:n_boosters
      max_doses = min_doses + 1 
      
      # Combine all of these properties into matrix form
      booster_df = cbind(booster_schedule, min_doses, max_doses) %>%
        apply(MARGIN = 2, FUN = as.character) %>%
        matrix(ncol = 3)
      
      # Initialize vector replacement patterns
      booster_template = c()
      
      # Loop through each booster shot
      for(i in 1:n_boosters) {
        
        # Shortcut to booster i
        this_booster = booster_df[i, ]
        
        # Assign name to each column to ensure correct input format for str_replace_all
        names(this_booster) = c("@age@", "@min_dose@", "@max_dose@")
        
        # Replace placeholders in template
        booster_template = c(booster_template, str_replace_all(xml_booster, 
                                                               this_booster))
        
      }
      
      # Collapse vector of strings to single string
      xml_booster = booster_template %>% paste(collapse = "\n")
      
      # Replace values in the XML snippet using sed patterns
      xml_template = str_replace(xml_template, "@deployment_booster@", xml_booster)
      
    }
    
  }
  
  # Fix line breaks so that replacing snippets on the cluster works well
  snippet = gsub("\\n", "\\\\n", xml_template)
  
  # Return complete deployment snippet
  return(snippet)
  
}
