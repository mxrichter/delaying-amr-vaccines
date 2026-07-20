################################################################################
# ADAPTIVE SAMPLING
#
# Improve accuracy of Gaussian Process (GP) model emulator
#
################################################################################

# ------------------------------------------------------------------------------
# Parent function for all adaptive sampling processes.
# ----------------------------------------------------
adaptive_sampling = function(o, m) {
  
  # Skip process step if not specified in list of options
  if(!(8 %in% o$do_step)) return(m)
  
  message("* Running adaptive sampling algorithm")

  # Load YAML parameters
  p = load_parameters(o)
  
  # Shortcut to parameter grid
  fixed_parameter_grid = p$parameters$fixed_grid
  
  # List of unique scenarios
  fids = unique(fixed_parameter_grid$fid)
  
  # Only run process for 'selection coefficient'
  lhc_outcome = p$lhc_outcomes[1]
  
  # Loop through unique scenarios
  for(i in seq_along(fids)) {
    
    # Shortcut to latest emulator
    last_emulator = readRDS(file.path(h(m)$pth$gp_models, paste0(fids[i], ".", lhc_outcome, ".rds")))
    
    # Sample parameters from Latin hypercube
    new_samples = tgp::lhs(10000, p$parameters$continuous[, .(lb, ub)] %>% as.matrix()) %>%
               as.data.table()
    
    # Assign parameter names as column IDs
    colnames(new_samples) = p$parameters$continuous$parameter
    
    # Normalize parameters to scale [0, 1]
    X = new_samples %>%
      normalise_params(p$parameters$continuous$parameter, 
                       p$parameters$continuous[, .(lb, ub)]) %>% as.matrix()
    
    # Use fitted GP to predict selection coefficient
    predictions = predict(object = last_emulator$gp_model, X)
    
    # Append predicted variance to parameter list
    samples_weighted = cbind(new_samples, weight = predictions$sd2)
    
    # Take parameter sets with highest predicted variance
    new_parameters = samples_weighted[order(-weight)][1:o$n_adaptive_samples]
    
    # Remove weight column
    new_parameters[, weight := NULL]
   
    # Store new list of parameters
    m[[paste0("sample_", m$sample_num)]]$new_parameters[[fids[i]]] = new_parameters
    
  }

  # Return to the user all GP model and prediction
  return(m)
}

# ---------------------------------------------------------------------
# Concatenate sample folder name
# ---------------------------------------------------------------------
ps = function(s, i) {
  return(s[[paste0("sample_",i)]])
}