################################################################################
# SENSITIVITY ANALYSIS
#                                                                             
# Run sensitivity analysis based on Gaussian Process emulator
#
################################################################################

#-------------------------------------------------------------------------------
# Main function for sensitivity analysis
#-------------------------------------------------------------------------------
sensitivity_analysis = function(o) {

  # Set seed as specified in list of options (times 5)
  set.seed(as.integer(o$seed*5))
  
  # Load parameters
  p = load_parameters(o)
  
  # If scenario does not contain any LHC samples, skip this process step
  if(!hasLHC(p))
    return(NULL)
  
  message("* Running global sensitivity analysis")
  
  # Shortcut to names of parameters that are sampled from LHC
  vars = p$parameters$continuous$parameter
  
  # Parameter bounds in matrix form
  bounds = p$parameters$continuous[, .(lb, ub)] %>% 
    as.matrix()
  
  # List of scenarios that are simulated
  fids = p$parameters$fixed_grid$fid
  
  # List of directories that are in gp_models folder (shows number of samples from ASD)
  dirs = list.dirs(o$pth$gp_models, full.names = FALSE)
  idx = gsub(paste0(dirs[1], "sample_"), "", dirs) %>% as.numeric() %>% max(na.rm = TRUE)
    
  # Pre-allocate some space for results
  indices_list = plot_df_list = list()

  # Run sensitivity analysis for every LHC outcome, ...
  for(this_outcome in p$lhc_outcomes) {
    
    message(paste0("  ~ Outcome: ", this_outcome))
    
    pb = startProgressBar(length(fids))
    
    for(f in seq_along(fids)) {
      
      # Load latest Gaussian Process results file
      gp_model = readRDS(paste0(o$pth$gp_models, "/", "sample_", idx, "/", fids[f], ".", this_outcome, ".rds"))
      
      # ---- Variance-decomposition analysis (Sobol indices) ----
      
      # Calculate the Sobol indices based on model emulator
      sobol_indices = calc_sobol_idx(o, p, gp_model)
      
      # Format first-order and total Sobol indices for storing in data.table
      setDT(sobol_indices$S, keep.rownames = "parameter")
      sobol_indices$S = sobol_indices$S[, .(fid = fids[f],
                                            parameter, 
                                            idx_order = "first", 
                                            value = original,
                                            se = `std. error`,
                                            lower = `min. c.i.`,
                                            upper = `max. c.i.`)]
      
      setDT(sobol_indices$T, keep.rownames = "parameter")
      sobol_indices$T = sobol_indices$T[, .(fid = fids[f],
                                            parameter, 
                                            idx_order = "total", 
                                            value = original,
                                            se = `std. error`,
                                            lower = `min. c.i.`,
                                            upper = `max. c.i.`)]
      
      indices_list[[paste0(f, ".", this_outcome)]] = rbindlist(list(sobol_indices$S, sobol_indices$T))
      
      # Ensure Sobol indices are non-negative and sum to 1
      indices_list[[paste0(f, ".", this_outcome)]][, value := pmax(0, value)]
      
      indices_list[[paste0(f, ".", this_outcome)]][, norm_value := value/sum(value), by = list(fid, idx_order)]
      
      indices_list[[paste0(f, ".", this_outcome)]][, outcome := this_outcome]
      
      # ---- Global sensitivity analysis ----
      
      sobol_df = cbind(sobol_indices$X, output = sobol_indices$y)
      
      # Preallocate list for all quantiles
      quantiles_list = list()
      
      # Loop through all parameters
      for(var in vars) {
        
        # Cut the parameter space into bins
        var_bins = seq(from = min(sobol_df[[var]]), 
                       to   = max(sobol_df[[var]]), 
                       length = o$n_bins + 1)
        
        # Preallocate matrix for quantiles for this variable
        quantiles = matrix(NA, nrow = o$n_bins, ncol = 6)
        
        # Loop through the bins
        for (j in 1:o$n_bins) {
          
          # Outputs associated with this bin
          split_df = sobol_df %>%
            select(input = !!var, output) %>%
            filter(input > var_bins[j],
                   input < var_bins[j + 1])
          
          # Estimate the median output in each fraction of the parameter space
          quantiles[j, ] = c(quantile(split_df$output, c(o$quantiles[1], 0.5, o$quantiles[2])), 
                             min(split_df$output), max(split_df$output), mean(split_df$output))
        }
        
        # Transform back to original parameter scale
        x_transform = seq(from = bounds[match(var, vars), 1],
                          to   = bounds[match(var, vars), 2],
                          length.out = o$n_bins)
        
        # Story aggregated summary statistics for each bin
        quantiles_list[[var]] = quantiles %>%
          as_named_dt(c("lower", "median", "upper", "min", "max", "mean")) %>%
          mutate(variable = var,
                 x        = x_transform)
        
      }
      
      # Convert to data.table and assign scenario and LHC outcome
      plot_df_list[[paste0(fids[f], ".", this_outcome)]] = rbindlist(quantiles_list) 
      
      plot_df_list[[paste0(fids[f], ".", this_outcome)]][, fid := fids[f]]
      
      plot_df_list[[paste0(fids[f], ".", this_outcome)]][, outcome := this_outcome]
    
      setTxtProgressBar(pb, f)
     
    }
    
    close(pb)
    
  }
  
  # Concatenate all outputs
  indices_list = rbindlist(indices_list)
  
  plot_df_list = rbindlist(plot_df_list)
  
  output = list(sobol_indices    = indices_list,
                sensitivity_data = plot_df_list,
                opts             = o,
                parameters       = p)

  # Save output of sensitivity analysis
  saveRDS(output, paste0(o$pth$manuscript, "global_sensitivity_analysis.rds"))

  return(NULL)
}

#-------------------------------------------------------------------------------
# Calculate Sobol indices from trained model
#-------------------------------------------------------------------------------
calc_sobol_idx = function(o, p, trained_model) {

  # Shortcut to model emulator
  gp_model = trained_model$gp_model
  
  # Wrapper function for model prediction
  sobol_fn = function(X) {
    return(predict(x = as.matrix(X), object = gp_model)$mean)
  }
    
  # Parameter bounds
  bounds = p$parameters$continuous[, .(lb, ub)] %>% as.matrix() %>% unname()
  
  # Construct two LHC samples and normalize to [0,1]
  sample1 = lhs(o$sobol_n_samples, bounds) %>% as.data.frame()
  colnames(sample1) = p$parameters$continuous$parameter
  
  sample1 = sample1 %>%
    normalise_params(params = p$parameters$continuous$parameter,
                     bounds = p$parameters$continuous[, .(lb, ub)]) 

  sample2 = lhs(o$sobol_n_samples, bounds) %>% as.data.frame()
  colnames(sample2) = p$parameters$continuous$parameter

  sample2 = sample2 %>%
    normalise_params(params = p$parameters$continuous$parameter,
                     bounds = p$parameters$continuous[, .(lb, ub)])
  
  # Estimate Sobol indices
  sobol_idx = soboljansen(model = sobol_fn, 
                          as.data.frame(sample1), 
                          as.data.frame(sample2), 
                          nboot = o$sobol_n_bootstraps)
  
  # Return Sobol indices
  return(sobol_idx)
  
}