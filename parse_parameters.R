################################################################################
# PARSE PARAMETERS
#
# Build parameter tables for OpenMalaria simulations
#
################################################################################

# ------------------------------------------------------------------------------
# Generate the input table for the simulation
# ------------------------------------------------------------------------------
sample_parameters = function(o, p, m) {
  
  # Skip step if not specified in list of options
  if(!(1 %in% o$do_step)) return(NULL)

  message("  - Sampling from parameter space")
  
  # If this is the first iteration (or if adaptive sampling is switched off), 
  # use LHS to generate parameter table
  if (m$sample_num == 0) {
    
    # ---- Fixed parameters ----
    
    # Merge fixed parameters into one data.table
    eir_pth = paste0(o$pth$output, "/", p$eir_filter, "/8_manuscript/manuscript_data_emulator.rds")
    
    # Use stored access and EIR values to build pre-defined PfPR2-10 settings
    if(p$eir_filter != "" & file.exists(eir_pth)) {
      
      # Shortcut to parameter table
      eir_scenarios = readRDS(eir_pth)$eir_aligned[, .(access, eir = input_eir)]
      
      # Put together parameter table
      param_fitted = p$parameters$fixed
      param_fitted$eir = unique(eir_scenarios$eir) %>% sort()
      param_fitted$access = unique(eir_scenarios$access) %>% sort()
      
      # Cross-join to create all possible combinations
      param_fixed_df = do.call(CJ, param_fitted)
      
      # Only use pre-defined access and EIR combinations
      param_fixed_df = param_fixed_df %>%
        semi_join(eir_scenarios, by = c("access", "eir")) %>%
        as.data.table()
      
      # Output diagnostic messages
      message(paste0("    ! Replacing EIR and access values specified in scenario yaml by values from '", p$eir_filter, "'"))
      message(paste0("    ! Access: [", paste0(param_fitted$access, collapse = ", "), "]"))
      message(paste0("    ! EIR: [", paste0(param_fitted$eir, collapse = ", "), "]"))
      
    } else {
      
      # Otherwise, just use cross-join to construct parameter table
      param_fixed_df = do.call(CJ, p$parameters$fixed)
      
    }
    
    # If resistance scenarios are pre-defined (combinations of ART and LF resistance)
    if(p$resistance_scenarios != "") {
      
      # File path to resistance scenarios
      file_name = paste0(o$pth$input, "/scenarios/", p$resistance_scenarios)
      
      # If this file exists, ...
      if(file.exists(file_name)) {
        
        # Load resistance scenarios as data.table
        resistance_scenarios = fread(file_name)[, .(ar_mkr_resistant, lf_resistance_level)]
        
        # Restrict parameter table to scenarios that are specified in resistance scenario table
        param_fixed_df = param_fixed_df %>%
          semi_join(resistance_scenarios, by = c("ar_mkr_resistant", "lf_resistance_level")) %>%
          as.data.table()
        
      }
      
    }
    
    # If scenario is treatment failure scenario, use five-digit scenario codes
    if(grepl("treat", o$analysis_name)) {
      
      param_fixed_df[, fid := paste0("f", str_pad(.I, 5, pad = "0"))]
      
    } else {
      
      param_fixed_df[, fid := paste0("f", str_pad(.I, 4, pad = "0"))]
    }
    
    # Append identifier for set of fixed parameters
    setkey(param_fixed_df, fid)
    
    p$parameters$fixed_grid = param_fixed_df
    
    # ---- Continuous parameters ----
    
    if(hasLHC(p)) {
      
      # Generate LHS between parameter bounds
      param_continuous_df = lhs(o$lhc_samples, p$parameters$continuous[, .(lb, ub)] %>% as.matrix()) %>%
        as.data.table()
      
      # Assign names of continuous parameters to parameter data.table
      colnames(param_continuous_df) = p$parameters$continuous$parameter
      
      param_continuous_df[, cid := paste0("c", str_pad(.I, 4, pad = "0"))]
      
      setkey(param_continuous_df, cid)
      
      # ---- Finalize parameter data.table ----
      
      # Cross-join of fixed and continuous parameter sets (convoluted approach, find alternative if possible)
      param_df = Reduce(merge, list(param_fixed_df %>% as.data.frame(), param_continuous_df %>% as.data.frame()))
      
    } else {
      
      # Default value for CID column if no parameters are sampled
      param_fixed_df[, cid := "c0000"]
      
      param_df = param_fixed_df %>% as.data.frame()
      
    }
    
  } else {
    
    # Load parameter space created in previous adaptive sampling iteration
    new_param_space = h(m, m$sample_num - 1)$new_parameters
    
    # Pre-allocate empty list
    new_param_df = list()
    
    # Build new parameter table based on results from adaptive sampling iteration
    for(i in names(new_param_space)) {
      
      param_continuous_df = cbind(p$parameters$fixed_grid[fid == i], new_param_space[[i]])
      
      param_continuous_df[, cid := paste0("c", str_pad(.I, 4, pad = "0"))]  
      
      new_param_df[[i]] = param_continuous_df
      
    }
    
    param_df = rbindlist(new_param_df) %>% as.data.frame()
    
  }
  
  # Create data.table for seeds
  seeds = data.table(seed = 1:o$n_seeds)
  
  # Cross-join parameter data.table with seeds
  param_df = merge(param_df, seeds) %>% as.data.table()
  
  # Round value for access to treatment
  param_df[, access := round(access, digits = 10)]
  
  # Assign separate SID for sample iteration
  param_df[, "sid" := paste0("s", str_pad(m$sample_num, 4, pad = "0"))]
  
  # Create unique scenario name
  param_df[, "scenario" := paste(sid, fid, cid, seed, sep = "_")]
  
  # Save final parameter table
  saveRDS(param_df, paste0(h(m)$pth$parameters, "/", "all_scenarios.rds"))
  
  return(NULL)
  
}

# ------------------------------------------------------------------------------
# Check whether any parameters are sampled from Latin hypercube
# ------------------------------------------------------------------------------
hasLHC = function(p) {
  
  return(nrow(p$parameters$continuous) > 0)
  
}