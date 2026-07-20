################################################################################
# EMULATOR
#
# All functions related to fitting Gaussian processes and assessing emulator 
# performance
#
################################################################################

# ------------------------------------------------------------------------------
# Main function to Fit Gaussian process to simulation results
# ------------------------------------------------------------------------------
run_gp = function(o, p, m) {

  # Skip function if no emulator is needed or step not specified in options list
  if(!hasLHC(p))
    return(NULL)
  
  if(!(7 %in% o$do_step)) return(NULL)
  
  message("* Fitting GP emulators")
  
  # Load list of simulations
  scenarios_df = readRDS(file.path(h(m)$pth$parameters, "all_scenarios.rds"))
  
  # Shortcut to parameters that were sampled from LHC
  continuous_parameters = p$parameters$continuous$parameter
  
  # Extract list of settings for which emulator are fitted one after another
  fids = unique(scenarios_df$fid)
  
  # Check which emulator results are already available
  fids = check_existing(o, fids, folder = "gp_models", m = m)
  
  error_msg = c()
  
  # Loop through LHC outcomes (mainly used: selection_coefficient)
  for(this_outcome in p$lhc_outcomes) {
    
    message(paste0("  ~ Outcome: ", this_outcome))
    
    pb = startProgressBar(length(fids))
    
    # Loop over scenarios, ...
    for(i in seq_along(fids)) {
      
      # Extract paths to parameters for each sample and concatenate
      file_paths_params = paste0(o$pth$parameters, 
                                              "sample_", 0:m$sample_num, "/all_scenarios.rds")
      all_params = lapply(file_paths_params, FUN = readRDS) %>%
        rbindlist(use.names = TRUE)
      
      # Extract paths to summarized results and concatenate
      file_paths_results = paste0(o$pth$summarised, "sample_", 0:m$sample_num, "/", fids[i], ".rds")
      
      output_df = lapply(file_paths_results, FUN = readRDS) %>%
        rbindlist(use.names = TRUE)
      
      # Only include LHC outcome of this iteration
      output_df = output_df[outcome == this_outcome & !is.na(value)]
      
      # Skip LHC outcome if not available for this scenario
      if(nrow(output_df) == 0) {
        error_msg = c(error_msg, 
                      paste0("      ! Fittung emulator for scenario ", fids[i], " not possible due to missing data!"))
        
        setTxtProgressBar(pb, i)
        next
      }
      
      # Merge list of parameters to LHC results
      output_df = output_df %>%
        left_join(all_params, by = c("fid", "scenario")) %>%
        as.data.table()
      
      # Append combination of sample ID and continuous sample ID
      output_df[, sid_cid := paste0(sid, "_", cid)]
      
      # Select subset of columns
      cols = c("scenario", "sid_cid", continuous_parameters, "value")
      output_df = output_df[, ..cols]

      # Normalize parameters to range of [0, 1]
      norm_cols = c("sid_cid", continuous_parameters, "value")
      
      output_std_df = output_df[, ..norm_cols] %>% 
        normalise_params(continuous_parameters, p$parameters$continuous[, .(lb, ub)])
      
      # Throw error if something was transformed incorrectly
      if(!all(output_std_df$value == output_df$value))
        stop("! Value was transformed by accident.")
      
      # Make sure the order of parameters and results is not perturbed
      keys = unique(output_std_df$sid_cid) %>% unique()
      
      # Store order of continuous settings
      # order_keys = output_df$sid_cid
      
      # If it is the first step in adaptive sampling routine, ...
      if(m$sample_num == 0) {
        
        n_test = round(o$gp_test_split * length(keys))
        
        test_idx = sample.int(length(keys), n_test, replace = FALSE)
        test_keys = keys[test_idx]
        
        train_keys = keys[!(keys %in% test_keys)]
        
      } else {
        
        init_emulator = readRDS(paste0(m$sample_0$pth$gp_models, "/", fids[i], ".", this_outcome, ".rds"))
        
        test_keys = init_emulator$test_keys
        
        train_keys = keys[!(keys %in% test_keys)]
        
      }
      
      # ---- Train data ----
      
      # Standardize data for training and prepare for emulator fitting
      X_train = output_std_df[sid_cid %in% train_keys, ..continuous_parameters] %>% 
        as.matrix() %>% unname()
      
      Z_train = output_std_df[sid_cid %in% train_keys, value]
      
      train_data = find_reps(X_train, Z_train, rescale = FALSE, normalize = FALSE,
                             return.Zlist = FALSE)
      
      # ---- Test data ----
      
      # Prepare testing data for comparison
      test_data = output_std_df[sid_cid %in% test_keys] %>% 
        group_by(sid_cid, across(all_of(continuous_parameters))) %>% 
        summarize(value = mean(value), .groups = "drop") %>% 
        as.data.table()
      
      X_test = test_data[, ..continuous_parameters] %>% as.matrix() %>% unname()
      Z_test = test_data[, value]
      
      # Previously, R and Python were implemented for fitting GP emulator
      if(o$gp_engine == "R") {
        
        gp_model = mleHetGP(X = list(X0 = train_data$X0, Z0 = train_data$Z0, mult = train_data$mult),
                            Z = train_data$Z, covtype = o$gp_kernel, maxit = o$gp_max_iter)
        
        Z_train_predicted = predict(gp_model, X_train)
        
        Z_test_predicted = predict(object = gp_model, X_test)
        
        rmse = scores(gp_model, X_test, Z_test, return.rmse = TRUE)
        
      } else {
        
        optim_params = run_optuna_study(train_data, list(X0 = X_test, Z0 = Z_test))
        
        gp_model = fit_gp_python(train_data, optim_params)
        
        joblib = import("joblib")
        joblib$dump(gp_model, file.path(h(m)$pth$gp_models, paste0(fids[f], ".", this_outcome, ".pkl")))
        
        Z_train_predicted = list(mean = predict_gp_python(gp_model, X_train))
        
        Z_test_predicted = list(mean = predict_gp_python(gp_model, X_test))
        
      }
      
      # Collate training and testing results
      results_train = data.table(set = "train", 
                                 actual = Z_train, 
                                 predicted = Z_train_predicted$mean)
      
      results_test = data.table(set = "test", 
                                actual = Z_test, 
                                predicted = Z_test_predicted$mean)
      
      results = rbindlist(list(results_train, results_test), use.names = FALSE)
      
      # Ensure diagnostic plot is a perfect square
      lims_x = min(results$actual, results$predicted)
      lims_y = max(results$actual, results$predicted)
      
      # Collate all emulator outputs
      emulator = list(gp_model    = gp_model,
                     #  std_params  = X_std$std_params,
                      X_train     = X_train,
                      Z_train     = Z_train,
                      Z_train_predicted = Z_train_predicted$mean,
                      X_test      = X_test,
                      Z_test      = Z_test,
                      Z_test_predicted = Z_test_predicted$mean,
                      train_keys  = train_keys,
                      test_keys   = test_keys,
                      params      = all_params,
                      rmse        = rmse)
      
      # ---- Further diagnostic of extreme values ----
      
      diag_train = cbind(output_std_df[sid_cid %in% train_keys, .(sid_cid)], value = (predict(gp_model, X_train)$mean - Z_train)^2)
      
      diag_top5 = diag_train %>%
        group_by(sid_cid) %>%
        arrange(-value) %>%
        slice_head(n = 1) %>%
        ungroup() %>%
        arrange(-value) %>%
        slice_head(n = 10)
      
      input_train_reduced_df = output_df[sid_cid %in% train_keys, .(sid_cid, access, ar_mkr_resistant, eir, coverage, lf_resistance_level)] %>% unique()
      
      input_train_reduced_df[sid_cid %in% diag_top5$sid_cid]
      
      saveRDS(emulator, file.path(h(m)$pth$gp_models, paste0(fids[i], ".", this_outcome, ".rds")))
      
      results$set = factor(results$set, levels = c("train", "test"),
                           labels = c("Train", "Test"))
      
      # Calculate correlation in test dataset
      performance_text = paste0("Correlation in test set: ", round(cor(Z_test, Z_test_predicted$mean), digits = 3))
      
      # Put together diagnostic plot
      diagnostic_plot = ggplot(results, aes(x = actual, y = predicted, color = set)) +
        geom_point(size = 3, alpha = .8) + 
        scale_color_manual(name = "Dataset", values = o$pal[2:length(o$pal)]) +
        scale_x_continuous(limits = c(lims_x, lims_y), name = "Actual") +
        scale_y_continuous(limits = c(lims_x, lims_y), name = "Predicted") + 
        geom_abline(slope = 1) +
        annotate("text", label = performance_text, x = 0.9*lims_y, y = 0.9*lims_x, size = 4, hjust = 1) +
        theme_amr() +
        theme(legend.position = "bottom",
              axis.title.y = element_text(margin = margin(0, 15, 0, 5)),
              axis.title.x = element_text(margin = margin(15, 0, 0, 0)),
              plot.margin = margin(10, 15, 5, 5))
      
      # Save diagnostic plot
      ggsave(file.path(h(m)$pth$gp_models, paste0(fids[i], ".", this_outcome, "_", o$gp_engine, "_diagnostic.png")), 
             plot = diagnostic_plot, width = 8, height = 8)
      
      setTxtProgressBar(pb, i)
      
    }
    
    close(pb)
    
    if(length(error_msg) > 0)
      message(paste0(error_msg, sep = "\n"))
    
  }
  
  return(NULL)
}

# -----------------------------------------------------------
# Train a GP regression with hetGP given a training data set.
# -----------------------------------------------------------
gp_train = function(train_data) {

  # Define indices of parameters and response
  n_params = ncol(train_data) - 1
  param_idx = 1:n_params
  response_idx = n_params + 1

  # Find and remove duplicates
  prep_data = find_reps(
    X = as.matrix(train_data[, param_idx]),
    Z = as.matrix(train_data[, response_idx]),
    rescale = FALSE, normalize = FALSE)

  # Prepare input for mleHetGP model
  X = list(
    X0 = as.matrix(prep_data$X0),
    Z0 = as.matrix(prep_data$Z0),
    mult = prep_data$mult)

  # Run GP model
  gp_model = mleHetGP(
    X = X, Z = prep_data$Z,
    covtype = "Matern5_2")
  
  # Return gp model
  return(gp_model)
}

# -----------------------------------------------
# Run GP with test data to determine model error.
# -----------------------------------------------
gp_test = function(gp_model, train_data, test_data_2) {

  # Define indices of parameters and response
  n_params = ncol(train_data) - 1
  param_idx = 1:n_params
  response_idx = n_params + 1

  # Apply the model on the test data
  predict_test = predict(x = as.matrix(test_data_2[, param_idx]), object = gp_model)
  predict_test = predict_test$mean #### look at the mean of predictions

  # The actual response associated with the test data
  actual_test = test_data_2[, response_idx]

  # Repeat this for the training data
  predict_train = predict(x = as.matrix(train_data[, param_idx]), object = gp_model)
  predict_train = predict_train$mean

  # The actual response associated with the test data
  actual_train = train_data[, response_idx]

  # Concatenate output into list
  gp_output = list(
    actual_test = actual_test,
    predict_test = predict_test,
    actual_train = actual_train,
    predict_train = predict_train)
  
  # Return the predicted data and true value
  return(gp_output)
}

#--------------------------------------------------------------------------
# function to Reorganise the gp accuracy result from a list to a dataframe.
#--------------------------------------------------------------------------
assess_emulator_performance = function(o) {
  
  p = load_parameters(o)
  
  fids = p$parameters$fixed_grid$fid
  
  dirs = list.dirs(o$pth$gp_models)
  
  idx = gsub(paste0(dirs[1], "/sample_"), "", dirs) %>% as.numeric()
  
  n_samples = idx[!is.na(idx)]
  
  precision_list = list()
  
  for(f in seq_along(fids)) {
  
    precision_list[[f]] = list()
    
    for(i in 1:length(n_samples)) {
      
      results_gp = readRDS(paste0(o$pth$gp_models, paste0("sample_", n_samples[[i]], "/", fids[f], ".rds")))
      
      test_data = cbind(results_gp$Z_test, results_gp$Z_test_predicted)
      
      precision_list[[f]][[i]] = data.table(fid = fids[f],
                                       iteration = n_samples[i],
                                       precision = corr(test_data))
      
    }
      
    precision_list[[f]] = rbindlist(precision_list[[f]])
    
  }

  precision_df = rbindlist(precision_list)
  
  precision_plot = ggplot(precision_df, aes(x = factor(iteration), y=precision, fill = fid, label = round(precision, 4))) +
    geom_bar(stat = "identity", fill = o$pal[2]) +
    facet_wrap(~ fid, ncol = 4) +
    scale_x_discrete(name = "Adaptive sampling iteration") +
    scale_y_continuous(name = "Correlation in test set", expand = expansion(mult = c(0,0))) +
    geom_text(size = 3.5, position = position_stack(vjust = 0.5), color = "white") +
    theme_amr() +
    theme(axis.title.y = element_text(margin = margin(0, 10, 0, 0)))

  ggsave(paste0(o$pth$manuscript, "precision_emulators.png"), precision_plot, width = 8, height = 4)
    
}

fit_gp_python = function(train_data, optim_params) {
  
  # Import required modules
  sklearn = import("sklearn")
  
  # Access GaussianProcessRegressor and kernel options
  gaussian_process = sklearn$gaussian_process$GaussianProcessRegressor
  kernels = sklearn$gaussian_process$kernels
  
  if(optim_params$kernel == "Matern") {
    
    kernel = kernels$Matern(length_scale = optim_params$length_scale, nu = optim_params$nu)  
    
  } else if(optim_params$kernel == "RBF") {
    
    kernel = kernels$RBF(length_scale = optim_params$length_scale)
    
  }
  
  # Create a GaussianProcessRegressor instance
  gp_model = gaussian_process(kernel = kernel, alpha = optim_params$alpha)
  
  # Fit the Gaussian Process model
  gp_object = gp_model$fit(train_data$X0, train_data$Z0)
  
  return(gp_object)
  
}

predict_gp_python = function(gp_model, input_data) {
  
  input_predicted = gp_model$predict(input_data, return_std = TRUE)[[1]]
  
  return(input_predicted)
  
}

run_optuna_study = function(train_data, test_data) {
  
  optuna = import("optuna")
  
  sklearn = import("sklearn")
  
  # Access GaussianProcessRegressor and kernels
  gaussian_process = sklearn$gaussian_process$GaussianProcessRegressor
  
  kernels = sklearn$gaussian_process$kernels
  
  objective = function(trial) {
    # Choose a kernel type
    kernel_type = trial$suggest_categorical("kernel", c("RBF", "Matern"))
    
    # Define hyperparameters for each kernel type
    if (kernel_type == "RBF") {
      length_scale = trial$suggest_float("length_scale", 0.1, 10.0)
      kernel = kernels$RBF(length_scale = length_scale)
      
    } else if (kernel_type == "Matern") {
      length_scale = trial$suggest_float("length_scale", 0.1, 10.0)
      nu = trial$suggest_categorical("nu", c(0.5, 1.5, 2.5))  # Nu options for Matern kernel
      kernel = kernels$Matern(length_scale = length_scale, nu = nu)
      
    }
    
    # Suggest a noise level (alpha)
    alpha = trial$suggest_float("alpha", 1e-10, 1e-2, log = TRUE)
    
    # Create GaussianProcessRegressor with the selected kernel and alpha
    gp_model = gaussian_process(kernel = kernel, alpha = alpha)
    
    # Fit the model to the training data
    gp_model$fit(train_data$X0, train_data$Z0)
    
    # Evaluate the model on the validation data (R² score as an example objective)
    score = gp_model$score(test_data$X0, test_data$Z0)
    
    test_predicted = gp_model$predict(test_data$X0)
    
    rmse = sklearn$metrics$root_mean_squared_error(test_predicted, test_data$Z0)
    
    return(rmse)
  }
  
  # Create the Optuna study
  study = optuna$create_study(direction = "minimize")
  
  # Run the optimization
  study$optimize(objective, n_trials = 50)
  
  best_params = study$best_params
  
  return(best_params)
  
}
