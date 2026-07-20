################################################################################
# AUXILIARY
#
# A series of helpful R functions
#
# Author: Andrew J. Shattock(andrewjames.shattock@swisstph.ch)
#
# Adapted by: Thiery Masserey (thiery.masserey@swisstph.ch), and
#             Maximilian Richter (maximilian.richter@swisstph.ch)
#
################################################################################

# ---------------------------------------------------------
# Assert that a condition is true. Throw error if not.
# ---------------------------------------------------------
assert = function(condition, errMsg) {
  if (!condition) {
    stop(errMsg)
  }
}

# ---------------------------------------------------------
# Set as datatable and rename columns in one line (by AJ Shattock)
# ---------------------------------------------------------
as_named_dt = function(x, new_names) {
  
  # Convert to datatable
  dt = as.data.table(x)
  
  # Check new names are correct length
  old_names = names(dt)
  if (length(old_names) != length(new_names))
    stop("Inconsistent number of column names provided")
  
  # Set new column names
  named_dt = setnames(dt, old_names, new_names)
  
  return(named_dt)
}

# ---------------------------------------------------------
# Clear the console.
# ---------------------------------------------------------
clc = function() {
  cat("\014")
}

# ---------------------------------------------------------
# Clear all figures.
# ---------------------------------------------------------
clf = function() {
  graphics.off()
}

# ---------------------------------------------------------
# Clear environment (and also console and figures by default).
# ---------------------------------------------------------
clear = function(clearConsole = TRUE, clearFigs = TRUE) {
  if (clearConsole == TRUE) {
    clc()
  }
  if (clearFigs == TRUE) {
    clf()
  }
  rm(list = ls(envir = globalenv()), envir = globalenv())
}

# ---------------------------------------------------------
# Is an element contained within a vector? Returns boolean.
# ---------------------------------------------------------
contains = function(vec, el) {
  bool = is.element(el, vec) # Simply for readability
  return(bool)
}

# ---------------------------------------------------------
# Create a log file (see function wait_for_jobs)
# ---------------------------------------------------------
create_bash_log = function(log_file) {
  if (file.exists(log_file)) file.remove(log_file)
  Sys.sleep(0.1)
  file.create(log_file)
  Sys.sleep(0.1)
}

# ---------------------------------------------------------
# Create a log file (see function wait_for_jobs); by AJ Shattock
# ---------------------------------------------------------
create_bash_log2 = function(pth, log = NULL, err = NULL) {
  for (this_name in c(log, err)) {
    this_file = paste0(pth, this_name)
    if (file.exists(this_file)) file.remove(this_file)
    Sys.sleep(0.1)
    file.create(this_file)
    Sys.sleep(0.1)
  }
  return(paste0(pth, log))
}

# ---------------------------------------------------------
# Create dataframe with custom column headings in one line.
# ---------------------------------------------------------
dataFrame = function(x, colNames = NULL, rowNames = NULL, ...) {
  df = data.frame(x, row.names = rowNames, ...) # Create dataframe
  if (!is.null(colNames)) {
    colnames(df) = colNames
  } # Rename column headers if need be
  return(df)
}

# ---------------------------------------------------------
# Drop one or more columns from a dataframe - improved readability.
# ---------------------------------------------------------
df_drop = function(df, drop_cols) {
  drop_idx = (names(df) %in% drop_cols)
  df = df[, !drop_idx]
  return(df)
}

# ---------------------------------------------------------
# Convert dataframe to a list using a column as list fields.
# ---------------------------------------------------------
df2list = function(df, col_idx, append_list = NULL) {
  if (is.null(append_list)) append_list = list() # Initiate list if needed
  for (i in 1:nrow(df)) {
    vec = unlist(unname(df[i, -col_idx])) # Convert to vector
    vec = na.trim(vec, sides = "right") # Remove trailing NAs
    append_list[[df[i, col_idx]]] = vec # Append vector to list
  }
  return(append_list)
}

# ---------------------------------------------------------
# Create figure window and display plot.
# ---------------------------------------------------------
figure = function(figHandle) {
  operatingSystem = Sys.info()["sysname"]
  if (operatingSystem == "Linux") plot(figHandle)
  if (operatingSystem == "Windows") {
    windows() # Other techniques available
    print(figHandle)
  }
}

# ---------------------------------------------------------
# Platform specific file separator - for readability.
# ---------------------------------------------------------
file_sep = function() {
  platform_file_sep = .Platform$file.sep
  return(platform_file_sep)
}

# ---------------------------------------------------------
# Display blank line in console.
# ---------------------------------------------------------
leaveLine = function(nLines = 1) {
  for (i in 1:nLines) {
    message("")
  }
}

# ---------------------------------------------------------
# Convert a field of a list into a vector.
# ---------------------------------------------------------
list2vec = function(l, idx = NA) {
  vals = unlist(lapply(l, "[[", idx))
  return(vals)
}

# ---------------------------------------------------------
# Number of running and pending jobs on the cluster
# ---------------------------------------------------------
n_slurm_jobs = function(user) {

  # Base sq command for user
  sq = paste("squeue -u", user)

  # Concatenate full commands
  slurm_running = paste(sq, "-t running | wc -l")
  slurm_pending = paste(sq, "-t pending | wc -l")
  slurm_ondemand = paste(sq, "-q interactive | wc -l") # Interactive jobs

  get_jobs_fn = function(x) as.numeric(system(x, intern = TRUE)) - 1
  
  # System call to determine number of slurm processes
  n_running = get_jobs_fn(slurm_running)
  n_pending = get_jobs_fn(slurm_pending)
  n_ondemand = get_jobs_fn(slurm_ondemand)

  # Compile into list
  n_jobs = list(running = n_running - n_ondemand, pending = n_pending)

  return(n_jobs)
}

# ---------------------------------------------------------
# Sum whilst ignoring NAs - better readability.
# ---------------------------------------------------------
nasum = function(x) {
  y = sum(x, na.rm = TRUE)
  return(y)
}

# ---------------------------------------------------------
# Suppress output from a function call.
# ---------------------------------------------------------
quiet = function(x) {
  sinkCon = file("sink.txt")
  sink(sinkCon, type = "output")
  sink(sinkCon, type = "message")
  on.exit(sink(type = "output"))
  on.exit(sink(type = "message"), add = TRUE)
  on.exit(file.remove("sink.txt"), add = TRUE)
  invisible(force(x))
}

# ---------------------------------------------------------
# Concatenate strings with no seperator.
# ---------------------------------------------------------
strcat = function(...) {
  str = paste(..., sep = "")
  return(str)
}

# ---------------------------------------------------------
# Sumbit an array job to the SciCORE cluster. (From AJ Shattock)
# ---------------------------------------------------------
submit_cluster_jobs = function(o, n_jobs, bash_file, ...) {
 
  # Skip if no jobs to run
  if (n_jobs > 0) {
    
    # Check if user is currently running any cluster jobs (see auxiliary.R)
    # check_cluster_jobs(o, action = o$cluster_conflict_action)
    
    # Create a new log file for the cluster jobs (see auxiliary.R)
    log_file = create_bash_log2(paste0(o$pth$log, "/"), o$log_file)
    
    # Construct sbatch array command for running in parallel
    sbatch_array = paste0("--array=1-", n_jobs, "%", o$job_limit)
    
    # Format user-defined options into sbatch-interpretable options
    sbatch_options = c(paste0("--partition=", o$cluster_partition), 
                       paste0("--mem=",       o$job_memory), 
                       paste0("--qos=",       o$job_queue))
    
    # If using scicore partition, include run time option (unlimited on covid19 partition)
    if (o$cluster_partition == "scicore")
      sbatch_options = c(sbatch_options, paste0("--time=", o$job_time))
    
    # Collapse into a whitespace-seperated string
    sbatch_options = paste0(sbatch_options, collapse = " ")
    
    # Extract and format arguments to bash file (including log file)
    bash_inputs = paste(paste0(list(...), collapse = " "), log_file)
    
    # Concatenate system command
    sys_command = paste("sbatch", sbatch_options, sbatch_array, bash_file, bash_inputs)

    # Invoke this command
    system(sys_command)
    
    # Wait for all cluster jobs to complete (see auxiliary.R)
    wait_for_jobs(o, log_file, n_jobs)
  }
}

# ---------------------------------------------------------
# Check for cluster errors, stop if any found
# ---------------------------------------------------------
stop_if_errors = function(pth, err, err_tol = 0, msg = NULL) {
  
  # Set default error message
  if (is.null(msg))
    msg = "Fatal errors when running cluster jobs"
  
  # Check the error file exists
  err_file = paste0(pth, "/", err)
  if (file.exists(err_file)) {
    
    # Load errors from file and convert to vector
    errors = readLines(err_file)
    
    # Stop if any errors, and display them all
    if (length(errors) >= (err_tol + 1))
      stop(msg, " (", length(errors), " errors)\n\n", 
           paste(errors, collapse = "\n"))
  }
}

# ---------------------------------------------------------
# Summing function for higher dimensional arrays.
# ---------------------------------------------------------
sumMat = function(x, sumDim) {
  dims = dim(x)
  nDims = length(dims)
  if (nDims == 2) { # 2D matrices included for readability
    if (sumDim == 1) {
      y = rowSums(x) # Simply apply rowSums
    } else if (sumDim == 2) {
      y = colSums(x) # Simply apply colSums
    }
  }
  if (nDims > 2) { # Higher dimensional arrays are summed and squeezed
    applyDim = setdiff(1:nDims, sumDim)
    y = apply(x, applyDim, sum) # Use apply method
  }
  return(y)
}

# ---------------------------------------------------------
# Format a number with thousand mark separators.
# ---------------------------------------------------------
thou_sep = function(val) {
  format_val = format(val, big.mark = ",") # Automatically rounds
  return(format_val)
}

# ---------------------------------------------------------
# Repeat vector into a 2D matrix.
# ---------------------------------------------------------
vec2mat = function(v, nRow = 1, nCol = 1) {
  m = matrix(v, nrow = nRow, ncol = nCol * length(v), byrow = TRUE)
  return(m)
}

# ---------------------------------------------------------
# Wait until all cluster jobs have finished; by AJ Shattock
# ---------------------------------------------------------
wait_for_jobs = function(o, log_file, n_lines, wait_time = 1, pause_time = 5) {
  
  # Wait for log file to be created
  while (!file.exists(log_file)) Sys.sleep(wait_time)
  
  # Initiate a progress bar
  pb = startProgressBar(n_lines)
  
  # Wait for at least one line to be written
  while (file.info(log_file)$size == 0) Sys.sleep(wait_time)
  
  # ---- Continuously update progress bar ----
  
  # Wait for all jobs to write to log file
  while (nrow(read.table(log_file)) < n_lines) {
    
    # Update progress bar
    setTxtProgressBar(pb, nrow(read.table(log_file)))
    
    # Check number of running and pending jobs
    n_jobs = n_slurm_jobs(user = o$user)
    
    # Have all jobs already finished?
    if (sum(unlist(n_jobs)) == 0) {
      
      # Pause - let any recently completed jobs write to log
      Sys.sleep(pause_time)
      
      break # Break out of while loop
    }
    
    # Wait before testing again
    Sys.sleep(wait_time)
  }
  
  # Finalise progress bar
  setTxtProgressBar(pb, n_lines)
  close(pb)
  
  # ---- Report any batch job failures ----
  
  # Number of jobs not reported in log file
  n_missing = n_lines - nrow(read.table(log_file))
  
  # A problem if this is non-zero
  if (n_missing > 0) {
    
    # Append this error message to an error log file
    err_message = paste0(" ! ", n_missing, " batch job(s) did not finish")
    write(err_message, file = paste0(o$pth$log, "/", o$err_file), append = TRUE)
  }
}

# ---------------------------------------------------------
# Writes a line of text to file (handling file connection).
# ---------------------------------------------------------
write_lines = function(txt, file_name) {
  file_conn = file(file_name, blocking = TRUE)
  writeLines(txt, file_conn)
  close(file_conn)
}

#---------------------------------------------------------------------------
# Initialize progress bar
#---------------------------------------------------------------------------
startProgressBar = function(n_rows) {
  
  bar = txtProgressBar(min = 0, max = n_rows, initial = 0, width = 50, style = 3)
  
  return(bar)
  
}

# ---------------------------------------------------------
# Transform parameters by normalising to a unit cube
# ---------------------------------------------------------
normalise_params = function(df, params, bounds) {
  
  # Normalise all predictor variables
  for (param in params) {
    param_idx = which(params %in% param)
    
    # Use parameter bounds rather than bounds of sampled points
    df[[param]] = normalise_0to1(
      x = df[[param]], 
      x_min = bounds[[param_idx, 1]], 
      x_max = bounds[[param_idx, 2]])
  }
  
  return(df)
}

normalise_0to1 = function(x, x_min = NULL, x_max = NULL, direction = "forward") {
  
  if (!tolower(direction) %in% c("forward", "backward"))
    stop("Normalisation direction must be either 'forward' or 'backward'")
  
  # Forward normalisation
  if (tolower(direction) == "forward") {
    
    # Take bounds from data unless given
    if (is.null(x_min)) x_min = min(x)
    if (is.null(x_max)) x_max = max(x)
    
    # Normalisation equation
    y = (x - x_min) / (x_max - x_min)
    
    # Append original min and max values, needed to backtransform
    attributes(y)$x_min = x_min
    attributes(y)$x_max = x_max
  }
  
  # Backward normalisation
  if (tolower(direction) == "backward") {
    
    # Take bounds from attitubutes of pre-normalised data unless given
    if (is.null(x_min)) x_min = attributes(x)$x_min
    if (is.null(x_max)) x_max = attributes(x)$x_max
    
    # Rearrange equation to solve for x
    #
    # NOTE: as.vector removes all attributes
    y = as.vector(x * (x_max - x_min) + x_min)
  }
  
  return(y)
}

# ---------------------------------------------------------
# Recursively update (nested) lists
# ---------------------------------------------------------
update_list = function(list1, list2) {
  
  for (name in names(list2)) {
    
    if (is.list(list2[[name]]) && is.list(list1[[name]])) {
      # Recursive process if both entries are lists
      list1[[name]] = update_list(list1[[name]], list2[[name]])
      
    } else {
      # Otherwise, overwrite
      list1[[name]] = list2[[name]]
      
    }
  }
  return(list1)
}

# ---------------------------------------------------------
# Custom plotting theme
# ---------------------------------------------------------
theme_amr = function(base_size = 16) {
  
  theme_minimal() +
    theme(
      axis.text = element_text(size = base_size, color = "black"),
      axis.title = element_text(size = base_size, color = "black"),
      legend.text = element_text(size = base_size, color = "black"),
      legend.title = element_text(size = base_size, color = "black"),
      legend.position = "bottom",
      strip.background = element_rect(fill = NA, color = NA),
      strip.text = element_text(size = base_size),
      panel.grid.major  = element_line(color = "white"),
      panel.background = element_rect(fill = "#E7ECEF"),
      plot.background = element_rect(fill = "white", color = "white"))
  
}