################################################################################
# SUBMIT
#
# Submit cluster tasks from a cluster array.
#
################################################################################

suppressPackageStartupMessages(source("dependencies.R"))

# Extract input from bash file
args = commandArgs(trailingOnly = TRUE)

# Name arguments when provided from bash
if (length(args) > 0) {
  job_type = as.character(args[1])
  job_id   = as.numeric(args[2])
  
  if(length(args) > 2) {
    sample_id = as.numeric(args[3])
  } else {
    sample_id = NULL
  }
    
} else {  # Otherwise define some defaults - used for testing and debugging

  job_type  = "postprocess" # Options: postprocess
  job_id    = 1
  sample_id = 0 # Only used when postprocessing samples
  
}

# Reload simulation options (see options.R)
o = set_options()

# Run task defined by job_type for job_id and catch any error
tryCatch(
  run_cluster_job(o, job_type, job_id, sample_id),  # See cluster_jobs.R
  
  # Error handler
  error = function(e) {
    
    # Concatenate (ideally useful) error message
    err_message = paste0(" ! ", e$message, " (array ID: ", job_id, ")")
    
    message(err_message)
    
    # Append this error message to an error log file
    write(err_message, file = paste0(o$pth$log,"/", o$err_file), append = TRUE)
  },
  
  # Close up
  finally = message("Closing cluster job")
)