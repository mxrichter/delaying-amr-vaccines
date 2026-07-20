################################################################################
# POSTPROCESS
# 
# Start cluster job for postprocessing the raw output from OpenMalaria
# (job type 'postprocess', see cluster_jobs.R)
# 
################################################################################

#-------------------------------------------------------------------------------
# Launch function to postprocess the output from OpenMalaria
#-------------------------------------------------------------------------------
postprocess_results = function(o, p, m) {

  # Exit function if postprocessing is not part of current workflow
  if(!(5 %in% o$do_step)) return(NULL)
  
  message("* Processing results")

  # Load list of scenarios
  scenarios_df = readRDS(paste0(h(m)$pth$parameters, "/all_scenarios.rds"))
  
  # List of scenarios
  fids = scenarios_df$fid %>% unique()
  
  # Number of scenarios that are to be processed
  n_settings = length(fids)
  
  # Submit all jobs to the cluster (see auxiliary.R)
  submit_cluster_jobs(o, n_settings, "bash_submit.sh", "postprocess", 
                      m$sample_num)
  
  # Throw an error if any cluster jobs failed (see auxiliary.R)
  stop_if_errors(o$pth$log, o$err_file, err_tol = 1)

}