################################################################################
# RUN OPENMALARIA
#
# Launch file for OpenMalaria simulations
#
################################################################################

# ------------------------------------------------------------------------------
# Run OpenMalaria for all XML scenario files
# ------------------------------------------------------------------------------
run_openmalaria = function(o, p, m) {
 
  # Skip step if not specified in list of options 
  if(!(4 %in% o$do_step)) return(NULL)
  
  message("* Running OpenMalaria simulations") 
  
  # List of XML files in directory
  all_xml_files = list.files(h(m)$pth$sim_xml) 
  
  # List of raw input files
  scenarios_in = all_xml_files %>% sub(".xml", "", .)
  
  # List of raw output files
  scenarios_out = list.files(h(m)$pth$raw_out) %>% sub("_out.txt", "", .)
  
  # Filter which scenarios have already been simulated
  scenarios_simulated = intersect(scenarios_in, scenarios_out)
  
  # List of scenarios that have not yet been simulated
  scenarios_left = setdiff(scenarios_in, scenarios_out)
  
  # Output number of simulated/missing simulations
  message(paste0("    > Number of simulations: ", length(all_xml_files)))
  message(paste0("    > Skipping: ", length(scenarios_in) - length(scenarios_left)))
  
  # If any scenarios have been simulated already, ...
  if(length(scenarios_simulated) > 0) {

    # List of XML files of scenarios that have been simulated already
    scenarios_simulated_xml = file.path(h(m)$pth$sim_xml, paste0(scenarios_simulated, ".xml"))

    # Delete XML file
    unlink(scenarios_simulated_xml)

  }

  # Number of scenarios that are going to be simulated
  n_jobs = length(scenarios_left)

  # If there are any scenarios left, ...
  if(n_jobs > 0) {

    # Create a new log file for the cluster jobs
    log_file = file.path(o$pth$log_files, "scicore_log.txt")
    create_bash_log(log_file)

    # Number of batches that are passed to cluster
    n_batches = ceil(n_jobs/o$batch_size)

    # Construct slurm array command for running in parallel
    slurm_array = paste0("--array=1-", n_batches, "%", o$job_limit, " ",
                         "--cpus-per-task=", o$cpus_per_task, " ",
                         "--output=", o$pth$log_files, "/scicore_%A_%a.txt")

    # Add country sub directory to file paths
    xml_path = paste0(h(m)$pth$sim_xml, "/")
    sim_path = paste0(h(m)$pth$raw_out, "/")

    # Path to OpenMalaria files
    o$pth$om_files = paste0(o$pth$om_files, "/")

    # Concatenate system command
    sys_command = paste("sbatch", slurm_array, "bash_openmalaria.sh",
                        xml_path, sim_path, o$pth$om_files, log_file, o$batch_size)

    # Invoke this command
    system(sys_command)

    # Wait for all cluster jobs to complete
    wait_for_jobs(o, log_file, n_batches)

  }

}