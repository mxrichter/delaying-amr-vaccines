###############################################################################
# DIRECTORIES
#
# Create folders and store paths
###############################################################################

# -----------------------------------------------------------------------------
# Define paths for project inputs and outputs
# -----------------------------------------------------------------------------
set_directories = function(o) {

  # Initiate file path lists
  pth = out = list()
  
  # ---- Parent directories ----
  
  # Base path to code repositories
  base_path = file.path("/scicore/home/penny", o$user)

  # Path to project-related code
  pth$code = file.path(base_path, "delaying-amr") # Working directory
  
  # ---- Main folders within project repository ----
  
  # Path to input files (currently empty)
  pth$input = file.path(pth$code, "input")
  
  # Path to output folder
  pth$output = file.path(pth$code, "output")
  
  # Path to publication figures folder
  pth$publication_figures = file.path(pth$output, "publication_figures")
  
  # Path to any files related to configuration of the model
  pth$config = file.path(pth$code, "config")
  
  # Path to output from scicore cluster (log files)
  pth$log_files = file.path(pth$code, "log")

  # Path to OpenMalaria resource files (OM version 45)
  pth$om_files = file.path(base_path, "OpenMalaria", "OM_schema45") 

  # ---- Output directories and files ----

  # Parent path to all output files relating to this project
  output_stem = file.path(pth$output, o$analysis_name)
  
  # Path to parameter tables
  out$parameters = file.path(output_stem, "1_parameters")
  
  # Path to sed patterns to build XML files
  out$sed_patterns = file.path(output_stem, "2_sed-patterns")

  # Path to scenario XML files (input for OpenMalaria)
  out$sim_xml = file.path(output_stem, "3_xml-files")

  # Path to OpenMalaria raw outputs
  out$raw_out = file.path(output_stem, "4_raw-output")

  # Path to post processed outputs
  out$processed_out = file.path(output_stem, "5_processed-output")

  # Path to summarised output
  out$summarised = file.path(output_stem, "6_summarised-output")

  # Path to the optimized Gaussian process
  out$gp_models = file.path(output_stem, "7_gp_models")
  
  # Path to manuscript figures
  out$manuscript = file.path(output_stem, "8_manuscript")

  # Paths to config files
  pth$my_options = file.path(pth$code, "config", "my_options.yaml")
  
  # Make all output directories
  make_output_directories(out) # see function below

  # Append paths to o list
  o = append_directories(o, pth, out) # see function below

  return(o)
}

# -----------------------------------------------------------------------------
# Create output (sub)directories if they do not already exist
# -----------------------------------------------------------------------------
make_output_directories = function(out) {

  # Extract all path names in list
  pth_names = names(out)

  # Loop through these path names
  for (pth_name in pth_names) {
    this_pth = out[[pth_name]]

    # If it does not already exist, create it
    if (!dir.exists(this_pth)) {
      dir.create(this_pth, recursive = TRUE)
    }
  }
}

# -----------------------------------------------------------------------------
# Concatenate separators and append directories to list
# -----------------------------------------------------------------------------
append_directories = function(o, pth, out) {
  
  # Extract all path names in list
  pth_names = names(out)

  # Loop through these path names
  for (pth_name in pth_names) {
    this_pth = out[[pth_name]]

    # Add a file separator to end of output paths
    out[[pth_name]] = paste0(this_pth, file_sep())
  }

  # Concatenate lists
  o$pth = c(pth, out)

  return(o)
}