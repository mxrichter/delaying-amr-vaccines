################################################################################
# THE IMPACT OF VACCINES ON DRUG-RESISTANT PLASMODIUM FALCIPARUM
#
# Runs on R 4.3.2 and OpenMalaria v45
#
# Author: Maximilian Richter (maximilian.richter@swisstph.ch)
#
# Code adapted from, 
#  * Thiery Masserey (thiery.masserey@swisstph.ch), and
#  * Andrew James Shattock (andrew.shattock@thekids.org.au)
#
################################################################################

# Clear global environment
rm(list = ls()) 

# Set working directory
setwd(getSrcDirectory(function(){})[1])

# Source R scripts and libraries
source("dependencies.R") 

# Tidy up
if (interactive()) clc() # Clear console
if (interactive()) clf() # Close figures

message("[ THE IMPACT OF VACCINES ON DRUG-RESISTANT P. FALCIPARUM ]")
  
# Initialize simulation options
o = set_options() # see options.R

# Kick off OpenMalaria simulation workflow
# run_simulations(o) # see simulations.R

# Assess final model emulator performance
# assess_emulator_performance(o) # see emulator.R

# Run sensitivity analysis
# sensitivity_analysis(o) # see sensitivity_analysis.R

# Visualize results
prepare_manuscript(o) # see manuscript.R

message("[ Finished ]")