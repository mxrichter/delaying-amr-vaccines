###############################################################################
# DEPENDENCIES
#                                             
# Load all necessary R scripts and R libraries
#
###############################################################################

# Load R libraries
library(tgp)
library(hetGP)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(pracma)
library(stringr)
library(xlsx)
library(rlist)
library(zoo)
library(data.table)
library(sensitivity)
library(wrapr)
library(ggplot2)
library(reshape2)
library(yaml)
library(gridExtra)
library(lubridate)
library(ggpubr)
library(reticulate)
library(ggtext)

# Source R scripts
source("directories.R")
source("auxiliary.R")
source("simulations.R")
source("openmalaria_setup.R")
source("cluster_jobs.R")
source("openmalaria_run.R")
source("postprocess.R")
source("summarize_results.R")
source("options.R")
source("parse_parameters.R")
source("adaptive_sampling.R")
source("sensitivity_analysis.R")
source("emulator.R")
source("manuscript.R")