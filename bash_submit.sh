#!/bin/bash

#SBATCH --job-name=amr_postprocess
#SBATCH --account=penny
#SBATCH --cpus-per-task=1
#SBATCH --output=/dev/null
#SBATCH --error=log/scicore_%A_%a.txt

############################################################
# BASH SUBMIT
#
# Submit array of cluster tasks.
#
# Arguments:
#  JOB_TYPE: Type of job we're submitting - see cluster_jobs.R for options
#  LOG_FILE: Text log file that stores IDs of completed tasks
#
# Various sbatch options can be set in options.R (see 'cluster settings' 
# section). For example, cluster partition, job memory, and slurm queue.
#
# Written by A.J.Shattock
# Adapted by Maximilian Richter (maximilian.richter@swisstph.ch)
############################################################

# Load R
module purge
ml R/4.3.2-foss-2023a

# Extract inputs
job_type=$1
sample_id=$2
log_file=$3

# Extract array job ID
job_id=$(expr ${SLURM_ARRAY_TASK_ID})

# Run R script which calls simulation function
Rscript submit.R $job_type $job_id $sample_id

# Write job ID to log file
echo "[Postprocessing] Array task $SLURM_JOB_ID: Job $SLURM_ARRAY_TASK_ID finished ($job_id, $sample_id)" >> "$log_file"
