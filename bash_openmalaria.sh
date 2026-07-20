#!/bin/bash
#SBATCH --job-name=amr_openmalaria
#SBATCH --account=penny
#SBATCH --time=06:00:00
#SBATCH --qos=6hours
#SBATCH --mem-per-cpu=2GB
###############################################################################
# Script for submitting OpenMalaria simulations as jobs to the cluster.
#                                                                      
# Input : XML file for each scenario                                   
#                                                                      
# Task: Run each scenario with OpenMalaria on on the cluster            
#                                                                      
# Output: One scenario output file for each simulation                 
#
# Arguments:
#		INPUT_DIR:    Directory containing the scenario XML files
#		DEST_DIR:     Directory where two output files will be created per 
#                 OpenMalaria simulation                                                  
#	    RESOURCE_DIR: Directory of OpenMalaria resources                            
#                                                                               
# Created by,
#   * Monica Golumbeanu (monica.golumbeanu@swisstph.ch),
#   * Andrew James Shattock (andrewjames.shattock@swisstph.ch), and
#   * Thiery Masserey (thiery.masserey@swisstph.ch)
#
# Adapted by Maximilian Richter (maximilian.richter@swisstph.ch)
###############################################################################

# Load OpenMalaria module
module purge

export LMOD_DISABLE_SAME_NAME_AUTOSWAP="no"

ml parallel

ml OpenMalaria/45.0-intel-compilers-2023.1.0

# Define input variables
INPUT_DIR=$1
DEST_DIR=$2
RESOURCE_DIR=$3
LOG_FILE=$4

# Define batch size and parse input as integer
BATCH_SIZE=$5
BATCH_SIZE=$((BATCH_SIZE))

# Extract scenario XML files
SCENARIO_FILES=(${INPUT_DIR}*.xml)

CALLS=()  # Initialize an empty array

# Loop over files and append commands to the text file
for file in "${SCENARIO_FILES[@]}"; do

    SCENARIO_FILE=$file
    
    basename=$(basename "$file")
    
    SCENARIO_NAME="${basename%.xml}"      # Removes .xml: my_file
    
    OUTPUT1=$DEST_DIR$SCENARIO_NAME"_out.txt"
    OUTPUT2=$DEST_DIR$SCENARIO_NAME"_cts.txt"
    
    # Builds function calls
    CALLS+=("openMalaria --scenario $SCENARIO_FILE --output $OUTPUT1")

done

# Calculate the start and end line numbers for this job array task
START=$(( (SLURM_ARRAY_TASK_ID - 1) * BATCH_SIZE))

# Step 3: Slice the CALLS array from START to END
BUNDLE_CALLS=("${CALLS[@]:$START:$BATCH_SIZE}")  # Slice the CALLS array for this bundle

JOBLOG_FILE="/scicore/home/penny/richte0002/delaying-amr/log/scicore_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}_joblog.txt"

cd $RESOURCE_DIR

# Step 4: Run the selected calls in parallel
parallel --joblog "$JOBLOG_FILE" -j $((BATCH_SIZE)) ::: "${BUNDLE_CALLS[@]}"  # Run up to 8 tasks in parallel per bundle

echo "[OpenMalaria] Array task $SLURM_JOB_ID: Job $SLURM_ARRAY_TASK_ID finished" >> "$LOG_FILE"
