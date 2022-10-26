#!/bin/bash

METHOD=$1  # e.g. "grnboost2"|"genie3"
CELL_TYPE=$2  # e.g. "all_data" or "Macrophage" - the aggregated cell type data identifier (cell type-specific or just all data)
PROJ_HOME="/gpfs/projects/bsc08/shared_projects/scGRN_analysis"

##################################################
#                 List of tasks                  #
##################################################
TASK_FILE=${PROJ_HOME}/sbatch/greasy/greasy_tasks_agg_${CELL_TYPE}_${METHOD}

##################################################
#                 Greasy log file                #
##################################################
export GREASY_LOGFILE=${PROJ_HOME}/sbatch/greasy/logs/agg_${CELL_TYPE}_${METHOD}.log

##################################################
#                   Run greasy!                  #
##################################################
mkdir -p ${PROJ_HOME}/sbatch/greasy/logs  # create folder for greasy run
/apps/GREASY/latest/INTEL/IMPI/bin/greasy $TASK_FILE