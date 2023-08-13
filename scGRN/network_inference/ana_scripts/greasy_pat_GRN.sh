#!/bin/bash

METHOD=$1  # e.g. "grnboost2"|"genie3"
PATIENT=$2  # e.g. C141
PROJ_HOME="/gpfs/projects/bsc08/shared_projects/scGRN_analysis"

##################################################
#                 List of tasks                  #
##################################################
TASK_FILE=${PROJ_HOME}/sbatch/greasy/greasy_tasks_pat_${PATIENT}_${METHOD}

##################################################
#                 Greasy log file                #
##################################################
export GREASY_LOGFILE=${PROJ_HOME}/sbatch/greasy/logs/pat_${PATIENT}_${METHOD}.log

##################################################
#                   Run greasy!                  #
##################################################
mkdir -p ${PROJ_HOME}/sbatch/greasy/logs  # create folder for greasy run
/apps/GREASY/latest/INTEL/IMPI/bin/greasy $TASK_FILE
