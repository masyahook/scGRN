#!/bin/bash

ALGO=$1  # either `leiden` or `louvain`
PATIENT=$2  # e.g. C141
PROJ_HOME="/gpfs/projects/bsc08/shared_projects/scGRN_analysis"

##################################################
#                 List of tasks                  #
##################################################
TASK_FILE=${PROJ_HOME}/sbatch/greasy/greasy_tasks_community_ana_pat_${PATIENT}_${ALGO}

##################################################
#                 Greasy log file                #
##################################################
export GREASY_LOGFILE=${PROJ_HOME}/sbatch/greasy/logs/community_ana_pat_${PATIENT}_${ALGO}.log

##################################################
#                   Run greasy!                  #
##################################################
mkdir -p ${PROJ_HOME}/sbatch/greasy/logs  # create folder for greasy run
/apps/GREASY/latest/INTEL/IMPI/bin/greasy $TASK_FILE
