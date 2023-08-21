#!/bin/bash

############################################################
############################################################
#######                                              #######
#######  Run community analysis on one patient data  #######
#######                                              #######
############################################################
############################################################

# Loading the modules
# module load python/3.7.4

# export PATH=/home/bsc08/bsc08890/.local/bin:/home/bsc08/bsc08890/bin:$PATH
# export PYTHONPATH=/home/bsc08/bsc08890/.local/lib/python3.7/site-packages:$PYTHONPATH

##########################################
############# Input params ###############
##########################################

ALGO=$1  # either `leiden` or `louvain`
CELL_TYPE=$2  # e.g. "all_data" or "Macrophage" - the aggregated cell type data identifier (cell type-specific or just all data)
PAT_TYPE=$3  # e.g. "all_patients" or ("C"|"M"|"S") - the aggregated patient type data identifier (patient type-specific or just all patients)
RUN_TYPE=$4  # e.g. "GREASY"|"SBATCH" - used to separate logs for different types of computation
TASK_NUM=$5  # e.g. 1 - just number of the tasks (useful to label when searching in logs)

# ----------------------------- #
#### USER-SPECIFIED CONFIGS #####

PROJ_HOME="/gpfs/projects/bsc08/shared_projects/scGRN_analysis"

#################################
# ----------------------------- #

# Adding task number label to the logs
if [ ! -z "$TASK_NUM" ]
then
    TASK_NUM=_${TASK_NUM}
fi

# Defining directory of logs
if [ "$RUN_TYPE" == "GREASY" ]
then
    LOG_OUTPUT="${PROJ_HOME}/sbatch/greasy/logs"
else
    LOG_OUTPUT="${PROJ_HOME}/sbatch/logs"
fi

LOG_OUT=${LOG_OUTPUT}/${PAT_TYPE}_type_${TASK_NUM}_${CELL_TYPE}_COMMUNITY_ANA_${ALGO}_${SLURM_JOBID}.out
LOG_ERR=${LOG_OUTPUT}/${PAT_TYPE}_type_${TASK_NUM}_${CELL_TYPE}_COMMUNITY_ANA_${ALGO}_${SLURM_JOBID}.err

##########################################
############### Running ##################
##########################################

# Creating log files
printf "Performing community analysis cell type specific, pat type aggregated '${PAT_TYPE}' (${CELL_TYPE}) data with ${ALGO} algorithm..\n\n\n" > $LOG_OUT
printf "" > $LOG_ERR

# Running community analysis based on input parameters
python ../community_ana.py -p $PAT_TYPE -c $CELL_TYPE -a $ALGO 2> $LOG_ERR 1> $LOG_OUT && \
printf "Finished with community analysis..\n\n" >> $LOG_OUT || \
printf "Failed community analysis..\n\n" >> $LOG_OUT
