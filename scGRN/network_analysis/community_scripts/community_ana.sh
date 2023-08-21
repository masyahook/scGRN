#!/bin/bash
##SBATCH --job-name="pySCENIC"
##SBATCH --workdir=.
##SBATCH --output=./logs/pySCENIC_%j.out
##SBATCH --error=./logs/pySCENIC_%j.err
##SBATCH --ntasks=1


# Loading the modules
module load python/3.7.4

export PATH=/home/bsc08/bsc08890/.local/bin:/home/bsc08/bsc08890/bin:$PATH
export PYTHONPATH=/home/bsc08/bsc08890/.local/lib/python3.7/site-packages:$PYTHONPATH

##########################################
############# Input params ###############
##########################################

PATIENT=$1  # C51
DATA=$2  # raw_data
ALGO=$3  # either `leiden` or `louvain`
RUN_TYPE=$4
TASK_NUM=$5

# Adding task number label to the logs
if [ ! -z "$TASK_NUM" ]
then
    TASK_NUM=_${TASK_NUM}
fi

# Defining directory of logs
if [ "$RUN_TYPE" == "GREASY" ]
then
    LOG_OUTPUT="/gpfs/projects/bsc08/bsc08890/sbatch/greasy/logs"
else
    LOG_OUTPUT="/gpfs/projects/bsc08/bsc08890/sbatch/logs"
fi

LOG_OUT=${LOG_OUTPUT}/${PATIENT}${TASK_NUM}_${DATA}_COMMUNITY_ANA_${ALGO}_${SLURM_JOBID}.out
LOG_ERR=${LOG_OUTPUT}/${PATIENT}${TASK_NUM}_${DATA}_COMMUNITY_ANA_${ALGO}_${SLURM_JOBID}.err

##########################################
##########################################
##########################################

# Creating log files
printf "Performing community analysis ${PATIENT} (${DATA}) data with ${ALGO} algorithm..\n\n\n" > $LOG_OUT
printf "" > $LOG_ERR

# Running community analysis based on input parameters
python /gpfs/home/bsc08/bsc08890/src/community_ana.py -p $PATIENT -d $DATA -a $ALGO 2> $LOG_ERR 1> $LOG_OUT && \
printf "Finished with community analysis..\n\n" >> $LOG_OUT || \
printf "Failed community analysis..\n\n" >> $LOG_OUT

