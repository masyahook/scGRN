#!/bin/bash
##SBATCH --job-name="GREASY_grnboost2"
##SBATCH --workdir=.
##SBATCH --output=./logs/GREASY_grnboost2_%j.out
##SBATCH --error=./logs/GREASY_grnboost2_%j.err
##SBATCH --cpus-per-task=12
##SBATCH --ntasks=18
##SBATCH --time=24:00:00

PATIENT=$1
METHOD=$2

##################################################
#                 List of tasks                  #  
##################################################
FILE=/gpfs/projects/bsc08/bsc08890/sbatch/greasy/greasy_tasks_community_ana_${PATIENT}_${METHOD}

##################################################
#                 Greasy log file                #  
##################################################
export GREASY_LOGFILE=/gpfs/projects/bsc08/bsc08890/sbatch/greasy/logs/commmunity_ana_${PATIENT}_${METHOD}.log

##################################################
#                   Run greasy!                  #  
##################################################
/apps/GREASY/latest/INTEL/IMPI/bin/greasy $FILE