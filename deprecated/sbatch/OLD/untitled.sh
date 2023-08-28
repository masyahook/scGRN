#!/bin/bash
#SBATCH --job-name="grnboost2"
#SBATCH --workdir=.
#SBATCH --output=./logs/grnboost2_%j.out
#SBATCH --error=./logs/grnboost2_%j.err
#SBATCH --cpus-per-task=16
#SBATCH --ntasks=1
#SBATCH --time=20:00:00


module load python/3.7.4

export PATH=/home/bsc08/bsc08890/.local/bin:/home/bsc08/bsc08890/bin:$PATH
export PYTHONPATH=/home/bsc08/bsc08890/.local/lib/python3.7/site-packages:$PYTHONPATH

##########################################
############# Input params ###############
##########################################

METHOD=$1  # grnboost2 or genie3
PATH2DATA=~/res/PilotWorkflow/data/Seurat/raw_data.tsv
PATH2TFLIST=~/data/TF_lists/lambert2018.txt
GRN_ADJ=~/res/PilotWorkflow/data/$METHOD/raw_data.tsv
GRN_ADJ_COR=~/res/PilotWorkflow/data/$METHOD/raw_data_cor.tsv

# Pass to the script either 'grnboost2' or 'genie3'
# For example:
# > pyscenic grn grnboost2

##########################################
##########################################
##########################################

pyscenic grn $PATH2DATA $PATH2TFLIST -m $METHOD -o $GRN_ADJ --transpose --num_workers 4
pyscenic add_cor $GRN_ADJ $PATH2DATA --output $GRN_ADJ_COR --transpose