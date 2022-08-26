#!/bin/bash
#SBATCH --job-name="merge_data_by_cell_type"
#SBATCH --workdir=.
#SBATCH --output=./logs/merge_data_by_cell_type_%j.out
#SBATCH --error=./logs/merge_data_by_cell_type_%j.err
#SBATCH --ntasks=1
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=48


# Loading the modules
module load python/3.7.4

export PATH=/home/bsc08/bsc08890/.local/bin:/home/bsc08/bsc08890/bin:$PATH
export PYTHONPATH=/home/bsc08/bsc08890/.local/lib/python3.7/site-packages:$PYTHONPATH

python /gpfs/home/bsc08/bsc08890/src/merge_data_by_cell_type.py