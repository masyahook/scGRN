#!/bin/bash
#SBATCH --job-name="all_adj_lists_to_process"
#SBATCH --workdir=.
#SBATCH --output=./logs/all_adj_lists_to_process_%j.out
#SBATCH --error=./logs/all_adj_lists_to_process_%j.err
#SBATCH --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=48
#SBATCH --qos=debug


# Loading the modules
module load python/3.7.4

export PATH=/home/bsc08/bsc08890/.local/bin:/home/bsc08/bsc08890/bin:$PATH
export PYTHONPATH=/home/bsc08/bsc08890/.local/lib/python3.7/site-packages:$PYTHONPATH

python /gpfs/home/bsc08/bsc08890/src/all_adj_lists_to_process.py