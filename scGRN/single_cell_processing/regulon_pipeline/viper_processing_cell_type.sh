#!/bin/bash
#SBATCH --job-name="dorothea_viper_processing"
#SBATCH --workdir=.
#SBATCH --output=./logs/dorothea_viper_processing_%j.out
#SBATCH --error=./logs/dorothea_viper_processing_%j.err
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --constraint=highmem
#SBATCH --qos=debug

module load gcc/7.2.0 mkl intel/2018.4 python/3.7.4 graphviz
module load pcre2/10.39 hdf5 szip R/4.1.2

export PATH=/home/bsc08/bsc08890/.local/bin:/home/bsc08/bsc08890/bin:$PATH
export PYTHONPATH=/home/bsc08/bsc08890/.local/lib/python3.7/site-packages:$PYTHONPATH

Rscript ~/R_pipeline/viper_processing_cell_type.R && \ 
python /gpfs/home/bsc08/bsc08890/src/all_viper_to_pickle.py -a True --regulon pyscenic_dorothea && \
echo "Finished running Viper.." || \
echo "Failed running Viper.."