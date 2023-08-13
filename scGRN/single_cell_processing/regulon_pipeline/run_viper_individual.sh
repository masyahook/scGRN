#!/bin/bash

############################################################
############################################################
#######                                              #######
#######  Run VIPER regulon pipeline for one patient  #######
#######                                              #######
############################################################
############################################################

##########################################
############# Input params ###############
##########################################

PAT=$1  # e.g. C51 - mandatory parameter
META=$2  # e.g. ../sc_metadata.tsv - mandatory parameter
OUT=$3  # e.g. /gpfs/projects/bsc08/shared_projects/scGRN_analysis/Data_home/res/covid_19 - mandatory parameter
REGULON=$4  # e.g. pyscenic|dorothea|pyscenic_dorothea - regulon type to use
Q_THRESH=$5  # e.g. 0.95 - quantile threshold - quantile threshold to use if REGULON=pyscenic
PLEIOT_CORR=$6  # e.g. T|F - apply VIPER pleiotropy correction
N_PROC=$7  # e.g. 6 - number of parallel processes
SOBJ=$8  # e.g. T|F - save as seurat object
VERB=$9  # e.g. T|F - verbosity

##########################################
############### Running ##################
##########################################

# remembering the current directory, move to the directory of the scripts
CURRENT_DIR=$(pwd)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
cd $SCRIPT_DIR

# running
Rscript run_viper_individual.R -p "$PAT" -m "$META" -o "$OUT" \
  -r "$REGULON" -q "$Q_THRESH" -c "$PLEIOT_CORR" -n "$N_PROC" -s "$SOBJ" \
  -v "$VERB"
  
# moving back
cd $CURRENT_DIR