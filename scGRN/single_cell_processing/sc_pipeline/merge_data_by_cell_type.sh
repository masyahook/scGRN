#!/bin/bash

#################################################
#################################################
#######                                   #######
#######  Merge patient data by cell type  #######
#######                                   #######
#################################################
#################################################

##########################################
############# Input params ###############
##########################################

META=$1  # e.g. ../sc_metadata.tsv - mandatory parameter
OUT=$2  # e.g. /gpfs/projects/bsc08/bsc08890/res/covid_19 - mandatory parameter
N_PROC=$3  # e.g. 6
PRE_MERGED=$4  # e.g. T|F
SOBJ=$5  # e.g. T|F
VERB=$6  # e.g. T|F

##########################################
############### Running ##################
##########################################

# remembering the current directory, move to the directory of the scripts
CURRENT_DIR=$(pwd)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
cd $SCRIPT_DIR

# running
Rscript merge_data_by_cell_type.R -m "$META" -o "$OUT" -n "$N_PROC" \
  --use_merged_precomputed "$PRE_MERGED" -s "$SOBJ" -v "$VERB"

# moving back
cd $CURRENT_DIR