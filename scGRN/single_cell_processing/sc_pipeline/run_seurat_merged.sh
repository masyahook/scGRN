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
OUT=$2  # e.g. /gpfs/projects/bsc08/shared_projects/scGRN_analysis/Data_home/res/covid_19 - mandatory parameter
N_PROC=$3  # e.g. 6
# first step of pipeline is data merging, result will be saved and could be re-used for future runs
# the flag below indicates whether to skip the first step
PRE_MERGED=$4  # e.g. T|F
SOBJ=$5  # e.g. T|F
VERB=$6  # e.g. T|F

##########################################
############### Running ##################
##########################################

# running
Rscript run_seurat_merged.R -m "$META" -o "$OUT" -n "$N_PROC" \
  --use_merged_precomputed "$PRE_MERGED" -s "$SOBJ" -v "$VERB"
