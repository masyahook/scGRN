#!/bin/bash

#########################################################################
#########################################################################
#######                                                           #######
#######  Run VIPER regulon pipeline for on merged cell-type data  #######
#######                                                           #######
#########################################################################
#########################################################################

##########################################
############# Input params ###############
##########################################

META_CTYPE=$1  # e.g. /gpfs/projects/bsc08/shared_projects/scGRN_analysis/Data_home/res/covid_19/cell_type_meta.tsv - mandatory parameter
OUT=$2  # e.g. /gpfs/projects/bsc08/shared_projects/scGRN_analysis/Data_home/res/covid_19 - mandatory parameter
REGULON=$3  # e.g. pyscenic|dorothea|pyscenic_dorothea - regulon type to use
Q_THRESH=$4  # e.g. 0.95 - quantile threshold to use if REGULON=pyscenic
PLEIOT_CORR=$5  # e.g. T|F - apply VIPER pleiotropy correction
N_PROC=$6  # e.g. 6 - number of parallel processes
SOBJ=$7  # e.g. T|F - save as seurat object
VERB=$8  # e.g. T|F - verbosity

##########################################
############### Running ##################
##########################################

# running
Rscript run_viper_merged.R -m "$META_CTYPE" -o "$OUT" -r "$REGULON" -q "$Q_THRESH"
  -c "$PLEIOT_CORR" -n "$N_PROC" -s "$SOBJ" -v "$VERB"
