#!/bin/bash

#################################################
#################################################
#######                                   #######
#######  Run full VIPER regulon pipeline  #######
#######                                   #######
#################################################
#################################################

##########################################
############# Input params ###############
##########################################

META=$1  # e.g. ../sc_metadata.tsv - mandatory parameter
META_CTYPE=$2  # e.g. /gpfs/projects/bsc08/shared_projects/scGRN_analysis/Data_home/res/covid_19/cell_type_meta.tsv - mandatory parameter
OUT=$3  # e.g. /gpfs/projects/bsc08/shared_projects/scGRN_analysis/Data_home/res/covid_19 - mandatory parameter
REGULON=$4  # e.g. pyscenic|dorothea|pyscenic_dorothea - regulon type to use
Q_THRESH=$5  # e.g. 0.95 - quantile threshold to use if REGULON=pyscenic
PLEIOT_CORR=$6  # e.g. T|F - apply VIPER pleiotropy correction
N_PROC=$7  # e.g. 6 - number of parallel processes
SOBJ=$8  # e.g. T|F - save as seurat object
VERB=$9  # e.g. T|F - verbosity

##########################################
############### Running ##################
##########################################

# running
Rscript run_viper.R -m "$META" -o "$OUT" -r "$REGULON" -q "$Q_THRESH" \
  -c "$PLEIOT_CORR" -n "$N_PROC" -s "$SOBJ" -v "$VERB" && \
Rscript run_viper_merged.R -m "$META_CTYPE" -o "$OUT" -r "$REGULON" \
  -q "$Q_THRESH" -c "$PLEIOT_CORR" -n "$N_PROC" -s "$SOBJ" -v "$VERB" && \
python all_viper_to_pickle.py -o $OUT -r $REGULON && \
printf "Finished running regulon pipeline..\n\n" || \
printf "Failed running regulon pipeline..\n\n"
