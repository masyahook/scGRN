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
OUT=$2  # e.g. /gpfs/projects/bsc08/shared_projects/scGRN_analysis/Data_home/res/covid_19 - mandatory parameter
REGULON=$3  # e.g. pyscenic|dorothea|pyscenic_dorothea
Q_THRESH=$4  # e.g. 0.95 - quantile threshold
PLEIOT_CORR=$5  # e.g. T|F
N_PROC=$6  # e.g. 6
SOBJ=$7  # e.g. T|F
VERB=$8  # e.g. T|F

##########################################
############### Running ##################
##########################################

# remembering the current directory, move to the directory of the scripts
CURRENT_DIR=$(pwd)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
cd $SCRIPT_DIR

# running
Rscript run_viper.R -m "$META" -o "$OUT" -a "$ANNO" \
  --annotation_folder "$ANNO_F" -n "$N_PROC" -s "$SOBJ" \
  --add_figure_objects "$SFOBJ" -v "$VERB" && \
Rscript run_viper_merged.R -m "$META" -o "$OUT" -n "$N_PROC" \
  --use_merged_precomputed "$PRE_MERGED" -s "$SOBJ" -v "$VERB" && \
python all_viper_to_pickle.py  # FIXME

# moving back
cd $CURRENT_DIR
