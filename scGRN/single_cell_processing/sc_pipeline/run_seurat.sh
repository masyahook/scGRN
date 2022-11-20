#!/bin/bash

######################################################
######################################################
#######                                        #######
#######  Run Seurat pipeline for all patients  #######
#######                                        #######
######################################################
######################################################

##########################################
############# Input params ###############
##########################################

META=$1  # e.g. ../sc_metadata.tsv - mandatory parameter
OUT=$2  # e.g. /gpfs/projects/bsc08/shared_projects/scGRN_analysis/Data_home/res/covid_19 - mandatory parameter
ANNO=$3  # e.g. HumanPrimaryCellAtlasData
ANNO_F=$4  # e.g. /gpfs/projects/bsc08/shared_projects/scGRN_analysis/Data_home/data/SingleR
N_PROC=$5  # e.g. 6
SOBJ=$6  # e.g. T|F
SFOBJ=$7  # e.g. T|F
VERB=$8  # e.g. T|F

##########################################
############### Running ##################
##########################################

# remembering the current directory, move to the directory of the scripts
CURRENT_DIR=$(pwd)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
cd $SCRIPT_DIR

# running
Rscript run_seurat.R -m "$META" -o "$OUT" -a "$ANNO" \
  --annotation_folder "$ANNO_F" -n "$N_PROC" -s "$SOBJ" \
  --add_figure_objects "$SFOBJ" -v "$VERB"

# moving back
cd $CURRENT_DIR