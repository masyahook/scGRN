#!/bin/bash

#####################################################
#####################################################
#######                                       #######
#######  Run Seurat pipeline for one patient  #######
#######                                       #######
#####################################################
#####################################################

##########################################
############# Input params ###############
##########################################

PAT=$1  # e.g. C51 - mandatory parameter
META=$2  # e.g. ../sc_metadata.tsv - mandatory parameter
OUT=$3  # e.g. /gpfs/projects/bsc08/shared_projects/scGRN_analysis/Data_home/res/covid_19 - mandatory parameter
ANNO=$4  # e.g. HumanPrimaryCellAtlasData
ANNO_F=$5  # e.g. /gpfs/projects/bsc08/shared_projects/scGRN_analysis/Data_home/data/SingleR
SOBJ=$6  # e.g. T|F
SFOBJ=$7  # e.g. T|F
VERB=$8  # e.g. T|F

##########################################
############### Running ##################
##########################################

# running
Rscript run_seurat_individual.R -p "$PAT" -m "$META" -o "$OUT" \
  -a "$ANNO" --annotation_folder "$ANNO_F" -s "$SOBJ" \
  --save_figure_objects "$SFOBJ" -v "$VERB"
