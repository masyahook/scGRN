#!/bin/bash

##########################################################
##########################################################
#######                                            #######
#######  Run full single-cell processing pipeline  #######
#######                                            #######
##########################################################
##########################################################

##########################################
############# Input params ###############
##########################################

META=$1  # e.g. ../sc_metadata.tsv - mandatory parameter
OUT=$2  # e.g. /gpfs/projects/bsc08/shared_projects/scGRN_analysis/Data_home/res/covid_19 - mandatory parameter
ANNO=$3  # e.g. HumanPrimaryCellAtlasData
ANNO_F=$4  # e.g. /gpfs/projects/bsc08/shared_projects/scGRN_analysis/Data_home/data/SingleR
N_PROC=$5  # e.g. 6
PRE_MERGED=$6  # T|F
SOBJ=$7  # e.g. T|F
SFOBJ=$8  # e.g. T|F
VERB=$9  # e.g. T|F

##########################################
############### Running ##################
##########################################

# running
Rscript run_seurat.R -m "$META" -o "$OUT" -a "$ANNO" \
  --annotation_folder "$ANNO_F" -n "$N_PROC" -s "$SOBJ" \
  --save_figure_objects "$SFOBJ" -v "$VERB" && \
Rscript run_seurat_merged.R -m "$META" -o "$OUT" -n "$N_PROC" \
  --use_merged_precomputed "$PRE_MERGED" -s "$SOBJ" -v "$VERB"
