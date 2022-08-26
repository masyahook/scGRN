#!/bin/bash

##########################################
############# Input params ###############
##########################################

META=$1  # e.g. ../sc_metadata.tsv - mandatory parameter
OUT=$2  # e.g. /gpfs/projects/bsc08/bsc08890/res/covid_19 - mandatory parameter
ANNO=$3  # e.g. HumanPrimaryCellAtlasData
ANNO_F=$4  # e.g. /gpfs/projects/bsc08/bsc08890/data/SingleR
SOBJ=$5  # e.g. T|F
SFOBJ=$6  # e.g. T|F
VERB=$7  # e.g. T|F

##########################################
############### Running ##################
##########################################

Rscript run_seurat.R -m "$META" -o "$OUT" -a "$ANNO" \
  --annotation_folder "$ANNO_F" -s "$SOBJ" --add_figure_objects "$SFOBJ" \
  -v "$VERB""