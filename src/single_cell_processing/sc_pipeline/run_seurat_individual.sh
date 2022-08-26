#!/bin/bash

##########################################
############# Input params ###############
##########################################

PAT=$1  # e.g. C51 - mandatory parameter
META=$2  # e.g. ../sc_metadata.tsv - mandatory parameter
OUT=$3  # e.g. /gpfs/projects/bsc08/bsc08890/res/covid_19 - mandatory parameter
ANNO=$4  # e.g. HumanPrimaryCellAtlasData
ANNO_F=$5  # e.g. /gpfs/projects/bsc08/bsc08890/data/SingleR
SOBJ=$6  # e.g. T|F
SFOBJ=$7  # e.g. T|F
VERB=$8  # e.g. T|F

##########################################
############### Running ##################
##########################################

Rscript run_seurat_individual.R -p "$PAT" -m "$META" -o "$OUT" \
  -a "$ANNO" --annotation_folder "$ANNO_F" -s "$SOBJ" \
  --add_figure_objects "$SFOBJ" -v "$VERB"
