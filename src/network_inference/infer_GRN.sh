#!/bin/bash

#############################################################
#############################################################
#######                                               #######
#######  Infer gene regulatory network on input data  #######
#######                                               #######
#############################################################
#############################################################

##########################################
############# Input params ###############
##########################################

METHOD=$1  # e.g. "grnboost2"|"genie3"
PATH2DATA=$2  # e.g. /gpfs/projects/bsc08/bsc08890/res/covid_19/C141/data/Seurat/raw_data.tsv
NUM_WORKERS=$3  # e.g. 36
Q_THRESH=$4  # e.g. 0.95 - quantile threshold
LOG_FOLDER=$5  # e.g. /gpfs/projects/bsc08/bsc08890/sbatch/logs - logging folder
PATH2TFLIST=$6  # e.g. /gpfs/projects/bsc08/bsc08890/data/TF_lists/lambert2018.txt
DB_NAMES=$7  # e.g. /gpfs/projects/bsc08/bsc08890/data/SCENIC/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather
MOTIF_ANNOTATION=$8  # e.g. /gpfs/projects/bsc08/bsc08890/data/SCENIC/motifs-v9-nr.hgnc-m0.001-o0.0.tbl
TASK_NUM=$9  # e.g. 1 - just number of the tasks (useful to label when searching in logs)

# Remembering the directories, moving to the directory of the scripts
CURRENT_DIR=$(pwd)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
DATA_DIR="$(dirname "$PATH2DATA")"
OUT_DIR="$(dirname "$DATA_DIR")"/$METHOD
FILE_NAME="$(basename "$PATH2DATA" .tsv)"
cd $SCRIPT_DIR

# Adding _TF suffix if running only on TFs
if [ -z "$PATH2TFLIST" ]
then
    TF_SUFFIX=""
else
    TF_SUFFIX="_TF"
fi

# Adding task number label to the logs
if [ -n "$TASK_NUM" ]
then
    TASK_NUM="_${TASK_NUM}"
fi

# Setting up outputs
GRN_ADJ=${OUT_DIR}/${FILE_NAME}${TF_SUFFIX}.tsv  # gene adjacency list
GRN_ADJ_COR=${OUT_DIR}/${FILE_NAME}${TF_SUFFIX}_cor.tsv  # gene adjacency list with Spearman correlation
GRN_ADJ_CTX=${OUT_DIR}/${FILE_NAME}${TF_SUFFIX}_ctx.tsv  # motif-enriched gene adjacency list

# Setting up logging
LOG_OUT=${LOG_FOLDER}/${FILE_NAME}_${METHOD}${TASK_NUM}${TF_SUFFIX}_${SLURM_JOBID}.out
LOG_ERR=${LOG_FOLDER}/${FILE_NAME}_${METHOD}${TASK_NUM}${TF_SUFFIX}_${SLURM_JOBID}.err

##########################################
############### Running ##################
##########################################

# Creating output folder, and folder for temporary dask files
mkdir -p $OUT_DIR
export DASK_TEMPORARY_DIRECTORY="/tmp/GRN_inference_dask_${SLURM_JOBID}"
mkdir -p $DASK_TEMPORARY_DIRECTORY

# Creating log files
printf "Processing ${PATH2DATA} data..\n" > $LOG_OUT
printf "The output will be stored in ${OUT_DIR}..\n\n\n" >> $LOG_OUT
printf "" > $LOG_ERR

# Running either on all genes, or only on TFs
# If running on all genes, using custom script
# If running on TFs, using pyscenic script
if [ -z "$PATH2TFLIST" ]; then
   
    printf "No TFs are provided - running $METHOD based on all genes..\n\n" >> $LOG_OUT
    printf "Computing adjacencies..\n" >> $LOG_OUT

    # Run pySCENIC on all genes
    python infer_GRN_on_all_genes.py -m $METHOD -i $PATH2DATA -o $GRN_ADJ --num_workers $NUM_WORKERS 2>> $LOG_ERR 1>> $LOG_OUT && \
    printf "Finished computing adjacencies..\n\n" >> $LOG_OUT || \
    printf "Failed computing adjacencies..\n\n" >> $LOG_OUT
    
    printf "Computing Spearman correlations..\n" >> $LOG_OUT

    # Computing the Spearman correlation for all adjacencies
    pyscenic add_cor $GRN_ADJ $PATH2DATA --output $GRN_ADJ_COR --transpose 2>> $LOG_ERR 1>> $LOG_OUT && \
    python post_proc_net_adj_list.py -f $GRN_ADJ_COR -q $Q_THRESH 2>> $LOG_ERR 1>> $LOG_OUT && \
    printf "Finished computing correlations..\n\n" >> $LOG_OUT || \
    printf "Failed computing correlations..\n\n" >> $LOG_OUT
    
else

    printf "TFs are provided - running $METHOD based on TFs..\n\n" >> $LOG_OUT
    printf "Computing adjacencies..\n" >> $LOG_OUT

    # Run pySCENIC on TF genes
    pyscenic grn $PATH2DATA $PATH2TFLIST -m $METHOD -o $GRN_ADJ --transpose --num_workers $NUM_WORKERS 2>> $LOG_ERR 1>> $LOG_OUT && \
    printf "Finished computing adjacencies..\n\n" >> $LOG_OUT || \
    printf "Failed computing adjacencies..\n\n" >> $LOG_OUT
    
    printf "Computing Spearman correlations..\n" >> $LOG_OUT

    # Computing the Spearman correlation for all adjacencies
    pyscenic add_cor $GRN_ADJ $PATH2DATA --output $GRN_ADJ_COR --transpose 2>> $LOG_ERR 1>> $LOG_OUT && \
    python post_proc_net_adj_list.py -f $GRN_ADJ_COR -q $Q_THRESH 2>> $LOG_ERR 1>> $LOG_OUT && \
    printf "Finished computing Spearman correlations..\n\n" >> $LOG_OUT || \
    printf "Failed computing Spearman correlations..\n\n" >> $LOG_OUT
    
    printf "Computing motif-enriched regulons..\n" >> $LOG_OUT

    # Filtering the list of adjacencies by picking only motif-enriched ones
    pyscenic ctx $GRN_ADJ $DB_NAMES --annotations_fname $MOTIF_ANNOTATION --expression_mtx_fname $PATH2DATA --output $GRN_ADJ_CTX --transpose --num_workers $NUM_WORKERS 2>> $LOG_ERR 1>> $LOG_OUT && \
    python post_proc_net_adj_list.py -f $GRN_ADJ_CTX -q $Q_THRESH 2>> $LOG_ERR 1>> $LOG_OUT && \
    printf "Finished computing motif-enriched regulons.." >> $LOG_OUT || \
    printf "Failed computing motif-enriched regulons.." >> $LOG_OUT
fi

# moving back
cd $CURRENT_DIR
