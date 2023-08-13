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
IN_PATH=$2  # e.g. /gpfs/projects/bsc08/shared_projects/res/covid_19/C141/data/Seurat/raw_data.tsv
OUT_PATH=$3  # e.g. /gpfs/projects/bsc08/shared_projects/res/covid_19/C141/data
NUM_WORKERS=$4  # e.g. 36
Q_THRESH=$5  # e.g. 0.95 - quantile threshold
LOG_FOLDER=$6  # e.g. /gpfs/projects/bsc08/shared_projects/sbatch/logs - logging folder
TF_LIST_PATH=$7  # e.g. /gpfs/projects/bsc08/shared_projects/scGRN_analysis/Data_home/data/TF_lists/lambert2018.txt - optional
DB_NAMES=$8  # e.g. /gpfs/projects/bsc08/shared_projects/scGRN_analysis/Data_home/data/SCENIC/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather - optional
MOTIF_ANNOTATION=$9  # e.g. /gpfs/projects/bsc08/shared_projects/scGRN_analysis/Data_home/data/SCENIC/motifs-v9-nr.hgnc-m0.001-o0.0.tbl - optional

# Extracting directory information
OUT_DIR="${OUT_PATH}/${METHOD}"
FILE_NAME="$(basename "$IN_PATH" .tsv)"

# Adding _TF suffix if running only on TFs
if [ -z "$TF_LIST_PATH" ]
then
    TF_SUFFIX=""
else
    TF_SUFFIX="_TF"
fi

# Setting up outputs
GRN_ADJ=${OUT_DIR}/${FILE_NAME}${TF_SUFFIX}.tsv  # gene adjacency list
GRN_ADJ_COR=${OUT_DIR}/${FILE_NAME}${TF_SUFFIX}_cor.tsv  # gene adjacency list with Spearman correlation
GRN_ADJ_CTX=${OUT_DIR}/${FILE_NAME}${TF_SUFFIX}_ctx.tsv  # motif-enriched gene adjacency list

# Setting up logging
LOG_OUT=${LOG_FOLDER}/${FILE_NAME}_${METHOD}${TF_SUFFIX}_${SLURM_JOBID}.out
LOG_ERR=${LOG_FOLDER}/${FILE_NAME}_${METHOD}${TF_SUFFIX}_${SLURM_JOBID}.err

##########################################
############### Running ##################
##########################################

# Creating output folder, and folder for temporary dask files
mkdir -p $OUT_DIR
export DASK_TEMPORARY_DIRECTORY="/tmp/GRN_inference_dask_${SLURM_JOBID}"
mkdir -p $DASK_TEMPORARY_DIRECTORY

# Creating log files
printf "Processing ${IN_PATH} data..\n" > $LOG_OUT
printf "The output will be stored in ${OUT_DIR}..\n\n\n" >> $LOG_OUT
printf "" > $LOG_ERR

# Running either on all genes, or only on TFs
# If running on all genes, using custom script
# If running on TFs, using pyscenic script
if [ -z "$TF_LIST_PATH" ]; then
   
    printf "No TFs are provided - running $METHOD based on all genes..\n\n" >> $LOG_OUT
    printf "Computing adjacencies..\n" >> $LOG_OUT

    # Run pySCENIC on all genes
    python infer_GRN_on_all_genes.py -m $METHOD -i $IN_PATH -o $GRN_ADJ --num_workers $NUM_WORKERS 2>> $LOG_ERR 1>> $LOG_OUT && \
    printf "Finished computing adjacencies..\n\n" >> $LOG_OUT || \
    printf "Failed computing adjacencies..\n\n" >> $LOG_OUT
    
    printf "Computing Spearman correlations..\n" >> $LOG_OUT

    # Computing the Spearman correlation for all adjacencies
    pyscenic add_cor $GRN_ADJ $IN_PATH --output $GRN_ADJ_COR --transpose 2>> $LOG_ERR 1>> $LOG_OUT && \
    python post_proc_net_adj_list.py -f $GRN_ADJ_COR -q $Q_THRESH 2>> $LOG_ERR 1>> $LOG_OUT && \
    printf "Finished computing correlations..\n\n" >> $LOG_OUT || \
    printf "Failed computing correlations..\n\n" >> $LOG_OUT
    
else

    printf "TFs are provided - running $METHOD based on TFs..\n\n" >> $LOG_OUT
    printf "Computing adjacencies..\n" >> $LOG_OUT

    # Run pySCENIC on TF genes
    pyscenic grn $IN_PATH $TF_LIST_PATH -m $METHOD -o $GRN_ADJ --transpose --num_workers $NUM_WORKERS 2>> $LOG_ERR 1>> $LOG_OUT && \
    printf "Finished computing adjacencies..\n\n" >> $LOG_OUT || \
    printf "Failed computing adjacencies..\n\n" >> $LOG_OUT
    
    printf "Computing Spearman correlations..\n" >> $LOG_OUT

    # Computing the Spearman correlation for all adjacencies
    pyscenic add_cor $GRN_ADJ $IN_PATH --output $GRN_ADJ_COR --transpose 2>> $LOG_ERR 1>> $LOG_OUT && \
    python post_proc_net_adj_list.py -f $GRN_ADJ_COR -q $Q_THRESH 2>> $LOG_ERR 1>> $LOG_OUT && \
    printf "Finished computing Spearman correlations..\n\n" >> $LOG_OUT || \
    printf "Failed computing Spearman correlations..\n\n" >> $LOG_OUT
    
    printf "Computing motif-enriched regulons..\n" >> $LOG_OUT

    # Filtering the list of adjacencies by picking only motif-enriched ones
    pyscenic ctx $GRN_ADJ $DB_NAMES --annotations_fname $MOTIF_ANNOTATION --expression_mtx_fname $IN_PATH --output $GRN_ADJ_CTX --transpose --num_workers $NUM_WORKERS 2>> $LOG_ERR 1>> $LOG_OUT && \
    python post_proc_net_adj_list.py -f $GRN_ADJ_CTX -q $Q_THRESH 2>> $LOG_ERR 1>> $LOG_OUT && \
    printf "Finished computing motif-enriched regulons.." >> $LOG_OUT || \
    printf "Failed computing motif-enriched regulons.." >> $LOG_OUT
fi
