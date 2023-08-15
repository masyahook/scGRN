#!/bin/bash

############################################################################
############################################################################
#######                                                              #######
#######  Infer gene regulatory network on cell type aggregated data  #######
#######                                                              #######
############################################################################
############################################################################


##########################################
############# Input params ###############
##########################################

METHOD=$1  # e.g. "grnboost2"|"genie3"
CELL_TYPE=$2  # e.g. "all_data" or "Macrophage" - the aggregated cell type data identifier (cell type-specific or just all data)
PAT_TYPE=$3  # e.g. "all_patients" or ("C"|"M"|"S") - the aggregated patient type data identifier (patient type-specific or just all patients)
NUM_WORKERS=$4  # e.g. 36
Q_THRESH=$5  # e.g. 0.95 - quantile threshold
TF_LIST_PATH=$6  # e.g. /gpfs/projects/bsc08/shared_projects/scGRN_analysis/Data_home/data/TF_lists/lambert2018.txt
RUN_TYPE=$7  # e.g. "GREASY"|"SBATCH" - used to separate logs for different types of computation
TASK_NUM=$8  # e.g. 1 - just number of the tasks (useful to label when searching in logs)

# ----------------------------- #
#### USER-SPECIFIED CONFIGS #####

PROJ_HOME="/gpfs/projects/bsc08/shared_projects/scGRN_analysis"
DATA_ROOT_DIR="${PROJ_HOME}/Data_home/res/covid_19"
DB_NAMES="${PROJ_HOME}/Data_home/data/SCENIC/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather"
MOTIF_ANNOTATION="${PROJ_HOME}/Data_home/data/SCENIC/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"

#################################
# ----------------------------- #

# Creating suffix for for output data file
if [ "$PAT_TYPE" == "all_patients" ]
then
    PAT_SUFFIX=""
else
    PAT_SUFFIX="_${PAT_TYPE}_type"
fi

# Adding _TF suffix if running only on TFs
if [ -z "$TF_LIST_PATH" ]
then
    TF_SUFFIX=
else
    TF_SUFFIX=_TF
fi

# Adding task number label to the logs
if [ -n "$TASK_NUM" ]
then
    TASK_NUM=_${TASK_NUM}
fi

# Defining directory of logs
if [ "$RUN_TYPE" == "GREASY" ]
then
    LOG_OUTPUT="${PROJ_HOME}/sbatch/greasy/logs"
else
    LOG_OUTPUT="${PROJ_HOME}/sbatch/logs"
fi

PATH2DATA=${DATA_ROOT_DIR}/cell_types/${CELL_TYPE}/data/Seurat/raw_data${PAT_SUFFIX}.tsv
OUT_DIR=${DATA_ROOT_DIR}/cell_types/${CELL_TYPE}/data/${METHOD}

GRN_ADJ=${OUT_DIR}/raw_data${PAT_SUFFIX}${TF_SUFFIX}.tsv
GRN_ADJ_COR=${OUT_DIR}/raw_data${PAT_SUFFIX}${TF_SUFFIX}_cor.tsv
GRN_ADJ_CTX=${OUT_DIR}/raw_data${PAT_SUFFIX}${TF_SUFFIX}_ctx.tsv

LOG_OUT=${LOG_OUTPUT}/${CELL_TYPE}${TASK_NUM}${TF_SUFFIX}${PAT_SUFFIX}_${METHOD}_${SLURM_JOBID}.out
LOG_ERR=${LOG_OUTPUT}/${CELL_TYPE}${TASK_NUM}${TF_SUFFIX}${PAT_SUFFIX}_${METHOD}_${SLURM_JOBID}.err

##########################################
############### Running ##################
##########################################

# Creating output folder, logging folder, and folder for temporary dask files
mkdir -p ${DATA_ROOT_DIR}/cell_types/${CELL_TYPE}/data/${METHOD}
mkdir -p ${LOG_OUTPUT}
export DASK_TEMPORARY_DIRECTORY="/tmp/GRN_inference_dask_${SLURM_JOBID}"
mkdir -p $DASK_TEMPORARY_DIRECTORY

# Creating log files
printf "Processing ${CELL_TYPE}${PAT_SUFFIX} data stored in ${PATH2DATA}..\n" > $LOG_OUT
printf "The output will be stored in ${OUT_DIR}..\n\n\n" >> $LOG_OUT
printf "" > $LOG_ERR

# Running either on all genes, or only on TFs
# If running on all genes, using custom script
# If running on TFs, using pyscenic script
if [ -z "$TF_LIST_PATH" ]; then

    printf "No TFs are provided - running $METHOD based on all genes..\n\n" >> $LOG_OUT
    printf "Computing adjacencies..\n" >> $LOG_OUT

    # Run pySCENIC on all genes
    python ../infer_GRN_on_all_genes.py -m $METHOD -i $PATH2DATA -o $GRN_ADJ --num_workers $NUM_WORKERS 2>> $LOG_ERR 1>> $LOG_OUT && \
    printf "Finished computing adjacencies..\n\n" >> $LOG_OUT || \
    printf "Failed computing adjacencies..\n\n" >> $LOG_OUT

    printf "Computing Spearman correlations..\n" >> $LOG_OUT

    # Computing the Spearman correlation for all adjacencies
    pyscenic add_cor $GRN_ADJ $PATH2DATA --output $GRN_ADJ_COR --transpose 2>> $LOG_ERR 1>> $LOG_OUT && \
    python ../post_proc_net_adj_list.py -f $GRN_ADJ_COR -q $Q_THRESH 2>> $LOG_ERR 1>> $LOG_OUT && \
    printf "Finished computing correlations..\n\n" >> $LOG_OUT || \
    printf "Failed computing correlations..\n\n" >> $LOG_OUT

else

    printf "TFs are provided - running $METHOD based on TFs..\n\n" >> $LOG_OUT
    printf "Computing adjacencies..\n" >> $LOG_OUT

    # Run pySCENIC on TF genes
    pyscenic grn $PATH2DATA $TF_LIST_PATH -m $METHOD -o $GRN_ADJ --transpose --num_workers $NUM_WORKERS 2>> $LOG_ERR 1>> $LOG_OUT && \
    printf "Finished computing adjacencies..\n\n" >> $LOG_OUT || \
    printf "Failed computing adjacencies..\n\n" >> $LOG_OUT

    printf "Computing Spearman correlations..\n" >> $LOG_OUT

    # Computing the Spearman correlation for all adjacencies
    pyscenic add_cor $GRN_ADJ $PATH2DATA --output $GRN_ADJ_COR --transpose 2>> $LOG_ERR 1>> $LOG_OUT && \
    python ../post_proc_net_adj_list.py -f $GRN_ADJ_COR -q $Q_THRESH 2>> $LOG_ERR 1>> $LOG_OUT && \
    printf "Finished computing Spearman correlations..\n\n" >> $LOG_OUT || \
    printf "Failed computing Spearman correlations..\n\n" >> $LOG_OUT

    printf "Computing motif-enriched regulons..\n" >> $LOG_OUT

    # Filtering the list of adjacencies by picking only motif-enriched ones
    pyscenic ctx $GRN_ADJ $DB_NAMES --annotations_fname $MOTIF_ANNOTATION --expression_mtx_fname $PATH2DATA --output $GRN_ADJ_CTX --transpose --num_workers $NUM_WORKERS 2>> $LOG_ERR 1>> $LOG_OUT && \
    python ../post_proc_net_adj_list.py -f $GRN_ADJ_CTX -q $Q_THRESH 2>> $LOG_ERR 1>> $LOG_OUT && \
    printf "Finished computing motif-enriched regulons.." >> $LOG_OUT || \
    printf "Failed computing motif-enriched regulons.." >> $LOG_OUT
fi
