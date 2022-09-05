#!/bin/bash

###################################################################
###################################################################
#######                                                     #######
#######  Infer gene regulatory network on one patient data  #######
#######                                                     #######
###################################################################
###################################################################


##########################################
############# Input params ###############
##########################################

METHOD=$1  # e.g. "grnboost2"|"genie3"
PATIENT=$2  # e.g. C141
CELL_TYPE=$3  # e.g. 'all_data' or 'Macrophage' - the cell type data identifier (cell type-specific or just all data)
NUM_WORKERS=$4  # e.g. 36
Q_THRESH=$5  # e.g. 0.95 - quantile threshold
PATH2TFLIST=$6  # e.g. /gpfs/projects/bsc08/bsc08890/data/TF_lists/lambert2018.txt
RUN_TYPE=$7  # e.g. "GREASY"|"SBATCH" - used to separate logs for different types of computation
TASK_NUM=$8  # e.g. 1 - just number of the tasks (useful to label when searching in logs)

# Remembering the directories, moving to the directory of the scripts
CURRENT_DIR=$(pwd)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
cd $SCRIPT_DIR

# Creating suffix for accessing the correct data file
if [ -z "$CELL_TYPE" ]
then
    CELL_TYPE_SUFFIX=""
else
    CELL_TYPE_SUFFIX="_${CELL_TYPE}"
fi

# Adding _TF suffix if running only on TFs
if [ -z "$PATH2TFLIST" ]
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
    LOG_OUTPUT="/gpfs/projects/bsc08/bsc08890/sbatch/greasy/logs"
else
    LOG_OUTPUT="/gpfs/projects/bsc08/bsc08890/sbatch/logs"
fi
# ----------------------------- #
#### USER-SPECIFIED CONFIGS #####

DATA_ROOT_DIR="/gpfs/projects/bsc08/bsc08890/res/covid_19"
DB_NAMES="/gpfs/projects/bsc08/bsc08890/data/SCENIC/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather"
MOTIF_ANNOTATION="/gpfs/projects/bsc08/bsc08890/data/SCENIC/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"

#################################
# ----------------------------- #

PATH2DATA=${DATA_ROOT_DIR}/${PATIENT}/data/Seurat/raw_data${CELL_TYPE_SUFFIX}.tsv
OUT_DIR=${DATA_ROOT_DIR}/${PATIENT}/data/${METHOD}

GRN_ADJ=${OUT_DIR}/raw_data${CELL_TYPE_SUFFIX}${TF_SUFFIX}.tsv
GRN_ADJ_COR=${OUT_DIR}/raw_data${CELL_TYPE_SUFFIX}${TF_SUFFIX}_cor.tsv
GRN_ADJ_CTX=${OUT_DIR}/raw_data${CELL_TYPE_SUFFIX}${TF_SUFFIX}_ctx.tsv

LOG_OUT=${LOG_OUTPUT}/${PATIENT}${TASK_NUM}${TF_SUFFIX}${CELL_TYPE_SUFFIX}_pySCENIC_${METHOD}_${SLURM_JOBID}.out
LOG_ERR=${LOG_OUTPUT}/${PATIENT}${TASK_NUM}${TF_SUFFIX}${CELL_TYPE_SUFFIX}_pySCENIC_${METHOD}_${SLURM_JOBID}.err

##########################################
############### Running ##################
##########################################

# Creating output folder
mkdir -p ${DATA_ROOT_DIR}/${PATIENT}/data/${METHOD}

# Creating log files
printf "Processing ${PATIENT}${CELL_TYPE_SUFFIX} data stored in ${PATH2DATA}..\n" > $LOG_OUT
printf "The output will be stored in ${OUT_DIR}..\n\n\n" >> $LOG_OUT
printf "" > $LOG_ERR

# Running either on all genes, or only on TFs
# If running on all genes, using custom script
# If running on TFs, using pyscenic
if [ -z "$PATH2TFLIST" ]; then
   
    printf "No TFs are provided - running $METHOD based on all genes..\n\n" >> $LOG_OUT
    printf "Computing adjacencies..\n" >> $LOG_OUT

    # Run pySCENIC on all genes
    python ../run_pyscenic_on_all_genes.py -m $METHOD -i $PATH2DATA -o $GRN_ADJ --num_workers $NUM_WORKERS 2> $LOG_ERR 1> $LOG_OUT && \
    printf "Finished computing adjacencies..\n\n" >> $LOG_OUT || \
    printf "Failed computing adjacencies..\n\n" >> $LOG_OUT
    
    printf "Computing Spearman correlations..\n" >> $LOG_OUT

    # Computing the Spearman correlation for all adjacencies
    pyscenic add_cor $GRN_ADJ $PATH2DATA --output $GRN_ADJ_COR --transpose 2> $LOG_ERR 1> $LOG_OUT && \
    python ../post_proc_net_adj_list.py -f $GRN_ADJ_COR -q $Q_THRESH 2>> $LOG_ERR 1>> $LOG_OUT && \
    printf "Finished computing corelations..\n\n" >> $LOG_OUT || \
    printf "Failed computing correlations..\n\n" >> $LOG_OUT
    
else

    printf "TFs are provided - running $METHOD based on TFs..\n\n" >> $LOG_OUT
    printf "Computing adjacencies..\n" >> $LOG_OUT

    # Run pySCENIC on TF genes
    pyscenic grn $PATH2DATA $PATH2TFLIST -m $METHOD -o $GRN_ADJ --transpose --num_workers $NUM_WORKERS 2> $LOG_ERR 1> $LOG_OUT && \
    printf "Finished computing adjacencies..\n\n" >> $LOG_OUT || \
    printf "Failed computing adjacencies..\n\n" >> $LOG_OUT
    
    printf "Computing Spearman correlations..\n" >> $LOG_OUT

    # Computing the Spearman correlation for all adjacencies
    pyscenic add_cor $GRN_ADJ $PATH2DATA --output $GRN_ADJ_COR --transpose 2> $LOG_ERR 1> $LOG_OUT && \
    python ../post_proc_net_adj_list.py -f $GRN_ADJ_COR -q $Q_THRESH 2>> $LOG_ERR 1>> $LOG_OUT && \
    printf "Finished computing Spearman correlations..\n\n" >> $LOG_OUT || \
    printf "Failed computing Spearman correlations..\n\n" >> $LOG_OUT
    
    printf "Computing motif-enriched regulons..\n" >> $LOG_OUT

    # Filtering the list of adjacencies by picking only motif-enriched ones
    pyscenic ctx $GRN_ADJ $DB_NAMES --annotations_fname $MOTIF_ANNOTATION --expression_mtx_fname $PATH2DATA --output $GRN_ADJ_CTX --transpose --num_workers $NUM_WORKERS 2> $LOG_ERR 1> $LOG_OUT && \
    python ../post_proc_net_adj_list.py -f $GRN_ADJ_CTX -q $Q_THRESH 2>> $LOG_ERR 1>> $LOG_OUT && \
    printf "Finished computing motif-enriched regulons.." >> $LOG_OUT || \
    printf "Failed computing motif-enriched regulons.." >> $LOG_OUT
fi

# moving back
cd $CURRENT_DIR
