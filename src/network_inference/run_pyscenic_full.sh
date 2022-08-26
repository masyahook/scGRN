#!/bin/bash
##SBATCH --job-name="pySCENIC"
##SBATCH --workdir=.
##SBATCH --output=./logs/pySCENIC_%j.out
##SBATCH --error=./logs/pySCENIC_%j.err
##SBATCH --ntasks=1


# Loading the modules
module load python/3.7.4

export PATH=/home/bsc08/bsc08890/.local/bin:/home/bsc08/bsc08890/bin:$PATH
export PYTHONPATH=/home/bsc08/bsc08890/.local/lib/python3.7/site-packages:$PYTHONPATH

##########################################
############# Input params ###############
##########################################

METHOD=$1
PATIENT=$2  # C141
SUFFIX=$3
PATH2TFLIST=$4  # /gpfs/projects/bsc08/bsc08890/data/TF_lists/lambert2018.txt
NUM_WORKERS=$5
RUN_TYPE=$6
Q_THRESH=$7
TASK_NUM=$8

# Adding _TF suffix if running only on TFs
if [ -z "$PATH2TFLIST" ]
then
    TF_SUFFIX=
else
    TF_SUFFIX=_TF
fi

# Adding task number label to the logs
if [ ! -z "$TASK_NUM" ]
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

PATH2DATA=/gpfs/projects/bsc08/bsc08890/res/covid_19/${PATIENT}/data/Seurat/raw_data${SUFFIX}.tsv

GRN_ADJ=/gpfs/projects/bsc08/bsc08890/res/covid_19/${PATIENT}/data/${METHOD}/raw_data${SUFFIX}${TF_SUFFIX}.tsv
GRN_ADJ_COR=/gpfs/projects/bsc08/bsc08890/res/covid_19/${PATIENT}/data/${METHOD}/raw_data${SUFFIX}${TF_SUFFIX}_cor.tsv
OUTPUT=/gpfs/projects/bsc08/bsc08890/res/covid_19/${PATIENT}/data/${METHOD}/raw_data${SUFFIX}${TF_SUFFIX}_ctx.tsv

POST_PROCESS=/gpfs/home/bsc08/bsc08890/src/adj_list_to_process.py

DB_NAMES=/gpfs/projects/bsc08/bsc08890/data/SCENIC/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather
MOTIF_ANNOTATION=/gpfs/projects/bsc08/bsc08890/data/SCENIC/motifs-v9-nr.hgnc-m0.001-o0.0.tbl

LOG_OUT=${LOG_OUTPUT}/${PATIENT}${TASK_NUM}${TF_SUFFIX}${SUFFIX}_pySCENIC_${METHOD}_${SLURM_JOBID}.out
LOG_ERR=${LOG_OUTPUT}/${PATIENT}${TASK_NUM}${TF_SUFFIX}${SUFFIX}_pySCENIC_${METHOD}_${SLURM_JOBID}.err

##########################################
##########################################
##########################################

# Creating output folder
mkdir -p /gpfs/projects/bsc08/bsc08890/res/covid_19/${PATIENT}/data/${METHOD}

# Creating log files
printf "Processing ${PATIENT}${SUFFIX} data with stored in $PATH2DATA..\n\n\n" > $LOG_OUT
printf "" > $LOG_ERR

# Running either on all genes, or only on TFs
# If running on all genes, using custom script
# If running on TFs, using pyscenic
if [ -z "$PATH2TFLIST" ]; then
   
    printf "No TFs are provided - running $METHOD based on all genes..\n\n" >> $LOG_OUT
    printf "Computing adjacencies..\n" >> $LOG_OUT

    # Run pySCENIC on all genes
    python /gpfs/home/bsc08/bsc08890/src/run_pyscenic_grn_all_genes.py -m $METHOD -i $PATH2DATA -o $GRN_ADJ --num_workers $NUM_WORKERS 2> $LOG_ERR 1> $LOG_OUT && \
    printf "Finished computing adjacencies..\n\n" >> $LOG_OUT || \
    printf "Failed computing adjacencies..\n\n" >> $LOG_OUT
    
    printf "Computing Spearman correlations..\n" >> $LOG_OUT

    # Computing the Spearman correlation for all adjacencies
    pyscenic add_cor $GRN_ADJ $PATH2DATA --output $GRN_ADJ_COR --transpose 2> $LOG_ERR 1> $LOG_OUT && \
    python $POST_PROCESS -f $GRN_ADJ_COR -q $Q_THRESH 2>> $LOG_ERR 1>> $LOG_OUT && \
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
    python $POST_PROCESS -f $GRN_ADJ_COR -q $Q_THRESH 2>> $LOG_ERR 1>> $LOG_OUT && \
    printf "Finished computing Spearman correlations..\n\n" >> $LOG_OUT || \
    printf "Failed computing Spearman correlations..\n\n" >> $LOG_OUT
    
    printf "Computing motif-enriched regulons..\n" >> $LOG_OUT

    # Filtering the list of adjacencies by picking only motif-enriched ones
    pyscenic ctx $GRN_ADJ $DB_NAMES --annotations_fname $MOTIF_ANNOTATION --expression_mtx_fname $PATH2DATA --output $OUTPUT --transpose --num_workers $NUM_WORKERS 2> $LOG_ERR 1> $LOG_OUT && \
    python $POST_PROCESS -f $OUTPUT -q $Q_THRESH 2>> $LOG_ERR 1>> $LOG_OUT && \
    printf "Finished computing motif-enriched regulons.." >> $LOG_OUT || \
    printf "Failed computing motif-enriched regulons.." >> $LOG_OUT
fi
