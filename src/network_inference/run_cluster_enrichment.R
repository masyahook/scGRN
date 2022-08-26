rm(list=ls())

suppressPackageStartupMessages({
  # if (!require("pacman")) install.packages("pacman")
  # list.of.packages <- c("BiocManager","dplyr","SingleR","Matrix","Seurat","future","pheatmap","ggplot2","optparse","hdf5r")
  # BiocManager::install("SingleR")
  # BiocManager::install('limma')
  # BiocManager::install('SingleCellExperiment')
  # pacman::p_load(list.of.packages, character.only = TRUE)
  library(optparse)
  library(reticulate)
})

pd <- import('pandas')

setwd('/gpfs/home/bsc08/bsc08890')
source('src/func.R')

get_gene_name <- function (x){
  out <- substr(x, 1, which(strsplit(x, "")[[1]]=="(") - 2)
  return(out)
}

get_score <- function(x){
  out <-(substr(x, (which(strsplit(x, "")[[1]]=="=") + 1), (nchar(x) - 1)))
  return(as.double(out))
}

get_top_n_genes <- function(x, max_n){
  out <- str_split(x, '; ')[[1]]
  gene_n <- min(c(max_n, length(out)))
  return(out[1:gene_n])
}

###################### INPUT

option_list = list(
  make_option(c("-d", "--data"), type="character", default='T_cells', help="The data to run enrichment analysis on", metavar="character"),
  make_option(c("-t", "--analysis_type"), type="character", default='ora', help="The type of analysis: either ora or gsea", metavar="character"),
  make_option(c("-a", "--algo"), type="character", default='leiden', help="The clustering algorithm", metavar="character"),
  make_option(c("-s", "--data_subset"), type="character", default='all', help="The data type", metavar="character")
)
opt_parser <- OptionParser(option_list=option_list, add_help_option = T)
opt <- parse_args(opt_parser)

if (opt$data_subset == 'all'){
  fn <- paste0(
    '/gpfs/projects/bsc08/bsc08890/res/covid_19/cell_types/', 
    sprintf('%s/data/grnboost2/%s_communities/raw_data_communities_info.pickle', opt$data, opt$algo)
  )
  
} else {
  fn <- paste0(
    '/gpfs/projects/bsc08/bsc08890/res/covid_19/cell_types/', 
    sprintf('%s/data/grnboost2/%s_communities/raw_data_%s_type_communities_info.pickle', opt$data, opt$algo, opt$data_subset)
  )
}

max_n <- 50

# T cells - all
# top_n_dotplot <- 3
# top_n_cnetplot <- 5
# xtick_dotplot_size <- 7
# ytick_dotplot_size <- 7

# Macrophages - all
top_n_dotplot <- 3
top_n_cnetplot <- 5
xtick_dotplot_size <- 9
ytick_dotplot_size <- 9

# NK cells - all
# top_n_dotplot <- 3
# top_n_cnetplot <- 5
# xtick_dotplot_size <- 7
# ytick_dotplot_size <- 6

# Monocytes - all
# top_n_dotplot <- 3
# top_n_cnetplot <- 5
# xtick_dotplot_size <- 6
# ytick_dotplot_size <- 5

# T cells - C
# top_n_dotplot <- 3
# top_n_cnetplot <- 5
# xtick_dotplot_size <- 7
# ytick_dotplot_size <- 7

# Macrophage - C
# top_n_dotplot <- 3
# top_n_cnetplot <- 5
# xtick_dotplot_size <- 7
# ytick_dotplot_size <- 7

clustering_res <- pd$read_pickle(fn)
gene_scores <- lapply(clustering_res$all_sorted_genes, function(x) get_top_n_genes(x, max_n))
genes <- lapply(gene_scores, function (x) unlist(lapply(x, get_gene_name)))
scores <- lapply(gene_scores, function (x) unlist(lapply(x, get_score)))

markers_df <- data.frame()

for (i in 1:length(genes)){
  curr_df <- data.frame(cluster=sprintf('cluster_%s', i), 
                        gene=genes[[i]],
                        centrality=scores[[i]])
  markers_df <- rbind(markers_df, curr_df)
}
colnames(markers_df) <- c('cluster', 'gene', 'centrality')

if (opt$analysis_type == 'ora'){
  out <- run_ora(markers_df, is_clusters = T, 
                 cell_type_for_community_ana=paste(opt$data, opt$data_subset),
                 top_n_dotplot=top_n_dotplot, top_n_cnetplot=top_n_cnetplot,
                 xtick_dotplot_size=xtick_dotplot_size, 
                 ytick_dotplot_size=ytick_dotplot_size)
} else if (opt$analysis_type == 'gsea') {
  out <- run_gsea(markers_df, is_clusters = T)
}
