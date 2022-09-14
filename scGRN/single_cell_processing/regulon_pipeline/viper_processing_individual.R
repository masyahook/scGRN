rm(list=ls())

suppressPackageStartupMessages({
  #if (!require("pacman")) install.packages("pacman")
  #list.of.packages <- c("BiocManager","dplyr","SingleR","Matrix","Seurat","future","pheatmap","ggplot2","optparse","hdf5r")
  #BiocManager::install("SingleR")
  #BiocManager::install('limma')
  #BiocManager::install('SingleCellExperiment')
  #pacman::p_load(list.of.packages, character.only = TRUE)
  library(parallel)
  library(SingleR)
  library(Matrix)
  library(Seurat)
  library(dplyr)
  library(future)
  library(pheatmap)
  library(ggplot2)
  library(optparse)
  library(dorothea)
  library(tibble)
  library(tidyr)
})

###################### INPUT

option_list = list(
  make_option(c("-p", "--pat"), type="character", default='C51', help="Patient ID", metavar="character"),
  make_option(c("-c", "--pleiotropy_correction"), type="logical", default=F, help='Pleitropy correction?', metavar="T|F")
)
opt_parser <- OptionParser(option_list=option_list, add_help_option = T)
opt <- parse_args(opt_parser)

###################### LOAD DATA

setwd('~/')

# read metadata
meta.path <- '/gpfs/projects/bsc08/bsc08890/data/GSE145926_RAW/metadata.tsv'
meta <- read.table(meta.path, header=T, sep="\t", stringsAsFactors = F)

opt$outdir <- paste0("/gpfs/projects/bsc08/bsc08890/res/covid_19/", opt$pat)
opt$outdat <- paste0(opt$outdir, "/data/Seurat")
opt$outfig <- paste0(opt$outdir, "/figs/Seurat")
opt$patpath <- meta[which(meta$id == opt$pat),]$file
opt$verbose <- T

# Delete all previously generated files
file.remove(list.files(opt$outdat, full.names = T, recursive = T, pattern='.*dorothea.*'))

cat("\n\n")
cat("***********************************\n")
cat("*** VIPER PROCESSING BASED ON DOROTHEA REGULONS ***\n")
cat("***********************************\n\n")
cat("Outdir: ", opt$outdir, "\n")
cat("Patient: ", opt$pat, "\n")
cat("Verbose: ", opt$verbose, "\n\n")

my_ggsave <- function(obj, filename){
  suppressMessages(ggsave(obj, filename = filename))
}

###################### INIT
cat("      - Init\n")

load(paste0(opt$outdat, "/seurat_object.RData"))

###################### VIPER BASED ON DOROTHEA REGULONS

# Getting dorothea regulons, subsetting only confident ones
dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))
regulon <- dorothea_regulon_human %>%
  dplyr::filter(confidence %in% c("A","B","C"))

# Running viper
cat("      - Running viper\n")
sobj <- run_viper(sobj, regulon,
                  options = list(method = "scale", minsize = 4, 
                                 pleiotropy = opt$pleiotropy_correction, 
                                 eset.filter = FALSE, 
                                 cores = detectCores(), verbose = FALSE))

# Running viper-based pipeline
DefaultAssay(object = sobj) <- "dorothea"

###################### SCALE DATA
cat("      - Scaling\n")
sobj <- ScaleData(object = sobj, verbose = F)

###################### DIM REDUCTION
cat("      - Dimensionality reduction\n")

## PCA
sobj <- RunPCA(object = sobj, features = rownames(sobj), verbose = F)
# loadings
sobj@misc$plots$dorothea_pca_loadings <- VizDimLoadings(sobj, dims = 1:4, reduction = "pca")
my_ggsave(sobj@misc$plots$dorothea_pca_loadings, filename = paste0(opt$outfig, "/dorothea_pca_loadings.png"))
sobj@misc$plots$dorothea_pca_loadings
# elbow
sobj@misc$plots$dorothea_pca_elbow <- ElbowPlot(sobj, ndims = 40)
my_ggsave(sobj@misc$plots$dorothea_pca_elbow, filename = paste0(opt$outfig, "/dorothea_pca_elbow.png"))
sobj@misc$plots$dorothea_pca_elbow
#
sobj@misc$plots$dorothea_pca_main <- DimPlot(object = sobj, reduction = "pca")
my_ggsave(sobj@misc$plots$dorothea_pca_main, filename = paste0(opt$outfig, "/pca_main.png"))
sobj@misc$plots$dorothea_pca_main
sobj@misc$plots$dorothea_pca_heatmap <- DimHeatmap(sobj, dims = 1:4, ncol=2, nfeatures = 30, cells = 1000, balanced = TRUE, fast = F)
my_ggsave(sobj@misc$plots$dorothea_pca_heatmap, filename = paste0(opt$outfig, "/dorothea_pca_heatmap.png"))
sobj@misc$plots$dorothea_pca_heatmap

## T-SNE
sobj <- RunTSNE(object = sobj, dims = 1:10, verbose = F)
sobj@misc$plots$dorothea_tsne_main <- DimPlot(object = sobj, reduction = "tsne")
my_ggsave(sobj@misc$plots$dorothea_tsne_main, filename = paste0(opt$outfig, "/dorothea_tsne_main.png"))
sobj@misc$plots$dorothea_tsne_main

## U-MAP
sobj <- RunUMAP(object = sobj, dims=1:10, verbose = F, umap.method = "uwot", metric = "cosine")
sobj@misc$plots$dorothea_umap_main <- DimPlot(object = sobj, reduction = "umap")
my_ggsave(sobj@misc$plots$dorothea_umap_main, filename = paste0(opt$outfig, "/dorothea_umap_main.png"))
sobj@misc$plots$dorothea_umap_main

###################### FIND CLUSTERS
cat("      - Find clusters\n")

sobj <- FindNeighbors(sobj, dims = 1:10, verbose = FALSE)
sobj <- FindClusters(sobj, resolution = 0.5, verbose = FALSE)

sobj@misc$plots$dorothea_pca_clusters <- DimPlot(object = sobj, reduction = "pca")
my_ggsave(sobj@misc$plots$dorothea_pca_clusters, filename = paste0(opt$outfig, "/dorothea_pca_clusters.png"))
sobj@misc$plots$dorothea_pca_clusters
sobj@misc$plots$dorothea_tsne_clusters <- DimPlot(object = sobj, reduction = "tsne")
my_ggsave(sobj@misc$plots$dorothea_tsne_clusters, filename = paste0(opt$outfig, "/dorothea_tsne_clusters.png"))
sobj@misc$plots$dorothea_tsne_clusters
sobj@misc$plots$dorothea_umap_clusters <- DimPlot(object = sobj, reduction = "umap")
my_ggsave(sobj@misc$plots$dorothea_umap_clusters, filename = paste0(opt$outfig, "/dorothea_umap_clusters.png"))
sobj@misc$plots$dorothea_umap_clusters

## We transform Viper scores, scaled by seurat, into a data frame to better 
## handling the results
viper_scores_df <- GetAssayData(sobj, 
                                slot = "scale.data", 
                                assay = "dorothea") %>%
  data.frame(check.names = F) %>%
  t()

## We create a data frame containing the cells and their clusters
CellsClusters <- data.frame(cell = names(Idents(sobj)), 
                            cell_type = sobj@misc$cell_type_pred$labels,
                            check.names = F)

## We create a data frame with the Viper score per cell and its clusters
viper_scores_clusters <- viper_scores_df  %>%
  data.frame() %>% 
  rownames_to_column("cell") %>%
  gather(tf, activity, -cell) %>%
  inner_join(CellsClusters)

## We summarize the Viper scores by cell population
summarized_viper_scores <- viper_scores_clusters %>% 
  group_by(tf, cell_type) %>%
  summarise(avg = mean(activity),
            std = sd(activity))

## We select the 20 most variable TFs. (20 * num_clusters)
num_clusters <- length(unique((Idents(sobj))))
highly_variable_tfs <- summarized_viper_scores %>%
  group_by(tf) %>%
  mutate(var = var(avg))  %>%
  ungroup() %>%
  top_n(20*num_clusters, var) %>%
  distinct(tf)

## We prepare the data for the plot
summarized_viper_scores_df <- summarized_viper_scores %>%
  semi_join(highly_variable_tfs, by = "tf") %>%
  dplyr::select(-std) %>%   
  spread(tf, avg) %>%
  data.frame(row.names = 1, check.names = FALSE) 

palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(min(summarized_viper_scores_df), 0, 
                   length.out=ceiling(palette_length/2) + 1),
               seq(max(summarized_viper_scores_df)/palette_length, 
                   max(summarized_viper_scores_df), 
                   length.out=floor(palette_length/2)))

sobj@misc$plots$dorothea_tf_activity <- pheatmap(t(summarized_viper_scores_df),fontsize=14, 
                       fontsize_row = 10, 
                       color=my_color, breaks = my_breaks, 
                       main = "DoRothEA TF activity", angle_col = 45,
                       treeheight_col = 0,  border_color = NA) 
my_ggsave(sobj@misc$plots$dorothea_tf_activity, filename = paste0(opt$outfig, "/dorothea_tf_activity.png"))
sobj@misc$plots$dorothea_tf_activity

###################### FINISH
cat("      - Saving\n")

save(sobj, file=paste0(opt$outdat, "/dorothea_seurat_object.RData"))

write.table(sobj@assays$dorothea@data, file=paste0(opt$outdat, "/dorothea_data.tsv"), sep="\t", row.names=T, col.names=colnames(sobj@assays$dorothea@data), quote=F)
write.table(sobj@assays$dorothea@scale.data, file=paste0(opt$outdat, "/dorothea_scaled_data.tsv"), sep="\t", row.names=T, col.names=colnames(sobj@assays$dorothea@scale.data), quote=F)

#   cell-specific
cell_types <- unique(sobj[['cell_type']][['cell_type']])
for (i in 1:length(cell_types))
{
  type <- cell_types[i]
  cell_type_sobj <- subset(x=sobj, subset = cell_type == type)
  write.table(cell_type_sobj@assays$dorothea@data, file=paste0(opt$outdat, sprintf("/dorothea_data_%s.tsv", type)), sep="\t", row.names=T, col.names=colnames(cell_type_sobj@assays$dorothea@data), quote=F)
}

cat("\n\n[finished]\n\n")
