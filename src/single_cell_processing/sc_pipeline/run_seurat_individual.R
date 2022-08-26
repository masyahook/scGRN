rm(list=ls())  # clear namespace
suppressPackageStartupMessages({
  library(SingleR)
  library(Matrix)
  library(Seurat)
  library(dplyr)
  library(pheatmap)
  library(ggplot2)
  library(optparse)
  library(tidyr)
})  # load all packages

# save figures without messages
my_ggsave <- function(obj, filename){
  suppressMessages(ggsave(obj, filename = filename))
}

# deal with duplicate slashes
file_path = function(..., fsep = .Platform$file.sep){
  gsub("//", "/", file.path(..., fsep = fsep))
}

# prevent Rplots.pdf from being generated
pdf(NULL)

###################### INPUT

option_list = list(
  make_option(c("-p", "--pat"), type="character", help="Patient ID", metavar="character"),
  make_option(c("-m", "--meta_file"), type="character", help="Metadata file", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", help="Output folder", metavar="character"),
  make_option(c('-a', '--annotation'), type='character', help='Annotation dataset used to label cells', default='HumanPrimaryCellAtlasData', metavar='character'),
  make_option(c('--annotation_folder'), type='character', help='Folder with saved annotation datasets', default=NULL, metavar='character'),
  make_option(c("-s", "--serialize"), type="logical", default=T, help="Save Seurat object", metavar="T|F"),
  make_option(c('--add_figure_objects'), type='logical', default=T, help='Whether to add figure objects to Seurat object', metavar='T|F'),
  make_option(c("-v", "--verbose"), type="logical", default=F, help="Verbose level", metavar="T|F")
)
opt_parser <- OptionParser(option_list=option_list, add_help_option = T)
opt <- parse_args(opt_parser)

# mandatory params
if (is.null(opt$pat) || opt$pat == ''){
  print_help(opt_parser)
  stop("No patient ID provided", call.=FALSE)
}
if (is.null(opt$meta_file) || opt$meta_file == ''){
  print_help(opt_parser)
  stop("No metadata file path provided", call.=FALSE)
}
if (is.null(opt$outdir) || opt$outdir == ''){
  print_help(opt_parser)
  stop("No output folder path provided", call.=FALSE)
}

# arbitrary params
if (opt$annotation == ''){
  opt$annotation = 'HumanPrimaryCellAtlasData'
}
if (is.na(opt$serialize)){
  opt$serialize = T
}
if (is.na(opt$add_figure_objects)){
  opt$add_figure_objects = T
}
if (is.na(opt$verbose)){
  opt$verbose = F
}

###################### LOAD DATA

# read metadata
meta <- read.table(opt$meta_file, header=T, sep="\t", stringsAsFactors = F)

# set up working directories
opt$outdir <- file_path(opt$outdir, opt$pat)
opt$outdat <- file_path(opt$outdir, "data/Seurat")
opt$outfig <- file_path(opt$outdir, "figs/Seurat")
opt$patpath <- meta[which(meta$id == opt$pat),]$file

# create working directories
dir.create(opt$outdir, recursive=T)
dir.create(opt$outdat, recursive=T)
dir.create(opt$outfig, recursive=T)

# delete previously generated files
invisible(file.remove(list.files(opt$outdat, full.names = T, recursive = T)))
invisible(file.remove(list.files(opt$outfig, full.names = T, recursive = T)))

cat("\n\n")
cat("**********************************************************\n")
cat("*** SINGLE-CELL DATA PROCESSING FOR INDIVIDUAL PATIENT ***\n")
cat("**********************************************************\n\n")
cat("Patient: ", opt$pat, "\n")
cat("Outdir: ", opt$outdir, "\n")
cat("Metadata: ", opt$meta_file, "\n")
cat("Annotation: ", opt$annotation, "\n")
cat("Annotation folder: ", opt$annotation_folder, "\n")
cat("Serialize output: ", opt$serialize, "\n")
cat("Save figure objects: ", opt$save_figure_objects, '\n')
cat("Verbose: ", opt$verbose, "\n")

###################### LOAD DATA
cat("      - Init\n")

# init Seurat object
raw_counts <- Read10X_h5(opt$patpath, use.names = T, unique.features = T)
sobj <- CreateSeuratObject(counts = raw_counts, min.cells = 3, min.features = 200, project = opt$pat)

sobj@misc$plots <- list()

###################### QC
cat("      - QC\n")

# calculate percentage of MT gene counts
sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")

# plots for QC features
sobj@misc$plots$qc_violin <- VlnPlot(sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = .25)
my_ggsave(sobj@misc$plots$qc_violin, filename = file_path(opt$outfig, "qc_violin.png"))
sobj@misc$plots$qc_violin

plot1 <- FeatureScatter(sobj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
sobj@misc$plots$qc_cross <- plot1 + plot2
my_ggsave(sobj@misc$plots$qc_cross, filename = file_path(opt$outfig, "qc_cross.png"))
sobj@misc$plots$qc_cross

###################### FILTERING
cat("      - Filtering\n")

# empty droplets or double droplets ?
sobj@misc$pre_filt_dim <- dim(sobj@assays$RNA)
sobj <- subset(sobj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
sobj@misc$post_filt_dim <- dim(sobj@assays$RNA)

###################### NORMALIZATION
cat("      - Normalization\n")

sobj <- NormalizeData(object = sobj, verbose = opt$verbose)

###################### VARIABLE FEATURES
cat("      - Find variable features\n")

# getting variable genes, plotting top-10 ones
sobj <- FindVariableFeatures(object = sobj, selection.method = "vst", nfeatures = 2000, verbose = opt$verbose)
top10 <- head(VariableFeatures(sobj), 10)
plot1 <- VariableFeaturePlot(sobj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
sobj@misc$plots$variable_features <- plot1 + plot2
sobj@misc$plots$variable_features
my_ggsave(sobj@misc$plots$variable_features, filename = file_path(opt$outfig, "variable_features.png"))

###################### SCALE DATA
cat("      - Scaling\n")

sobj <- ScaleData(object = sobj, verbose = opt$verbose)

###################### DIM REDUCTION
cat("      - Dimensionality reduction\n")

## PCA
sobj <- RunPCA(object = sobj, features = VariableFeatures(object = sobj), verbose = opt$verbose)
# loadings
sobj@misc$plots$pca_loadings <- VizDimLoadings(sobj, dims = 1:4, reduction = "pca")
my_ggsave(sobj@misc$plots$pca_loadings, filename = file_path(opt$outfig, "pca_loadings.png"))
sobj@misc$plots$pca_loadings
# elbow
sobj@misc$plots$pca_elbow <- ElbowPlot(sobj, ndims = 40)
my_ggsave(sobj@misc$plots$pca_elbow, filename = file_path(opt$outfig, "pca_elbow.png"))
sobj@misc$plots$pca_elbow
# PCA scatter plot and PC heatmap plot
sobj@misc$plots$pca_main <- DimPlot(object = sobj, reduction = "pca")
my_ggsave(sobj@misc$plots$pca_main, filename = file_path(opt$outfig, "pca_main.png"))
sobj@misc$plots$pca_main
sobj@misc$plots$pca_heatmap <- DimHeatmap(sobj, dims = 1:4, ncol=2, nfeatures = 30, cells = 1000, balanced = TRUE, fast = F)
my_ggsave(sobj@misc$plots$pca_heatmap, filename = file_path(opt$outfig, "pca_heatmap.png"))
sobj@misc$plots$pca_heatmap

# # number of components
# sobj <- JackStraw(sobj, num.replicate = 100, dims=20)
# sobj <- ScoreJackStraw(sobj, dims = 1:20)
# JackStrawPlot(sobj, dims = 1:20)

## T-SNE
sobj <- RunTSNE(object = sobj, dims = 1:20, verbose = opt$verbose)
sobj@misc$plots$tsne_main <- DimPlot(object = sobj, reduction = "tsne")
my_ggsave(sobj@misc$plots$tsne_main, filename = file_path(opt$outfig, "tsne_main.png"))
sobj@misc$plots$tsne_main

## U-MAP
sobj <- RunUMAP(object = sobj, dims=1:20, verbose = opt$verbose)
sobj@misc$plots$umap_main <- DimPlot(object = sobj, reduction = "umap")
my_ggsave(sobj@misc$plots$umap_main, filename = file_path(opt$outfig, "umap_main.png"))
sobj@misc$plots$umap_main

###################### FIND CLUSTERS
cat("      - Find clusters\n")

sobj <- FindNeighbors(object = sobj, verbose = opt$verbose)
sobj <- FindClusters(object = sobj, verbose = opt$verbose)

# plotting clusters using PCA, T-SNE, UMAP
sobj@misc$plots$pca_clusters <- DimPlot(object = sobj, reduction = "pca")
my_ggsave(sobj@misc$plots$pca_clusters, filename = file_path(opt$outfig, "pca_clusters.png"))
sobj@misc$plots$pca_clusters
sobj@misc$plots$tsne_clusters <- DimPlot(object = sobj, reduction = "tsne")
my_ggsave(sobj@misc$plots$tsne_clusters, filename = file_path(opt$outfig, "tsne_clusters.png"))
sobj@misc$plots$tsne_clusters
sobj@misc$plots$umap_clusters <- DimPlot(object = sobj, reduction = "umap")
my_ggsave(sobj@misc$plots$umap_clusters, filename = file_path(opt$outfig, "umap_clusters.png"))
sobj@misc$plots$umap_clusters

###################### GENE MARKERS
cat("      - Gene markers by cluster\n")

# find markers for every cluster compared to all remaining cells, report only the positive ones
sobj.markers <- FindAllMarkers(sobj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = opt$verbose)
# top-2 markers per cluster, print them
sobj.markers %>% group_by(cluster) %>% top_n(n = 2, wt = "avg_log2FC")
sobj@misc$markers <- sobj.markers
# plot top-10 markers per cluster
top10 <- sobj.markers %>% group_by(cluster) %>% top_n(n = 10, wt = "avg_log2FC")
sobj@misc$plots$markers <- DoHeatmap(sobj, features = top10$gene) + NoLegend()
my_ggsave(sobj@misc$plots$markers, filename = file_path(opt$outfig, "markers.png"))
sobj@misc$plots$markers

###################### CELL TYPE IDENTIFICATION
cat("      - Cell type identification\n")

# Bioconductor obj
sobj_se <- as.SingleCellExperiment(sobj)

# loading reference data
if (!is.null(opt$annotation_folder)){
  load(file=file_path(opt$annotation_folder, sprintf('%s.RData', opt$annotation)))
} else {
  # download reference dataset for annotation
  if (opt$annotation == 'HumanPrimaryCellAtlasData') {
    ref <- SingleR::HumanPrimaryCellAtlasData()
  } else if (opt$annotation == 'BlueprintEncodeData') {
    ref <- SingleR::DatabaseImmuneCellExpressionData()
  } else if (opt$annotation == 'DatabaseImmuneCellExpressionData') {
    ref <- SingleR::BlueprintEncodeData()
  } else {
    ref <- SingleR::MonacoImmuneData()
  }
}

# prediction and saving to object
pred <- SingleR(test=sobj_se, ref=ref, labels=ref$label.main)
sobj@misc$cell_type_ref <- ref
sobj@misc$cell_type_pred <- pred
sobj@misc$plots$cell_type_pred_dist <- plotScoreDistribution(pred, size = .25)
my_ggsave(sobj@misc$plots$cell_type_pred_dist, filename = file_path(opt$outfig, "cell_type_pred_dist.png"))
sobj@misc$plots$cell_type_pred_heatmap <- plotScoreHeatmap(pred, show.pruned = T)
my_ggsave(sobj@misc$plots$cell_type_pred_heatmap, filename = file_path(opt$outfig, "cell_type_pred_heatmap.png"))

# cell type freqs
png(file_path(opt$outfig, "cell_type_freqs.png"), width=800, height=700, pointsize = 20)
par(mar=c(5,10,3,3))
freqs <- sort(table(pred$labels), decreasing = F)
barplot(freqs, horiz=T, las=2)

# cell type to clusters
pheatmap(log(10+table(pred$labels, sobj$seurat_clusters)), filename = file_path(opt$outfig, "cell_type_to_cluster.png"))

# dim reduction plots
sobj[["cell_type"]] <- pred$labels
sobj@misc$plots$pca_cell_types <- DimPlot(object = sobj, reduction = "pca", group.by = "cell_type")
my_ggsave(sobj@misc$plots$pca_cell_types, filename = file_path(opt$outfig, "pca_cell_types.png"))
sobj@misc$plots$tsne_cell_types <- DimPlot(object = sobj, reduction = "tsne", group.by = "cell_type")
my_ggsave(sobj@misc$plots$tsne_cell_types, filename = file_path(opt$outfig, "tsne_cell_types.png"))
sobj@misc$plots$umap_cell_types <- DimPlot(object = sobj, reduction = "umap", group.by = "cell_type")
my_ggsave(sobj@misc$plots$umap_cell_types, filename = file_path(opt$outfig, "umap_cell_types.png"))

# summarize by cell type
summ_fun <- function(x) mean(x, trim=.005)
summ_cells <- function(x, fun=summ_fun){
  apply(x,2,fun)
}
sobj@misc$cell_type_assay <- do.call("cbind", by(t(sobj@assays$RNA@data), pred$labels, summ_cells))

###################### FINISH
cat("      - Saving\n")

if (opt$add_figure_objects == F) {
  sobj_combined@misc$plots <- NULL
}

if (opt$serialize == T) save(sobj, file=file_path(opt$outdat, "seurat_object.RData"))

################# Load Seurat object from previous run ####################
# load(file_path(opt$outdat, "seurat_object.RData"))
###########################################################################

# prune labels - keep only cells with high confidence annotation
cell_types <- unique(pred$labels)
sobj_to_save <- sobj

# data
#   general
write.table(sobj_to_save@assays$RNA@data, file=file_path(opt$outdat, "norm_data.tsv"), sep="\t", row.names=T, col.names=colnames(sobj_to_save@assays$RNA@data), quote=F)
write.table(sobj_to_save@assays$RNA@counts, file=file_path(opt$outdat, "raw_data.tsv"), sep="\t", row.names=T, col.names=colnames(sobj_to_save@assays$RNA@counts), quote=F)
write.table(sobj_to_save@assays$RNA@scale.data, file=file_path(opt$outdat, "scaled_data.tsv"), sep="\t", row.names=T, col.names=colnames(sobj_to_save@assays$RNA@scale.data), quote=F)
#   cell info
write.table(sobj_to_save@meta.data, file=file_path(opt$outdat, "cells_metadata.tsv"), sep="\t", row.names=T, col.names=T, quote=F)

#   cell type-specific
for (i in 1:length(cell_types))
{
  type <- cell_types[i]
  cell_type_sobj_to_save <- subset(x=sobj_to_save, subset = cell_type == type)
  write.table(cell_type_sobj_to_save@assays$RNA@counts, file=file_path(opt$outdat, sprintf("raw_data_%s.tsv", type)), sep="\t", row.names=T, col.names=colnames(cell_type_sobj_to_save@assays$RNA@counts), quote=F)
  # write.table(cell_type_sobj_to_save@assays$RNA@data, file=file_path(opt$outdat, sprintf("norm_data_%s.tsv", type)), sep="\t", row.names=T, col.names=colnames(cell_type_sobj_to_save@assays$RNA@data), quote=F)
}

cat("\n\n[finished]\n\n")
