rm(list=ls())  # clear namespace

suppressPackageStartupMessages({
  library(SingleR)
  library(dplyr)
  library(Matrix)
  library(Seurat)
  library(future)
  library(pheatmap)
  library(ggplot2)
  library(optparse)
  library(tibble)
  library(tidyr)
  library(parallel)
  library(stringr)
  library(SeuratData)
  library(patchwork)
  library(future)
  library(data.table)
  library(grid)
})  # load all packages

# save figures without messages
my_ggsave <- function(obj, filename){
  suppressMessages(ggsave(obj, filename = filename))
}

# deal with duplicate slashes
file_path = function(..., fsep = .Platform$file.sep){
  gsub("//", "/", file.path(..., fsep = fsep))
}

# prevent Rplots.pdf from being generated (turn off to visualize in RStudio Server)
pdf(NULL)

###################### INPUT

option_list = list(
  make_option(c("-m", "--meta_file"), type="character", help="Metadata file", metavar="character"),  
  make_option(c("-o", "--outdir"), type="character", help="Output folder", metavar="character"),
  make_option(c('-n', '--num_proc'), type='integer', help='Number of processes run in parallel', default=4, metavar='integer'),
  make_option(c('--use_merged_precomputed'), type='logical', default=F, help='Whether to use pre-computed merged Seurat object', metavar='T|F'),
  make_option(c("-s", "--serialize"), type="logical", default=T, help="Save Seurat object", metavar="T|F"),
  make_option(c("-v", "--verbose"), type="logical", default=T, help="Verbose", metavar="T|F")
)
opt_parser <- OptionParser(option_list=option_list, add_help_option = T)
opt <- parse_args(opt_parser)

# mandatory params
if (is.null(opt$meta_file) || opt$meta_file == ''){
  print_help(opt_parser)
  stop("No metadata file path provided", call.=FALSE)
}
if (is.null(opt$outdir) || opt$outdir == ''){
  print_help(opt_parser)
  stop("No output folder path provided", call.=FALSE)
}

# arbitrary params
if (is.na(opt$num_proc)){
  opt$num_proc = 4
}
if (is.na(opt$serialize)){
  opt$serialize = T
}
if (is.na(opt$verbose)){
  opt$verbose = F
}

cat("\n\n")
cat("***********************************\n")
cat("*** MERGE DATASETS BY CELL TYPE ***\n")
cat("***********************************\n\n")
cat("Outdir: ", opt$outdir, "\n")
cat("Metadata: ", opt$meta_file, "\n")
cat("Parallelized by: ", opt$num_proc, '\n')
cat("Use pre-computed merged Seurat object: ", opt$use_merged_precomputed, "\n")
cat("Serialize output: ", opt$serialize, "\n")
cat("Verbose: ", opt$verbose, "\n")

###################### LOAD DATA

# read metadata
meta <- read.table(opt$meta_file, header=T, sep="\t", stringsAsFactors = F)
rownames(meta) <- meta$id

pat_types <- meta$group
names(pat_types) <- meta$id
pat_levels <- meta$id
pat_type_levels <- c('C', 'M', 'S')

# create directory where we will save cell-type specific data
cell_type_dir <- file_path(opt$outdir, 'cell_types')
dir.create(cell_type_dir, showWarnings = F)

# create sample directory for the whole dataset
full_dir <- file_path(cell_type_dir, 'all_data')
dir.create(full_dir, showWarnings = F)
full_fig_dir <- file_path(full_dir, 'figs/Seurat')
dir.create(full_fig_dir, showWarnings = F, recursive=T)
full_data_dir <- file_path(full_dir, 'data/Seurat')
dir.create(full_data_dir, showWarnings = F, recursive=T)

# define colors
colors <- list(green='#39B600', yellow='#D89000', red='#F8766D', blue='#00B0F6', 
               purple='#9590FF', cyan='#00BFC4', pink='E76BF3', 
               light_pink='#FF62BC', saturated_green='#00BF7D')

# Setting up parallelization
options(future.globals.maxSize = 2700 * 1024^2)
plan("multiprocess", workers = opt$num_proc)

############# DATA INTEGRATION
cat("      - Integrating data\n")

if (opt$use_merged_precomputed == T){
  
  # load previously merged Seurat object
  load(file=file_path(full_data_dir, "unprocessed_seurat_object.RData"))
  
} else {
  
  # loading patient objects
  pat_objs <- list()
  for (pat_id in names(pat_types)){
    
    # define paths to files
    curr_fn <- file_path(opt$outdir, pat_id, 'data/Seurat/raw_data.tsv')
    curr_meta_fn <- file_path(dirname(curr_fn), 'cells_metadata.tsv')
    
    # read the data and metadata
    curr_counts <- read.table(file=curr_fn, sep="\t", check.names=F)
    curr_meta <- read.table(file=curr_meta_fn, sep='\t', check.names=F)
    
    # rename orig.dent in metadata to patient ID
    curr_meta$pat_id <- pat_id
    
    # adding patient type data to metadata
    curr_meta$pat_type <- pat_types[pat_id]
    
    # load everything to Seurat object
    curr_obj <- CreateSeuratObject(counts = curr_counts, min.cells = 0, min.genes = 0, project = pat_id, meta.data = curr_meta)
    
    # save into the list
    pat_objs[pat_id] <- curr_obj
    
  }
  
  # normalize and identify variable features for each dataset independently
  pat_objs <- lapply(X = pat_objs, FUN = function(x) {
    x <- NormalizeData(x, verbose = opt$verbose)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = opt$verbose)
  })
  
  # select features that are repeatedly variable across datasets for integration
  features <- SelectIntegrationFeatures(object.list = pat_objs, verbose = opt$verbose)
  anchors <- FindIntegrationAnchors(object.list = pat_objs, anchor.features = features, verbose = opt$verbose)
  sobj_combined <- IntegrateData(anchorset = anchors, verbose = opt$verbose)
  
  # saving data
  save(sobj_combined, file=file_path(full_data_dir, "unprocessed_seurat_object.RData"))
  
  # clearing namespace
  rm(list=c("pat_objs", 'anchors', 'features', 'curr_obj'))
}

# refactor the patients
sobj_combined$pat_id <- factor(x = sobj_combined$pat_id, levels=pat_levels)
sobj_combined$pat_type <- factor(x = sobj_combined$pat_type, level=pat_type_levels)

# run typical analysis on the integrated assay
DefaultAssay(sobj_combined) <- "integrated"

###################### SCALE DATA
cat("      - Scaling\n")

sobj_combined <- ScaleData(object = sobj_combined, verbose = opt$verbose)

###################### DIM REDUCTION
cat("      - Dimensionality reduction\n")

# PCA
sobj_combined <- RunPCA(object = sobj_combined, features = VariableFeatures(object = sobj_combined), verbose = opt$verbose)
# loadings
sobj_combined@misc$plots$pca_loadings <- VizDimLoadings(sobj_combined, dims = 1:4, reduction = "pca")
my_ggsave(sobj_combined@misc$plots$pca_loadings, filename = file_path(full_fig_dir, "pca_loadings.png"))
# elbow
sobj_combined@misc$plots$pca_elbow <- ElbowPlot(sobj_combined, ndims = 40)
my_ggsave(sobj_combined@misc$plots$pca_elbow, filename = file_path(full_fig_dir, "pca_elbow.png"))
# PCA scatter plot and PC heatmap plot
sobj_combined@misc$plots$pca_main <- DimPlot(object = sobj_combined, reduction = "pca")
my_ggsave(sobj_combined@misc$plots$pca_main, filename = file_path(full_fig_dir, "pca_main.png"))
sobj_combined@misc$plots$pca_heatmap <- DimHeatmap(sobj_combined, dims = 1:4, ncol=2, nfeatures = 30, cells = 1000, balanced = TRUE, fast = F)
my_ggsave(sobj_combined@misc$plots$pca_heatmap, filename = file_path(full_fig_dir, "pca_heatmap.png"))

# # number of components
# sobj_combined <- JackStraw(sobj_combined, num.replicate = 100, dims=20)
# sobj_combined <- ScoreJackStraw(sobj_combined, dims = 1:20)
# JackStrawPlot(sobj_combined, dims = 1:20)

## T-SNE
sobj_combined <- RunTSNE(object = sobj_combined, dims = 1:20, verbose = opt$verbose)
sobj_combined@misc$plots$tsne_main <- DimPlot(object = sobj_combined, reduction = "tsne")
my_ggsave(sobj_combined@misc$plots$tsne_main, filename = file_path(full_fig_dir, "tsne_main.png"))

## U-MAP
sobj_combined <- RunUMAP(object = sobj_combined, dims=1:20, verbose = opt$verbose)
sobj_combined@misc$plots$umap_main <- DimPlot(object = sobj_combined, reduction = "umap")
my_ggsave(sobj_combined@misc$plots$umap_main, filename = file_path(full_fig_dir, "umap_main.png"))


###################### FIND CLUSTERS
cat("      - Find clusters\n")

sobj_combined <- FindNeighbors(object = sobj_combined, verbose = opt$verbose)
sobj_combined <- FindClusters(object = sobj_combined, verbose = opt$verbose)

# grouping by clusters
sobj_combined@misc$plots$pca_clusters <- DimPlot(object = sobj_combined, reduction = "pca")
my_ggsave(sobj_combined@misc$plots$pca_clusters, filename = file_path(full_fig_dir, "pca_clusters.png"))
sobj_combined@misc$plots$tsne_clusters <- DimPlot(object = sobj_combined, reduction = "tsne")
my_ggsave(sobj_combined@misc$plots$tsne_clusters, filename = file_path(full_fig_dir, "tsne_clusters.png"))
sobj_combined@misc$plots$umap_clusters <- DimPlot(object = sobj_combined, reduction = "umap")
my_ggsave(sobj_combined@misc$plots$umap_clusters, filename = file_path(full_fig_dir, "umap_clusters.png"))

# grouping by patient and patient type
p1 <- DimPlot(sobj_combined, reduction = "pca", group.by = "pat_type")
p2 <- DimPlot(sobj_combined, reduction = "pca", group.by = 'pat_id', label = TRUE, repel = TRUE)
sobj_combined@misc$plots$pca_pat <- p1 + p2
my_ggsave(sobj_combined@misc$plots$pca_pat, filename = file_path(full_fig_dir, "pca_pat.png"))

p2 <- DimPlot(sobj_combined, reduction = "tsne", group.by = "pat_type")
p2 <- DimPlot(sobj_combined, reduction = "tsne", group.by = 'pat_id', label = TRUE, repel = TRUE)
sobj_combined@misc$plots$tsne_pat <- p1 + p2
my_ggsave(sobj_combined@misc$plots$tsne_pat, filename = file_path(full_fig_dir, "tsne_pat.png"))

p1 <- DimPlot(sobj_combined, reduction = "umap", group.by = "pat_type")
p2 <- DimPlot(sobj_combined, reduction = "umap", group.by = 'pat_id', label = TRUE, repel = TRUE)
sobj_combined@misc$plots$umap_pat <- p1 + p2
my_ggsave(sobj_combined@misc$plots$umap_pat, filename = file_path(full_fig_dir, "umap_pat.png"))

sobj_combined@misc$plots$umap_pat_split <- 
    DimPlot(sobj_combined, reduction = "umap", split.by = "pat_id", group.by = 'pat_type', cols=c(colors$green, colors$yellow, colors$red)) + 
    theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7)) + labs(title="Patient ID")
my_ggsave(sobj_combined@misc$plots$umap_pat_split, filename = file_path(full_fig_dir, "umap_pat_split.png"))

###################### GENE MARKERS
cat("      - Gene markers by patient and patient type\n")

# next steps are computationally heave, will execute sequentially
plan("sequential")

# changing assay to unmodified data
DefaultAssay(sobj_combined) <- "RNA"
sobj_combined <- ScaleData(object = sobj_combined, verbose = opt$verbose)

# changing the identities to pat_id for DEA between pat ids
cluster_idents <- Idents(object = sobj_combined)
Idents(object = sobj_combined) <- sobj_combined@meta.data$pat_id

# find markers for every patient compared to all remaining cells, report only the positive ones
sobj.markers <- FindAllMarkers(sobj_combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = opt$verbose)
# top-2 markers per cluster, print them
sobj.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
sobj_combined@misc$markers_pat_id <- sobj.markers
# plot top-10 markers per cluster
top10 <- sobj.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
sobj_combined@misc$plots$markers_pat_id <- DoHeatmap(subset(sobj_combined, downsample = 1000), features = top10$gene) + NoLegend()
my_ggsave(sobj_combined@misc$plots$markers_pat_id, filename = file_path(full_fig_dir, "markers_pat_id.png"))

# changing the identities to pat type for DEA between pat types
Idents(object = sobj_combined) <- sobj_combined@meta.data$pat_type

# find markers for every patient group compared to all remaining cells, report only the positive ones
sobj.markers <- FindAllMarkers(sobj_combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = opt$verbose)
# top-2 markers per cluster, print them
sobj.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
sobj_combined@misc$markers_pat_type <- sobj.markers
# plot top-10 markers per cluster
top10 <- sobj.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) %>% arrange(factor(cluster, levels = pat_type_levels))
sobj_combined@misc$plots$markers_pat_type <- DoHeatmap(subset(sobj_combined, downsample = 1000), features = top10$gene, group.by='pat_type', group.colors=c(colors$green, colors$yellow, colors$red), angle=0) + NoLegend()
my_ggsave(sobj_combined@misc$plots$markers_pat_type, filename = file_path(full_fig_dir, "markers_pat_type.png"))

# Setting cluster identities back
Idents(object = sobj_combined) <- cluster_idents

###################### FINISH
cat("      - Saving\n")

# decreasing the object size for subsequent RAM intense computations
sobj_combined@misc$plots <- NULL
sobj_combined@misc$markers_pat_type <- NULL
sobj_combined@misc$markers_pat_id <- NULL

if(opt$serialize==T) save(sobj_combined, file=file_path(full_data_dir, "/seurat_object.RData"))

# unmodified data
write.table(sobj_combined@assays$RNA@counts, file=file_path(full_data_dir, "/raw_data.tsv"), sep="\t", row.names=T, col.names=colnames(sobj_combined@assays$RNA@counts), quote=F)

# integrated data
write.table(sobj_combined@assays$integrated@data, file=file_path(full_data_dir, "/corrected_data.tsv"), sep="\t", row.names=T, col.names=colnames(sobj_combined@assays$integrated@data), quote=F)

# cell info
write.table(sobj_combined@meta.data, file=file_path(full_data_dir, "/cells_metadata.tsv"), sep="\t", row.names=T, col.names=T, quote=F)

# patient type merged info
for (i in 1:3){
  
  type <- pat_type_levels[i]
  
  curr_sobj <- subset(sobj_combined, subset = pat_type == type)
  
  if(opt$serialize==T) save(curr_sobj, file=file_path(full_data_dir, sprintf("/seurat_object_%s_type.RData", type)))
  
  # unmodified data
  write.table(curr_sobj@assays$RNA@counts, file=file_path(full_data_dir, sprintf("/raw_data_%s_type.tsv", type)), sep="\t", row.names=T, col.names=colnames(curr_sobj@assays$RNA@counts), quote=F)
  
  # cell info
  write.table(curr_sobj@meta.data, file=file_path(full_data_dir, sprintf("/cells_metadata_%s_type.tsv", type)), sep="\t", row.names=T, col.names=T, quote=F)
  
  
}

##################### CELL-TYPE SPECIFIC PROCESSING

for (i in 1:ncol(meta)){
  
  curr_type <- colnames(meta)[i]
  cat("  > Processing cell type: ", curr_type, " (", i, " of ", ncol(meta), ")\n", sep="")
  
  tryCatch({
    
    # create sample directory
    sample_dir <- file_path(cell_type_dir, curr_type)
    dir.create(sample_dir, showWarnings = F)
    sample_fig_dir <- file_path(sample_dir, 'figs/Seurat')
    dir.create(sample_fig_dir, showWarnings = F, recursive=T)
    sample_data_dir <- file_path(sample_dir, 'data/Seurat')
    dir.create(sample_data_dir, showWarnings = F, recursive=T)
    
    # subsetting only data with a cell type
    curr_sobj <- subset(x=sobj_combined, subset = cell_type == curr_type)
    
    # run typical analysis on the integrated assay
    DefaultAssay(curr_sobj) <- "integrated"
    
    ###################### SCALE DATA
    
    curr_sobj <- ScaleData(object = curr_sobj, verbose = opt$verbose)
    
    ###################### DIM REDUCTION
    
    ## PCA
    curr_sobj <- RunPCA(object = curr_sobj, features = VariableFeatures(object = curr_sobj), verbose = F)
    # loadings
    curr_sobj@misc$plots$pca_loadings <- VizDimLoadings(curr_sobj, dims = 1:4, reduction = "pca")
    my_ggsave(curr_sobj@misc$plots$pca_loadings, filename = file_path(sample_fig_dir, "pca_loadings.png"))
    # elbow
    curr_sobj@misc$plots$pca_elbow <- ElbowPlot(curr_sobj, ndims = 40)
    my_ggsave(curr_sobj@misc$plots$pca_elbow, filename = file_path(sample_fig_dir, "pca_elbow.png"))
    # PCA scatter plot and PC heatmap plot
    curr_sobj@misc$plots$pca_main <- DimPlot(object = curr_sobj, reduction = "pca")
    my_ggsave(curr_sobj@misc$plots$pca_main, filename = file_path(sample_fig_dir, "pca_main.png"))
    curr_sobj@misc$plots$pca_heatmap <- DimHeatmap(curr_sobj, dims = 1:4, ncol=2, nfeatures = 30, cells = 1000, balanced = TRUE, fast = F)
    my_ggsave(curr_sobj@misc$plots$pca_heatmap, filename = file_path(sample_fig_dir, "pca_heatmap.png"))
    
    # # number of components
    # curr_sobj <- JackStraw(curr_sobj, num.replicate = 100, dims=20)
    # curr_sobj <- ScoreJackStraw(curr_sobj, dims = 1:20)
    # JackStrawPlot(curr_sobj, dims = 1:20)
    
    ## T-SNE
    curr_sobj <- RunTSNE(object = curr_sobj, dims = 1:20, verbose = opt$verbose)
    curr_sobj@misc$plots$tsne_main <- DimPlot(object = curr_sobj, reduction = "tsne")
    my_ggsave(curr_sobj@misc$plots$tsne_main, filename = file_path(sample_fig_dir, "/tsne_main.png"))
    
    ## U-MAP
    curr_sobj <- RunUMAP(object = curr_sobj, dims=1:20, verbose = opt$verbose)
    curr_sobj@misc$plots$umap_main <- DimPlot(object = curr_sobj, reduction = "umap")
    my_ggsave(curr_sobj@misc$plots$umap_main, filename = file_path(sample_fig_dir, "/umap_main.png"))
    
    
    ###################### FIND CLUSTERS
    
    curr_sobj <- FindNeighbors(object = curr_sobj, verbose = opt$verbose)
    curr_sobj <- FindClusters(object = curr_sobj, verbose = opt$verbose)
    
    # grouping by clusters
    curr_sobj@misc$plots$pca_clusters <- DimPlot(object = curr_sobj, reduction = "pca")
    my_ggsave(curr_sobj@misc$plots$pca_clusters, filename = file_path(sample_fig_dir, "pca_clusters.png"))
    curr_sobj@misc$plots$tsne_clusters <- DimPlot(object = curr_sobj, reduction = "tsne")
    my_ggsave(curr_sobj@misc$plots$tsne_clusters, filename = file_path(sample_fig_dir, "tsne_clusters.png"))
    curr_sobj@misc$plots$umap_clusters <- DimPlot(object = curr_sobj, reduction = "umap")
    my_ggsave(curr_sobj@misc$plots$umap_clusters, filename = file_path(sample_fig_dir, "umap_clusters.png"))
    
    # grouping by patient and patient type
    p1 <- DimPlot(curr_sobj, reduction = "pca", group.by = "pat_type")
    p2 <- DimPlot(curr_sobj, reduction = "pca", group.by = 'pat_id', label = F, repel = F)
    curr_sobj@misc$plots$pca_pat <- p1 + p2
    my_ggsave(curr_sobj@misc$plots$pca_pat, filename = file_path(sample_fig_dir, "pca_pat.png"))
    
    p1 <- DimPlot(curr_sobj, reduction = "tsne", group.by = "pat_type")
    p2 <- DimPlot(curr_sobj, reduction = "tsne", group.by = 'pat_id', label = F, repel = F)
    curr_sobj@misc$plots$tsne_pat <- p1 + p2
    my_ggsave(curr_sobj@misc$plots$tsne_pat, filename = file_path(sample_fig_dir, "tsne_pat.png"))
    
    p1 <- DimPlot(curr_sobj, reduction = "umap", group.by = "pat_type")
    p2 <- DimPlot(curr_sobj, reduction = "umap", group.by = 'pat_id', label = F, repel = F)
    curr_sobj@misc$plots$umap_pat <- p1 + p2
    my_ggsave(curr_sobj@misc$plots$umap_pat, filename = file_path(sample_fig_dir, "UMAP_pat.png"))
    
    curr_sobj@misc$plots$umap_pat_split <- DimPlot(sobj_combined, reduction = "umap", split.by = "pat_id", group.by = 'pat_type', cols=c(colors$green, colors$yellow, colors$red)) + 
      theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7)) + labs(title="Patient ID")
    my_ggsave(curr_sobj@misc$plots$umap_pat_split, filename = file_path(sample_fig_dir, "umap_pat_split.png"))
    
    ###################### GENE MARKERS
    
    # changing assay to unmodified data
    DefaultAssay(curr_sobj) <- "RNA"
    
    curr_sobj <- ScaleData(object = curr_sobj, verbose = opt$verbose)
    
    # changing the identities to pat_id
    cluster_idents <- Idents(object = curr_sobj)
    Idents(object = curr_sobj) <- curr_sobj@meta.data$pat_id
    
    # find markers for every patient compared to all remaining cells, report only the positive ones
    sobj.markers <- FindAllMarkers(curr_sobj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = opt$verbose)
    # top-2 markers per cluster, print them
    sobj.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
    curr_sobj@misc$markers_pat_id <- sobj.markers
    # plot top-10 markers per cluster
    top10 <- sobj.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
    curr_sobj@misc$plots$markers_pat_id <- DoHeatmap(subset(curr_sobj, downsample = 1000), features = top10$gene) + NoLegend()
    my_ggsave(curr_sobj@misc$plots$markers_pat_id, filename = file_path(sample_fig_dir, "markers_pat_id.png"))
    
    # changing the identities to pat type
    Idents(object = curr_sobj) <- curr_sobj@meta.data$pat_type
    
    # find markers for every patient group compared to all remaining cells, report only the positive ones
    sobj.markers <- FindAllMarkers(curr_sobj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = opt$verbose)
    # top-2 markers per cluster, print them
    sobj.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
    curr_sobj@misc$markers_pat_type <- sobj.markers
    # plot top-10 markers per cluster
    top10 <- sobj.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) %>% arrange(factor(cluster, levels = pat_type_levels))
    curr_sobj@misc$plots$markers_pat_type <- DoHeatmap(subset(curr_sobj, downsample = 1000), features = top10$gene) + NoLegend()
    my_ggsave(curr_sobj@misc$plots$markers_pat_type, filename = file_path(sample_fig_dir, "markers_pat_type.png"))
    
    # setting cluster identities back
    Idents(object = curr_sobj) <- cluster_idents
    
    ###################### FINISH
    cat("      - Saving\n")
    
    ################# Load Seurat object from previous run ####################
    # load(file_path(sample_data_dir, "seurat_object.RData"))
    ###########################################################################
    
    if (opt$serialize==T) save(curr_sobj, file=file_path(sample_data_dir, "seurat_object.RData"))
    
    # unmodified data
    write.table(curr_sobj@assays$RNA@counts, file=file_path(sample_data_dir, "raw_data.tsv"), sep="\t", row.names=T, col.names=colnames(curr_sobj@assays$RNA@counts), quote=F)
    
    # integrated data
    write.table(curr_sobj@assays$integrated@data, file=file_path(sample_data_dir, "corrected_data.tsv"), sep="\t", row.names=T, col.names=colnames(curr_sobj@assays$integrated@data), quote=F)
    
    # cell info
    write.table(curr_sobj@meta.data, file=file_path(sample_data_dir, "cells_metadata.tsv"), sep="\t", row.names=T, col.names=T, quote=F)
    
    # patient type merged info (only control/mild/severe patients)
    for (k in 1:3){
      
      type <- pat_type_levels[k]
      tryCatch({
        curr_subsobj <- subset(curr_sobj, subset = pat_type == type)
        
        if(opt$serialize==T) save(curr_subsobj, file=file_path(sample_data_dir, sprintf("seurat_object_%s_type.RData", type)))
        
        # unmodified data
        write.table(curr_subsobj@assays$RNA@counts, file=file_path(sample_data_dir, sprintf("raw_data_%s_type.tsv", type)), sep="\t", row.names=T, col.names=colnames(curr_subsobj@assays$RNA@counts), quote=F)
        
        # cell info
        write.table(curr_subsobj@meta.data, file=file_path(sample_data_dir, sprintf("cells_metadata_%s_type.tsv", type)), sep="\t", row.names=T, col.names=T, quote=F)
      },
        error = function(e) {
          cat(sprintf('      - No "%s" cells from patient type: "%s"..\n', curr_type, type))
      }
      )
    }
  },
   
  error = function(e) {
    cat(paste0("Encountered error: '", e, "' when processing cell type: ", curr_type, '\n'))
    }
  )
}

cat("\n\n[finished]\n\n")
