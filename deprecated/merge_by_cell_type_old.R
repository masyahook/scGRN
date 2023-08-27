rm(list=ls())

suppressPackageStartupMessages({
  # if (!require("pacman")) install.packages("pacman")
  # list.of.packages <- c("BiocManager","dplyr","SingleR","Matrix","Seurat","future","pheatmap","ggplot2","optparse","hdf5r")
  # BiocManager::install("SingleR")
  # BiocManager::install('limma')
  # BiocManager::install('SingleCellExperiment')
  # pacman::p_load(list.of.packages, character.only = TRUE)
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
})

###################### INPUT

option_list = list(
  make_option(c("-m", "--meta_file"), type="character", default='/gpfs/projects/bsc08/bsc08890/res/covid_19/cell_type_meta.tsv', help="Metadata file with cell types", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default='/gpfs/projects/bsc08/bsc08890/res/covid_19/cell_types', help="Output folder with cell-type specific data", metavar="character"),
  make_option(c('-t', "--patient_type"), type='character', default='/gpfs/projects/bsc08/bsc08890/data/GSE145926_RAW/metadata.tsv', help='Metadata file with patient disease types', metavar='character'),
  make_option(c("-v", "--verbose"), type="logical", default=T, help="Verbose", metavar="T|F"),
  make_option(c("-s", "--serialize"), type="logical", default=F, help="Save Seurat object", metavar="T|F")
)
opt_parser <- OptionParser(option_list=option_list, add_help_option = T)
opt <- parse_args(opt_parser)

# EVALUATE CLI
# if (is.null(opt$meta_file)){
#   print_help(opt_parser)
#   stop("No metadata file with cell types provided", call.=FALSE)
# }
# if (is.null(opt$outdir)){
#   print_help(opt_parser)
#   stop("No output folder provided", call.=FALSE)
# }
# if (is.null(opt$patient_type)){
#   print_help(opt_parser)
#   stop("No metadata file with patient disease types provided", call.=FALSE)
# }

# # DEBUG
# setwd("~")
# opt <- list()
# opt$meta_file = "/gpfs/projects/bsc08/bsc08890/res/covid_19/cell_type_meta.tsv"
# opt$outdir = "/gpfs/projects/bsc08/bsc08890/res/covid_19/cell_types"
# opt$patient_type = '/gpfs/projects/bsc08/bsc08890/data/GSE145926_RAW/metadata.tsv'


cat("\n\n")
cat("***********************************\n")
cat("*** MERGE DATASETS BY CELL TYPE ***\n")
cat("***********************************\n\n")
cat("metadata file: ", opt$meta_file, "\n")
cat("outdir: ", opt$outdir, "\n")

# read metadata
meta <- read.table(opt$meta_file, header=T, sep="\t", stringsAsFactors = F, check.names=F, na.strings = '')

# read patient disease types
pat_meta <- read.table(opt$patient_type, header=T, sep='\t', stringsAsFactors = F)
pat_types <- pat_meta$group
names(pat_types) <- pat_meta$id

# create directory where we will save cell-type specific data
dir.create(opt$outdir, showWarnings = F)

my_ggsave <- function(obj, filename){
  suppressMessages(ggsave(obj, filename = filename))
}

for(i in 1:ncol(meta)){
  
  cell_type <- colnames(meta)[i]
  cat("  > Processing cell type: ",cell_type," (",i," of ", ncol(meta),")\n", sep="")
  
  ###################### INIT
  cat("      - Init\n")
  
  # create sample directory
  sample_dir <- paste0(opt$outdir, "/", cell_type)
  dir.create(sample_dir, showWarnings = F)
  sample_fig_dir <- paste0(sample_dir, '/figs/Seurat')
  dir.create(sample_fig_dir, showWarnings = F, recursive=T)
  sample_data_dir <- paste0(sample_dir, '/data/Seurat')
  dir.create(sample_data_dir, showWarnings = F, recursive=T)
  
  # Find patients that have this type of cells
  pat_fns <- na.omit(meta[cell_type])
  
  # Reading data from all patients with this type of cells
  pat_objs <- list()
  num_cells_per_pat <- c()
  for (pat_id in rownames(pat_fns)){
    
    # Define paths to files
    curr_fn <- pat_fns[pat_id, cell_type]
    curr_meta_fn <- paste0(dirname(curr_fn), '/cells_metadata.tsv')
    
    # Read the data and metadata
    curr_counts <- read.table(file=curr_fn, sep="\t", check.names=F)
    curr_meta <- read.table(file=curr_meta_fn, sep='\t', check.names=F)
    
    # Select in metadata only cell-type specific cells
    curr_meta <- curr_meta[curr_meta$cell_type == cell_type,]
    
    # Rename orig.dent in metadata to patient ID
    curr_meta$orig.ident <- pat_id
    
    # Load everything to Seurat object
    curr_obj <- CreateSeuratObject(counts = curr_counts, min.cells = 0, min.genes = 0, project = pat_id, meta.data = curr_meta)
    
    # Save into the list
    pat_objs[pat_id] <- curr_obj
    
    # Getting number of cells per patient
    num_cells_per_pat[pat_id] <- nrow(curr_meta)
  }
  
  # If the number of cells is too low, then we will drop these patients
  if (any(num_cells_per_pat < 200) == T){
    small_pats <- names(num_cells_per_pat[num_cells_per_pat < 30])
    pat_objs <- pat_objs[!(names(pat_objs) %in% small_pats)]
    num_cells_per_pat <- num_cells_per_pat[!(names(num_cells_per_pat) %in% small_pats)]
  }
  
  ############# DATA INTEGRATION
  cat("      - Data integration\n")
  # normalize and identify variable features for each dataset independently
  pat_objs <- lapply(X = pat_objs, FUN = function(x) {
    x <- NormalizeData(x, verbose = F)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = F)
  })
  
  # select features that are repeatedly variable across datasets for integration
  features <- SelectIntegrationFeatures(object.list = pat_objs)
  anchors <- FindIntegrationAnchors(object.list = pat_objs, anchor.features = features, verbose = F, dims=1:15)
  sobj_combined <- IntegrateData(anchorset = anchors, verbose = F, dims = 1:15, k.weight = 50)
  
  # specify that we will perform downstream analysis on the corrected data note that the
  # original unmodified data still resides in the 'RNA' assay
  DefaultAssay(sobj_combined) <- "integrated"
  
  ###################### VARIABLE FEATURES
  cat("      - Find variable features\n")
  
  sobj_combined <- FindVariableFeatures(object = sobj_combined, selection.method = "vst", nfeatures = 2000, verbose = F)
  top10 <- head(VariableFeatures(sobj_combined), 10)
  plot1 <- VariableFeaturePlot(sobj_combined)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  sobj_combined@misc$plots$variable_features <- plot1 + plot2
  my_ggsave(sobj_combined@misc$plots$variable_features, filename = paste0(sample_fig_dir, "/variable_features.png"))
  
  
  ###################### SCALE DATA
  cat("      - Scaling\n")
  
  sobj_combined <- ScaleData(object = sobj_combined, verbose = F)
  
  
  ###################### DIM REDUCTION
  cat("      - Dimensionality reduction\n")
  
  ## PCA
  sobj_combined <- RunPCA(object = sobj_combined, features = VariableFeatures(object = sobj_combined), verbose = F)
  # loadings
  sobj_combined@misc$plots$pca_loadings <- VizDimLoadings(sobj_combined, dims = 1:4, reduction = "pca")
  my_ggsave(sobj_combined@misc$plots$pca_loadings, filename = paste0(sample_fig_dir, "/pca_loadings.png"))
  # elbow
  sobj_combined@misc$plots$pca_elbow <- ElbowPlot(sobj_combined, ndims = 40)
  my_ggsave(sobj_combined@misc$plots$pca_elbow, filename = paste0(sample_fig_dir, "/pca_elbow.png"))
  #
  sobj_combined@misc$plots$pca_main <- DimPlot(object = sobj_combined, reduction = "pca")
  my_ggsave(sobj_combined@misc$plots$pca_main, filename = paste0(sample_fig_dir, "/pca_main.png"))
  sobj_combined@misc$plots$pca_heatmap <- DimHeatmap(sobj_combined, dims = 1:4, ncol=2, nfeatures = 30, cells = 1000, balanced = TRUE, fast = F)
  my_ggsave(sobj_combined@misc$plots$pca_heatmap, filename = paste0(sample_fig_dir, "/pca_heatmap.png"))
  
  # # number of components
  # sobj_combined <- JackStraw(sobj_combined, num.replicate = 100, dims=20)
  # sobj_combined <- ScoreJackStraw(sobj_combined, dims = 1:20)
  # JackStrawPlot(sobj_combined, dims = 1:20)
  
  ## T-SNE
  sobj_combined <- RunTSNE(object = sobj_combined, dims = 1:20, verbose = F)
  sobj_combined@misc$plots$tsne_main <- DimPlot(object = sobj_combined, reduction = "tsne")
  sobj_combined@misc$plots$tsne_pat <- DimPlot(object = sobj_combined, reduction = "tsne", group.by = 'orig.ident')
  my_ggsave(sobj_combined@misc$plots$tsne_main, filename = paste0(sample_fig_dir, "/tsne_main.png"))
  my_ggsave(sobj_combined@misc$plots$tsne_pat, filename = paste0(sample_fig_dir, "/tsne_pat.png"))
  
  ## U-MAP
  sobj_combined <- RunUMAP(object = sobj_combined, dims=1:20, verbose = F)
  sobj_combined@misc$plots$umap_main <- DimPlot(object = sobj_combined, reduction = "umap")
  sobj_combined@misc$plots$umap_pat <- DimPlot(object = sobj_combined, reduction = "umap", group.by = 'orig.ident')
  my_ggsave(sobj_combined@misc$plots$umap_main, filename = paste0(sample_fig_dir, "/umap_main.png"))
  my_ggsave(sobj_combined@misc$plots$umap_pat, filename = paste0(sample_fig_dir, "/umap_pat.png"))
  
  
  
  ###################### FIND CLUSTERS
  cat("      - Find clusters\n")
  
  sobj_combined <- FindNeighbors(object = sobj_combined, verbose = F)
  sobj_combined <- FindClusters(object = sobj_combined, verbose = F)
  
  sobj_combined@misc$plots$pca_clusters <- DimPlot(object = sobj_combined, reduction = "pca")
  my_ggsave(sobj_combined@misc$plots$pca_clusters, filename = paste0(sample_fig_dir, "/pca_clusters.png"))
  sobj_combined@misc$plots$tsne_clusters <- DimPlot(object = sobj_combined, reduction = "tsne")
  my_ggsave(sobj_combined@misc$plots$tsne_clusters, filename = paste0(sample_fig_dir, "/tsne_clusters.png"))
  sobj_combined@misc$plots$umap_clusters <- DimPlot(object = sobj_combined, reduction = "umap")
  my_ggsave(sobj_combined@misc$plots$umap_clusters, filename = paste0(sample_fig_dir, "/umap_clusters.png"))
  
  
  ###################### GENE MARKERS
  cat("      - Gene markers by cluster\n")
  
  # find markers for every cluster compared to all remaining cells, report only the positive ones
  sobj_combined.markers <- FindAllMarkers(sobj_combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = F)
  # XXX
  # sobj_combined.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
  sobj_combined.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
  sobj_combined@misc$markers <- sobj_combined.markers
  # XXX
  # top10 <- sobj_combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  top10 <- sobj_combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  sobj_combined@misc$plots$markers <- DoHeatmap(sobj_combined, features = top10$gene) + NoLegend()
  my_ggsave(sobj_combined@misc$plots$markers, filename = paste0(sample_fig_dir, "/markers.png"))
  
  ###################### GENE MARKERS
  cat("      - Gene markers by patient\n")
  
  # Changing the identities to pat_id
  curr_idents <- Idents(object = sobj_combined)
  Idents(object = sobj_combined) <- sobj_combined@meta.data$orig.ident
  
  # find markers for every cluster compared to all remaining cells, report only the positive ones
  sobj_combined.markers <- FindAllMarkers(sobj_combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = F)
  # XXX
  # sobj_combined.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
  sobj_combined.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
  sobj_combined@misc$markers_pat <- sobj_combined.markers
  # XXX
  # top10 <- sobj_combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  top10 <- sobj_combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  sobj_combined@misc$plots$markers_pat <- DoHeatmap(sobj_combined, features = top10$gene) + NoLegend()
  my_ggsave(sobj_combined@misc$plots$markers_pat, filename = paste0(sample_fig_dir, "/markers_pat.png"))
  
  # Setting back previous identities
  Idents(object = sobj_combined) <- curr_idents
  
  ###################### FINISH
  cat("      - Saving\n")
  
  # TODO: WHY FAILS? (JCB)
  if (opt$serialize==T) save(sobj_combined, file=paste0(sample_data_dir,"/seurat_object.RData"))
  
  # data
  # general
  write.table(sobj_combined@assays$integrated@counts, file=paste0(sample_data_dir, "/raw_data.tsv"), sep="\t", row.names=T, col.names=colnames(sobj_combined@assays$integrated@counts), quote=F)
  # cell info
  write.table(sobj_combined@meta.data, file=paste0(sample_data_dir, "/cells_metadata.tsv"), sep="\t", row.names=T, col.names=T, quote=F)
}

cat("\n\n[finished]\n\n")