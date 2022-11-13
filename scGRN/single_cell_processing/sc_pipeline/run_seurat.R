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
  library(foreach)
  library(doParallel)
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

# per-patient execution function
run <- function(i){
  
  cat("  > Processing sample ", meta$id[i], " (", i, " of ", nrow(meta), ")\n", sep="")
  
  ###################### INIT
  cat("      - Init\n")
  
  # init Seurat object
  raw_counts <- Read10X_h5(meta$file[i], use.names = T, unique.features = T)
  sobj <- CreateSeuratObject(counts = raw_counts, min.cells = 3, min.features = 200, project = 'GRN inference')
  
  # adding some meta data about patients
  sobj@meta.data$pat_id <- meta$id[i]
  sobj@meta.data$pat_type <- pat_types[meta$id[i]]
  
  # create sample directory
  sample_dir <- file_path(opt$outdir, meta$id[i])
  dir.create(sample_dir, recursive = T, showWarnings = F)
  sample_fig_dir <- file_path(sample_dir, 'figs/Seurat')
  dir.create(sample_fig_dir, recursive = T, showWarnings = F)
  sample_data_dir <- file_path(sample_dir, 'data/Seurat')
  dir.create(sample_data_dir, recursive = T, showWarnings = F)
  
  # Delete previously generated files
  # file.remove(list.files(sample_data_dir, full.names = T, recursive = T))
  # file.remove(list.files(sample_fig_dir, full.names = T, recursive = T))
  
  sobj@misc$plots <- list()
  
  ###################### QC
  cat("      - QC\n")
  
  # calculate percentage of MT gene counts
  sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")
  
  # plots for QC features
  sobj@misc$plots$qc_violin <- VlnPlot(sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = .25)
  my_ggsave(sobj@misc$plots$qc_violin, filename = file_path(sample_fig_dir, "qc_violin.png"))
  
  plot1 <- FeatureScatter(sobj, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(sobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  sobj@misc$plots$qc_cross <- plot1 + plot2
  my_ggsave(sobj@misc$plots$qc_cross, filename = file_path(sample_fig_dir, "qc_cross.png"))
  
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
  
  sobj <- FindVariableFeatures(object = sobj, selection.method = "vst", nfeatures = 2000, verbose = opt$verbose)
  top10 <- head(VariableFeatures(sobj), 10)
  plot1 <- VariableFeaturePlot(sobj)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  sobj@misc$plots$variable_features <- plot1 + plot2
  my_ggsave(sobj@misc$plots$variable_features, filename = file_path(sample_fig_dir, "variable_features.png"))
  
  ###################### SCALE DATA
  cat("      - Scaling\n")
  
  sobj <- ScaleData(object = sobj, verbose = opt$verbose)
  
  ###################### DIM REDUCTION
  cat("      - Dimensionality reduction\n")
  
  # PCA
  sobj <- RunPCA(object = sobj, features = VariableFeatures(object = sobj), verbose = opt$verbose)
  # loadings
  sobj@misc$plots$pca_loadings <- VizDimLoadings(sobj, dims = 1:4, reduction = "pca")
  my_ggsave(sobj@misc$plots$pca_loadings, filename = file_path(sample_fig_dir, "pca_loadings.png"))
  # elbow
  sobj@misc$plots$pca_elbow <- ElbowPlot(sobj, ndims = 40)
  my_ggsave(sobj@misc$plots$pca_elbow, filename = file_path(sample_fig_dir, "pca_elbow.png"))
  # PCA scatter plot and PC heatmap plot
  sobj@misc$plots$pca_main <- DimPlot(object = sobj, reduction = "pca")
  my_ggsave(sobj@misc$plots$pca_main, filename = file_path(sample_fig_dir, "pca_main.png"))
  sobj@misc$plots$pca_heatmap <- DimHeatmap(sobj, dims = 1:4, ncol=2, nfeatures = 30, cells = 1000, balanced = TRUE, fast = F)
  my_ggsave(sobj@misc$plots$pca_heatmap, filename = file_path(sample_fig_dir, "pca_heatmap.png"))
  
  # # number of components
  # sobj <- JackStraw(sobj, num.replicate = 100, dims=20)
  # sobj <- ScoreJackStraw(sobj, dims = 1:20)
  # JackStrawPlot(sobj, dims = 1:20)
  
  ## T-SNE
  sobj <- RunTSNE(object = sobj, dims = 1:20, verbose = opt$verbose)
  sobj@misc$plots$tsne_main <- DimPlot(object = sobj, reduction = "tsne")
  my_ggsave(sobj@misc$plots$tsne_main, filename = file_path(sample_fig_dir, "tsne_main.png"))
  
  ## U-MAP
  sobj <- RunUMAP(object = sobj, dims=1:20, verbose = opt$verbose)
  sobj@misc$plots$umap_main <- DimPlot(object = sobj, reduction = "umap")
  my_ggsave(sobj@misc$plots$umap_main, filename = file_path(sample_fig_dir, "umap_main.png"))
  
  
  ###################### FIND CLUSTERS
  cat("      - Find clusters\n")
  
  sobj <- FindNeighbors(object = sobj, verbose = opt$verbose)
  sobj <- FindClusters(object = sobj, verbose = opt$verbose)
  
  # plotting clusters using PCA, T-SNE, UMAP
  sobj@misc$plots$pca_clusters <- DimPlot(object = sobj, reduction = "pca")
  my_ggsave(sobj@misc$plots$pca_clusters, filename = file_path(sample_fig_dir, "pca_clusters.png"))
  sobj@misc$plots$tsne_clusters <- DimPlot(object = sobj, reduction = "tsne")
  my_ggsave(sobj@misc$plots$tsne_clusters, filename = file_path(sample_fig_dir, "tsne_clusters.png"))
  sobj@misc$plots$umap_clusters <- DimPlot(object = sobj, reduction = "umap")
  my_ggsave(sobj@misc$plots$umap_clusters, filename = file_path(sample_fig_dir, "umap_clusters.png"))
  
  
  ###################### GENE MARKERS
  cat("      - Gene markers by cluster\n")
  
  # find markers for every cluster compared to all remaining cells, report only the positive ones
  sobj.markers <- FindAllMarkers(sobj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = opt$verbose)
  # top-2 markers per cluster, print them
  sobj.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
  sobj@misc$markers <- sobj.markers
  # plot top-10 markers per cluster
  top10 <- sobj.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  sobj@misc$plots$markers <- DoHeatmap(sobj, features = top10$gene) + NoLegend()
  my_ggsave(sobj@misc$plots$markers, filename = file_path(sample_fig_dir, "markers.png"))
  
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
  
  # Prediction
  pred <- SingleR(test=sobj_se, ref=ref, labels=ref$label.main)
  sobj@misc$cell_type_ref <- ref
  sobj@misc$cell_type_pred <- pred
  sobj@misc$plots$cell_type_pred_dist <- plotScoreDistribution(pred, size = .25)
  my_ggsave(sobj@misc$plots$cell_type_pred_dist, filename = file_path(sample_fig_dir, "cell_type_pred_dist.png"))
  sobj@misc$plots$cell_type_pred_heatmap <- plotScoreHeatmap(pred,show.pruned = T)
  my_ggsave(sobj@misc$plots$cell_type_pred_heatmap, filename = file_path(sample_fig_dir, "cell_type_pred_heatmap.png"))
  
  # Cell type freqs
  png(file_path(sample_fig_dir, "cell_type_freqs.png"), width=800, height=700, pointsize = 20)
  par(mar=c(5,10,3,3))
  freqs <- sort(table(pred$labels), decreasing = F)
  barplot(freqs, horiz=T, las=2)
  
  # cell type to clusters
  pheatmap(log(10+table(pred$labels, sobj$seurat_clusters)), filename = file_path(sample_fig_dir, "cell_type_to_cluster.png"))
  
  # dim reduction plots
  sobj[["cell_type"]] <- pred$labels
  sobj@misc$plots$pca_cell_types <- DimPlot(object = sobj, reduction = "pca", group.by = "cell_type")
  my_ggsave(sobj@misc$plots$pca_cell_types, filename = file_path(sample_fig_dir, "pca_cell_types.png"))
  sobj@misc$plots$tsne_cell_types <- DimPlot(object = sobj, reduction = "tsne", group.by = "cell_type")
  my_ggsave(sobj@misc$plots$tsne_cell_types, filename = file_path(sample_fig_dir, "tsne_cell_types.png"))
  sobj@misc$plots$umap_cell_types <- DimPlot(object = sobj, reduction = "umap", group.by = "cell_type")
  my_ggsave(sobj@misc$plots$umap_cell_types, filename = file_path(sample_fig_dir, "umap_cell_types.png"))
  
  
  # summarize by cell type
  summ_fun <- function(x) mean(x, trim=.005)
  summ_cells <- function(x, fun=summ_fun){
    # print(dim(x))
    apply(x,2,fun)
  }
  sobj@misc$cell_type_assay <- do.call("cbind",by(t(sobj@assays$RNA@data), pred$labels, summ_cells))
  
  
  ###################### FINISH
  cat("      - Saving\n")
  
  if(opt$serialize==T) save(sobj, file=file_path(sample_data_dir,"seurat_object.RData"))

  # prune labels - keep only cells with high confidence annotation
  cell_types <- unique(pred$labels)
  sobj_to_save <- sobj
  
  # data
  #   general
  write.table(sobj_to_save@assays$RNA@data, file=file_path(sample_data_dir, "norm_data.tsv"), sep="\t", row.names=T, col.names=colnames(sobj_to_save@assays$RNA@data), quote=F)
  write.table(sobj_to_save@assays$RNA@counts, file=file_path(sample_data_dir, "raw_data.tsv"), sep="\t", row.names=T, col.names=colnames(sobj_to_save@assays$RNA@counts), quote=F)
  write.table(sobj_to_save@assays$RNA@scale.data, file=file_path(sample_data_dir, "scaled_data.tsv"), sep="\t", row.names=T, col.names=colnames(sobj_to_save@assays$RNA@scale.data), quote=F)
  #   cell info
  write.table(sobj_to_save@meta.data, file=file_path(sample_data_dir, "cells_metadata.tsv"), sep="\t", row.names=T, col.names=T, quote=F)
  
  #   gene-specific
  HVG_norm_data = sobj_to_save@assays$RNA@data[rownames(sobj_to_save@assays$RNA@data) %in% VariableFeatures(sobj_to_save),]  # saving norm_data with HVG only  
  write.table(HVG_norm_data, file=file_path(sample_data_dir, "norm_data_HVG.tsv"), sep="\t", row.names=T, col.names=colnames(HVG_norm_data), quote=F)
  HVG_norm_counts = sobj_to_save@assays$RNA@counts[rownames(sobj_to_save@assays$RNA@counts) %in% VariableFeatures(sobj_to_save),]  # saving norm_data with HVG only  
  write.table(HVG_norm_counts, file=file_path(sample_data_dir, "raw_data_HVG.tsv"), sep="\t", row.names=T, col.names=colnames(HVG_norm_counts), quote=F)
  
  #   cell type-specific
  paths_out <- c()  # saving paths to cell type-specific matrices
  for (k in 1:length(cell_types)){
    type <- cell_types[k]
    type_fn <- file_path(sample_data_dir, sprintf("raw_data_%s.tsv", type))
    cell_type_sobj_to_save <- subset(x=sobj_to_save, subset = cell_type == type)
    write.table(cell_type_sobj_to_save@assays$RNA@counts, file=type_fn, sep="\t", row.names=T, col.names=colnames(cell_type_sobj_to_save@assays$RNA@counts), quote=F)
    # write.table(cell_type_sobj_to_save@assays$RNA@data, file=file_path(sample_data_dir, sprintf("norm_data_%s.tsv", type)), sep="\t", row.names=T, col.names=colnames(cell_type_sobj_to_save@assays$RNA@data), quote=F)
    paths_out[type] <- type_fn
  }
  
  # Returning metadata with file paths
  paths_out
}

###################### INPUT

option_list = list(
  make_option(c("-m", "--meta_file"), type="character", help="Metadata file", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", help="Output folder", metavar="character"),
  make_option(c('-a', '--annotation'), type='character', help='Annotation dataset used to label cells', default='HumanPrimaryCellAtlasData', metavar='character'),
  make_option(c('--annotation_folder'), type='character', help='Folder with saved annotation datasets', default=NULL, metavar='character'),
  make_option(c('-n', '--num_proc'), type='integer', help='Number of processes run in parallel', default=6, metavar='integer'),
  make_option(c("-s", "--serialize"), type="logical", default=T, help="Save Seurat object", metavar="T|F"),
  make_option(c('--add_figure_objects'), type='logical', default=T, help='Whether to add figure objects to Seurat object', metavar='T|F'),
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
  opt$num_proc = 6
}
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

cat("\n\n")
cat("***********************************\n")
cat("*** SINGLE-CELL DATA PROCESSING ***\n")
cat("***********************************\n\n")
cat("Outdir: ", opt$outdir, "\n")
cat("Metadata: ", opt$meta_file, "\n")
cat("Annotation: ", opt$annotation, "\n")
cat("Annotation folder: ", opt$annotation_folder, "\n")
cat("Parallelized by: ", opt$num_proc, '\n')
cat("Serialize output: ", opt$serialize, "\n")
cat("Save figure objects: ", opt$save_figure_objects, '\n')
cat("Verbose: ", opt$verbose, "\n")

###################### LOAD DATA

# read metadata
meta <- read.table(opt$meta_file, header=T, sep="\t", stringsAsFactors = F)
rownames(meta) <- meta$id
pat_types <- meta$group
names(pat_types) <- meta$id

# create working directory
dir.create(opt$outdir, recursive = T, showWarnings = F)

# process samples
sobjs <- list()
cell_type_path_meta <- data.frame(row.names=meta$id)

# setting up the cluster for parallel execution
cl <- makeCluster(opt$num_proc)
registerDoParallel(cl)

# run everything in parallel, obtain cell types
result <- foreach(i=1:nrow(meta), .packages=c('SingleR', 'dplyr', 'Matrix', 'Seurat', 'future', 'pheatmap', 'ggplot2')) %dopar% run(i)

stopCluster(cl)

# filling up the cell type meta data (with pathname info)
for (i in 1:nrow(meta)){
  for (type in names(result[[i]])){
    type_fn <- result[[i]][[type]]
    cell_type_path_meta[meta$id[i], type] = type_fn
  }
}

# saving meta file with cell types file paths
write.table(cell_type_path_meta, file=file_path(opt$outdir, 'cell_type_meta.tsv'), row.names=T, sep='\t', col.names=T, quote=F)

cat("\n\n[finished]\n\n")
