rm(list=ls())  # clear namespace

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
  library(dorothea)
  library(tibble)
  library(tidyr)
  library(rlang)
  library(reticulate)
})

pd <- import('pandas')

# save figures without messages
my_ggsave <- function(obj, filename){
  suppressMessages(ggsave(obj, filename = filename))
}

# deal with duplicate slashes
file_path = function(..., fsep = .Platform$file.sep){
  gsub("//", "/", file.path(..., fsep = fsep))
}

###################### INPUT

option_list = list(
  make_option(c("-m", "--meta_file"), type="character", default='/gpfs/projects/bsc08/bsc08890/res/covid_19/cell_type_meta.tsv', help="Metadata file with cell types", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default='/gpfs/projects/bsc08/bsc08890/res/covid_19', help="Output folder with cell-type specific data", metavar="character"),
  make_option(c('-r', '--regulon'), type='character', default='pyscenic', help='Regulons to use - either dorothea or pyscenic', metavar='character'),
  make_option(c('-q', '--quantile'), type='character', default='', help='Quantile threshold to search for the network', metavar='character'),
  make_option(c("-c", "--pleiotropy_correction"), type="logical", default=T, help='Pleitropy correction?', metavar="T|F"),
  make_option(c('-n', '--num_proc'), type='integer', help='Number of processes run in parallel', default=6, metavar='integer'),
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
if (opt$regulon == ''){
  opt$regulon = 'pyscenic'
}
if (opt$pleiotropy_correction == ''){
  opt$pleiotropy_correction = T
}
if (is.na(opt$num_proc)){
  opt$num_proc = 6
}
if (is.na(opt$serialize)){
  opt$serialize = T
}
if (is.na(opt$verbose)){
  opt$verbose = F
}

cat("\n\n")
cat("*************************************\n")
cat("*** VIPER PROCESSING BY CELL TYPE ***\n")
cat("*************************************\n\n")
cat("Outdir: ", opt$outdir, "\n")
cat("Metadata: ", opt$meta_file, "\n")
cat("Regulon type: ", opt$regulon, "\n")
cat("Network quantile threshold filter: ", opt$quantile, "\n")
cat("Apply pleiotropy correction: ", opt$pleiotropy_correction, "\n")
cat("Verbose: ", opt$verbose, "\n")


###################### LOAD DATA

# read metadata
meta <- read.table(opt$meta_file, header=T, sep="\t", stringsAsFactors = F, check.names=F, na.strings = '')
pat_type_levels <- c("C", "M", "S")

# create directory where we will save cell-type specific data
cell_type_dir <- paste0(opt$outdir, '/cell_types')
dir.create(cell_type_dir, showWarnings = F)

# getting DoRothEA
dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))

# getting dorothea regulons if chosen
if (opt$regulon == 'dorothea'){
  regulon <- dorothea_regulon_human %>%
    dplyr::filter(confidence %in% c("A","B","C"))
}

colors <- list(green='#39B600', yellow='#D89000', red='#F8766D', blue='#00B0F6', 
               purple='#9590FF', cyan='#00BFC4', pink='E76BF3', light_pink='#FF62BC',
               saturated_green='#00BF7D')

##################### CELL-TYPE SPECIFIC PROCESSING

for(i in 1:ncol(meta)){
  
  tryCatch({
    curr_type <- colnames(meta)[i]
    cat("  > Processing cell type: ",curr_type," (",i," of ", nrow(meta),")\n", sep="")
    
    ###################### INIT
    cat("      - Init\n")
    
    # Specifying data/fig folders
    sample_dir <- paste0(cell_type_dir, "/", curr_type)
    sample_fig_dir <- paste0(sample_dir, '/figs/Seurat')
    sample_data_dir <- paste0(sample_dir, '/data/Seurat')
    
    if (opt$regulon == 'pyscenic' | opt$regulon == 'pyscenic_dorothea'){
      thresh <- 0.1
      pyscenic_regulon_dir <- paste0(sample_dir, '/data/grnboost2/pickle')
      if (opt$quantile != ''){
        regulon <- pd$read_pickle(paste0(pyscenic_regulon_dir, sprintf('/raw_data_TF_ctx_filtered_%s.pickle', opt$quantile)))
      } else {
        regulon <- pd$read_pickle(paste0(pyscenic_regulon_dir, '/raw_data_TF_ctx.pickle'))
      }
      regulon <- regulon[,c('TF', 'importance', 'target', 'rho')]
      colnames(regulon) <- c('tf', 'importance', 'target', 'rho')
      if (opt$regulon == 'pyscenic') {
        regulon$mor <- regulon$rho
        regulon$mor[regulon$mor > thresh] <- 1
        regulon$mor[regulon$mor < -thresh] <- -1
        regulon <- regulon[regulon$mor %in% c(-1, 1),] %>% as_tibble()
      } else if (opt$regulon == 'pyscenic_dorothea') {
        regulon <- regulon %>%
          merge(dorothea_regulon_human, by.x=c('tf', 'target'), by.y=c('tf', 'target')) %>%
          as_tibble()
      }
    }
    
    # Delete all previously generated files
    # file.remove(list.files(sample_data_dir, full.names = T, recursive = T, pattern='.*dorothea.*'))
    
    ###################### VIPER BASED ON REGULONS
    
    load(paste0(sample_data_dir, "/seurat_object.RData"))
  
    # Running viper  
    cat("      - Running viper\n")
    
    if (nrow(curr_sobj) > 10000) {
      cores = 46
      
    } else {
      cores = 46
    }
    curr_sobj <- run_viper(curr_sobj, regulon,
                      options = list(method = "scale", minsize = 4, 
                                     pleiotropy = opt$pleiotropy_correction, 
                                     eset.filter = FALSE, 
                                     cores = cores, verbose = FALSE),
                      )
    
    # Running viper-based pipeline
    if (opt$regulon == 'pyscenic' | opt$regulon == 'pyscenic_dorothea'){
      curr_sobj <- RenameAssays(object = curr_sobj, 'dorothea'=opt$regulon)
    }
    DefaultAssay(object = curr_sobj) <- opt$regulon
    
    ###################### SCALE DATA
    cat("      - Scaling\n")
    
    curr_sobj <- ScaleData(object = curr_sobj, verbose = F)
    
    ###################### DIM REDUCTION
    cat("      - Dimensionality reduction\n")
    
    ## PCA
    curr_sobj <- RunPCA(object = curr_sobj, features = rownames(curr_sobj), verbose = F)
    # loadings
    curr_sobj@misc$plots[[sprintf('%s_pca_loadings', opt$regulon)]] <- VizDimLoadings(curr_sobj, dims = 1:4, reduction = "pca")
    my_ggsave(curr_sobj@misc$plots[[sprintf('%s_pca_loadings', opt$regulon)]], filename = paste0(sample_fig_dir, sprintf("/%s_pca_loadings.png", opt$regulon)))
    curr_sobj@misc$plots[[sprintf('%s_pca_loadings', opt$regulon)]]
    # elbow
    curr_sobj@misc$plots[[sprintf('%s_pca_elbow', opt$regulon)]] <- ElbowPlot(curr_sobj, ndims = 40)
    my_ggsave(curr_sobj@misc$plots[[sprintf('%s_pca_elbow', opt$regulon)]], filename = paste0(sample_fig_dir, sprintf("/%s_pca_elbow.png", opt$regulon)))
    curr_sobj@misc$plots[[sprintf('%s_pca_elbow', opt$regulon)]]
    #
    curr_sobj@misc$plots[[sprintf('%s_pca_main', opt$regulon)]] <- DimPlot(object = curr_sobj, reduction = "pca")
    my_ggsave(curr_sobj@misc$plots[[sprintf('%s_pca_main', opt$regulon)]], filename = paste0(sample_fig_dir, sprintf("/%s_pca_main.png", opt$regulon)))
    curr_sobj@misc$plots[[sprintf('%s_pca_main', opt$regulon)]]
    curr_sobj@misc$plots[[sprintf('%s_pca_heatmap', opt$regulon)]] <- DimHeatmap(curr_sobj, dims = 1:4, ncol=2, nfeatures = 30, cells = 1000, balanced = TRUE, fast = F)
    my_ggsave(curr_sobj@misc$plots[[sprintf('%s_pca_heatmap', opt$regulon)]], filename = paste0(sample_fig_dir, sprintf("/%s_pca_heatmap.png", opt$regulon)))
    curr_sobj@misc$plots[[sprintf('%s_pca_heatmap', opt$regulon)]]
    
    ## T-SNE
    curr_sobj <- RunTSNE(object = curr_sobj, dims = 1:10, verbose = F)
    curr_sobj@misc$plots[[sprintf('%s_tsne_main', opt$regulon)]] <- DimPlot(object = curr_sobj, reduction = "tsne")
    my_ggsave(curr_sobj@misc$plots[[sprintf('%s_tsne_main', opt$regulon)]], filename = paste0(sample_fig_dir, sprintf("/%s_tsne_main.png", opt$regulon)))
    curr_sobj@misc$plots[[sprintf('%s_tsne_main', opt$regulon)]]
    
    ## U-MAP
    curr_sobj <- RunUMAP(object = curr_sobj, dims=1:10, verbose = F, umap.method = "uwot", metric = "cosine")
    curr_sobj@misc$plots[[sprintf('%s_umap_main', opt$regulon)]] <- DimPlot(object = curr_sobj, reduction = "umap")
    my_ggsave(curr_sobj@misc$plots[[sprintf('%s_umap_main', opt$regulon)]], filename = paste0(sample_fig_dir, sprintf("/%s_umap_main.png", opt$regulon)))
    curr_sobj@misc$plots[[sprintf('%s_umap_main', opt$regulon)]]
    
    ###################### FIND CLUSTERS
    cat("      - Find clusters\n")
    
    curr_sobj <- FindNeighbors(curr_sobj, dims = 1:10, verbose = FALSE)
    curr_sobj <- FindClusters(curr_sobj, resolution = 0.5, verbose = FALSE)
    
    curr_sobj@misc$plots[[sprintf('%s_pca_clusters', opt$regulon)]] <- DimPlot(object = curr_sobj, reduction = "pca")
    my_ggsave(curr_sobj@misc$plots[[sprintf('%s_pca_clusters', opt$regulon)]], filename = paste0(sample_fig_dir, sprintf("/%s_pca_clusters.png", opt$regulon)))
    curr_sobj@misc$plots[[sprintf('%s_pca_clusters', opt$regulon)]]
    curr_sobj@misc$plots[[sprintf('%s_tsne_clusters', opt$regulon)]] <- DimPlot(object = curr_sobj, reduction = "tsne")
    my_ggsave(curr_sobj@misc$plots[[sprintf('%s_tsne_clusters', opt$regulon)]], filename = paste0(sample_fig_dir, sprintf("/%s_tsne_clusters.png", opt$regulon)))
    curr_sobj@misc$plots[[sprintf('%s_tsne_clusters', opt$regulon)]]
    curr_sobj@misc$plots[[sprintf('%s_umap_clusters', opt$regulon)]] <- DimPlot(object = curr_sobj, reduction = "umap")
    my_ggsave(curr_sobj@misc$plots[[sprintf('%s_umap_clusters', opt$regulon)]], filename = paste0(sample_fig_dir, sprintf("/%s_umap_clusters.png", opt$regulon)))
    curr_sobj@misc$plots[[sprintf('%s_umap_clusters', opt$regulon)]]
    
    # grouping by patient and patient type
    p1 <- DimPlot(curr_sobj, reduction = "umap", group.by = "pat_type")
    p2 <- DimPlot(curr_sobj, reduction = "umap", group.by = 'pat_id', label = TRUE, repel = TRUE)
    curr_sobj@misc$plots$umap_pat <- p1 + p2
    my_ggsave(curr_sobj@misc$plots$umap_pat, filename = paste0(sample_fig_dir, sprintf("/%s_UMAP_pat.png", opt$regulon)))
    
    curr_sobj@misc$plots$umap_pat_split <- DimPlot(curr_sobj, reduction = "umap", split.by = "pat_id", group.by = 'pat_type')
    my_ggsave(curr_sobj@misc$plots$umap_pat_split, filename = paste0(sample_fig_dir, sprintf("/%s_umap_pat_split.png", opt$regulon)))
    
    # DimPlot(curr_sobj, reduction = "umap", group.by = 'pat_type', cols=c(colors$green, colors$yellow, colors$red)) + theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7)) + labs(title='T cells')
    
    ###################### TF MARKERS
    ####### BY PAT_ID
    
    # Changing the identities to pat_id
    cluster_idents <- Idents(object = curr_sobj)
    if (dim(curr_sobj)[2] > 10000){
      plan('multiprocess', workers=46)
    }
    # curr_sobj$dis_type <- ifelse(curr_sobj$pat_type == 'C', 'Healthy', 'Covid')
    # curr_sobj$new_type <- ifelse(curr_sobj$pat_type == 'S', 'S', 'NS')
    # sobj.markers <- FindAllMarkers(curr_sobj, min.pct = 0, logfc.threshold = 0, return.thresh=1.01, verbose = F)
    # write.table(sobj.markers, file=sprintf('tmp/tf_markers_df_%s.tsv', curr_type), sep='\t')
    # C_markers <- sobj.markers[sobj.markers$cluster == 'C',]
    # M_markers <- sobj.markers[sobj.markers$cluster == 'M',]
    # S_markers <- sobj.markers[sobj.markers$cluster == 'S',]
    # C_best <- C_markers[(C_markers$p_val_adj < 0.05) & (C_markers$avg_log2FC > 1),]
    # M_best <- M_markers[(M_markers$p_val_adj < 0.05) & (M_markers$avg_log2FC > 1),]
    # S_best <- S_markers[(S_markers$p_val_adj < 0.05) & (S_markers$avg_log2FC > 1),]
    # cat('C:\n')
    # cat(C_best$gene, sep='\n')
    # cat('M:\n')
    # cat(M_best$gene, sep='\n')
    # cat('S:\n')
    # cat(S_best$gene, sep='\n')
    
    # source('src/func.R')
    # markers_df <- read.table(file='tmp/tf_markers_df.tsv', sep='\t')
    # out_ora <- run_ora(markers_df, p_value_cutoff=1)
    # out_gsea <- run_gsea(markers_df)
    
    
    Idents(object = curr_sobj) <- curr_sobj@meta.data$pat_id
    
    # find markers for every patient compared to all remaining cells, report only the positive ones
    sobj.markers <- FindAllMarkers(curr_sobj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = F)
    # XXX
    # curr_sobj.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
    sobj.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
    curr_sobj@misc$markers_pat_id <- sobj.markers
    # XXX
    # top10 <- sobj.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
    top10 <- sobj.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) %>% arrange(factor(cluster, levels = c('C', 'M', 'S')))
    curr_sobj@misc$plots$markers_pat_id <- DoHeatmap(subset(curr_sobj, downsample = 1000), features = top10$gene) + NoLegend()
    my_ggsave(curr_sobj@misc$plots$markers_pat_id, filename = paste0(sample_fig_dir, sprintf("/%s_markers_pat_id.png", opt$regulon)))
    
    # DoHeatmap(subset(curr_sobj, downsample = 1000), features = top10$gene, group.by='pat_id', angle=45, size=3) + NoLegend() + theme(axis.text.y = element_text(size = 6))
    # DoMultiBarHeatmap(subset(curr_sobj, downsample = 1000), features = top10$gene, group.by='pat_type', angle=45, size=3, additional.group.by='pat_id', label=F) + NoLegend() + theme(axis.text.y = element_text(size = 6))
    
    ####### BY PAT_TYPE
    
    # Changing the identities to pat type  m
    Idents(object = curr_sobj) <- curr_sobj@meta.data$pat_type
    
    # find markers for every patient group compared to all remaining cells, report only the positive ones
    sobj.markers <- FindAllMarkers(curr_sobj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = F)
    # XXX
    # sobj.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
    sobj.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
    curr_sobj@misc$markers_pat_type <- sobj.markers
    # XXX
    # top10 <- sobj.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
    top10 <- sobj.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) %>% arrange(factor(cluster, levels = c('C', 'M', 'S')))
    curr_sobj@misc$plots$markers_pat_type <- DoHeatmap(subset(curr_sobj, downsample = 1000), features = top10$gene) + NoLegend()
    my_ggsave(curr_sobj@misc$plots$markers_pat_type, filename = paste0(sample_fig_dir, sprintf("/%s_markers_pat_type.png", opt$regulon)))
    
    # out <- DoHeatmap(subset(curr_sobj, downsample = 1000), features = top10$gene, group.by='pat_type', group.colors=c(colors$green, colors$yellow, colors$red), angle=0, label=F) + NoLegend()
    # ggsave(out, filename=('tmp/tf_activity_pat_type.pdf'), width=8, height=5)
    
    #### Compare viper scores with TF expression
    # load(file=paste0(sample_data_dir, "/seurat_object.RData"))
    # Idents(object = curr_sobj) <- curr_sobj@meta.data$pat_type
    # DoHeatmap(subset(curr_sobj, downsample = 1000), features = top10$gene, group.by='pat_type', group.colors=c(colors$green, colors$yellow, colors$red), angle=0) + NoLegend()
    
    # Setting cluster identities back
    Idents(object = curr_sobj) <- cluster_idents
    plan('sequential')
    
    ###################### FINISH
    cat("      - Saving\n")
    
    curr_sobj@misc$plots <- NULL
    curr_sobj@misc$markers_pat_type <- NULL
    curr_sobj@misc$markers_pat_id <- NULL
    
    # Reload previously saved data
    # load(file=paste0(sample_data_dir, sprintf("/%s_seurat_object.RData", opt$regulon)))  # load(file=paste0(sample_data_dir, "/viper_seurat_object.RData"))
    
    # TODO: WHY FAILS? (JCB)
    if(opt$serialize==T) save(curr_sobj, file=paste0(sample_data_dir, sprintf("/%s_seurat_object.RData", opt$regulon)))
    
    write.table(curr_sobj@assays[[opt$regulon]]@data, file=paste0(sample_data_dir, sprintf("/%s_raw_data.tsv", opt$regulon)), sep="\t", row.names=T, col.names=colnames(curr_sobj@assays[[opt$regulon]]@data), quote=F)
    write.table(curr_sobj@assays[[opt$regulon]]@scale.data, file=paste0(sample_data_dir, sprintf("/%s_scaled_data.tsv", opt$regulon)), sep="\t", row.names=T, col.names=colnames(curr_sobj@assays[[opt$regulon]]@scale.data), quote=F)
    
    # patient type merged info
    for (k in 1:3){
      type <- pat_type_levels[k]
      tryCatch({
        curr_subsobj <- subset(curr_sobj, subset = pat_type == type)
        
        if(opt$serialize==T) save(curr_subsobj, file=paste0(sample_data_dir, sprintf("/%s_seurat_object_%s_type.RData", opt$regulon, type)))
        write.table(curr_subsobj@assays[[opt$regulon]]@data, file=paste0(sample_data_dir, sprintf("/%s_raw_data_%s_type.tsv", opt$regulon, type)), sep="\t", row.names=T, col.names=colnames(curr_subsobj@assays[[opt$regulon]]@data), quote=F)
      
      },
      error = function(e) {
        cat(sprintf('      - No "%s" cells from patient type: "%s.."\n', curr_type, type))
      }
      )
    }
    
    rm('curr_sobj')
    gc()
  }, 
  
  error = function(e) {
    cat(paste0("Encountered error: '", e, "' when processing cell type: ", curr_type, '\n'))
  }
  )
}

cat("\n\n[finished]\n\n")
