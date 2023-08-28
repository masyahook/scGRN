############## Production of thesis figures

###### Data description section

### UMAP visualization of data cell distribution (merge_by_cell_type.R script)
sobj_combined <- RunUMAP(object = sobj_combined, dims=1:20, verbose = F)
p1 <- DimPlot(sobj_combined, reduction = "umap", group.by = "pat_id") + labs(title="")
p2 <- DimPlot(sobj_combined, reduction = "umap", group.by = "pat_type", cols=c(colors$green, colors$yellow, colors$red)) + 
  labs(title="")
res_p <- p1 + p2
res_p
ggsave(res_p, filename = paste0('figs', "/umap_all.pdf"), width = 8, height = 4)

### UMAP visualization of different cell types (merge_by_cell_type.R script)

# T_cells + Macrophages
i <- 1
curr_type <- colnames(meta)[i]
cat("  > Processing cell type: ",curr_type," (",i," of ", ncol(meta),")\n", sep="")
sample_dir <- paste0(cell_type_dir, "/", curr_type)
dir.create(sample_dir, showWarnings = F)
sample_fig_dir <- paste0(sample_dir, '/figs/Seurat')
dir.create(sample_fig_dir, showWarnings = F, recursive=T)
sample_data_dir <- paste0(sample_dir, '/data/Seurat')
dir.create(sample_data_dir, showWarnings = F, recursive=T)
load(file=paste0(sample_data_dir, "/seurat_object.RData"))

p1 <- DimPlot(curr_sobj, reduction = "umap", split.by = "pat_id", group.by = 'pat_type', 
              cols=c(colors$green, colors$yellow, colors$red)) + 
  theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7)) +
  labs(title="", x="")

i <- 2
curr_type <- colnames(meta)[i]
cat("  > Processing cell type: ",curr_type," (",i," of ", ncol(meta),")\n", sep="")
sample_dir <- paste0(cell_type_dir, "/", curr_type)
dir.create(sample_dir, showWarnings = F)
sample_fig_dir <- paste0(sample_dir, '/figs/Seurat')
dir.create(sample_fig_dir, showWarnings = F, recursive=T)
sample_data_dir <- paste0(sample_dir, '/data/Seurat')
dir.create(sample_data_dir, showWarnings = F, recursive=T)
load(file=paste0(sample_data_dir, "/seurat_object.RData"))

p2 <- DimPlot(curr_sobj, reduction = "umap", split.by = "pat_id", group.by = 'pat_type', 
              cols=c(colors$green, colors$yellow, colors$red)) + 
  theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7)) +
  labs(title="")
res_p <- p1 + p2 + theme(panel.margin=unit(0, "cm"))
res_p
ggsave(res_p, filename = paste0('figs', sprintf("/umap_T_cells_Macrophage.pdf", curr_type)), width = 12, height = 5)

# NK_cell + Monocyte
i <- 6
curr_type <- colnames(meta)[i]
cat("  > Processing cell type: ",curr_type," (",i," of ", ncol(meta),")\n", sep="")
sample_dir <- paste0(cell_type_dir, "/", curr_type)
dir.create(sample_dir, showWarnings = F)
sample_fig_dir <- paste0(sample_dir, '/figs/Seurat')
dir.create(sample_fig_dir, showWarnings = F, recursive=T)
sample_data_dir <- paste0(sample_dir, '/data/Seurat')
dir.create(sample_data_dir, showWarnings = F, recursive=T)
load(file=paste0(sample_data_dir, "/seurat_object.RData"))

p1 <- DimPlot(curr_sobj, reduction = "umap", split.by = "pat_id", group.by = 'pat_type', 
              cols=c(colors$green, colors$yellow, colors$red)) + 
  theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7)) +
  labs(title="", x="")

i <- 7
curr_type <- colnames(meta)[i]
cat("  > Processing cell type: ",curr_type," (",i," of ", ncol(meta),")\n", sep="")
sample_dir <- paste0(cell_type_dir, "/", curr_type)
dir.create(sample_dir, showWarnings = F)
sample_fig_dir <- paste0(sample_dir, '/figs/Seurat')
dir.create(sample_fig_dir, showWarnings = F, recursive=T)
sample_data_dir <- paste0(sample_dir, '/data/Seurat')
dir.create(sample_data_dir, showWarnings = F, recursive=T)
load(file=paste0(sample_data_dir, "/seurat_object.RData"))

p2 <- DimPlot(curr_sobj, reduction = "umap", split.by = "pat_id", group.by = 'pat_type', 
              cols=c(colors$green, colors$yellow, colors$red)) + 
  theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7)) +
  labs(title="")
res_p <- p1 + p2 + theme(panel.margin=unit(0, "cm"))
res_p
ggsave(res_p, filename = paste0('figs', sprintf("/umap_NK_cell_Monocyte.pdf", curr_type)), width = 12, height = 5)

### UMAP visualization of the cell types

library(gridExtra)
i <- 1
curr_type <- colnames(meta)[i]

# Specifying data/fig folders
sample_dir <- paste0(cell_type_dir, "/", curr_type)
sample_fig_dir <- paste0(sample_dir, '/figs/Seurat')
sample_data_dir <- paste0(sample_dir, '/data/Seurat')

load(file=paste0(sample_data_dir, "/seurat_object.RData"))
p1 <- DimPlot(curr_sobj, reduction = "umap", group.by = 'pat_type', cols=c(colors$green, colors$yellow, colors$red)) + theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7)) + labs(title='T cells')

i <- 2
curr_type <- colnames(meta)[i]

# Specifying data/fig folders
sample_dir <- paste0(cell_type_dir, "/", curr_type)
sample_fig_dir <- paste0(sample_dir, '/figs/Seurat')
sample_data_dir <- paste0(sample_dir, '/data/Seurat')

load(file=paste0(sample_data_dir, "/seurat_object.RData"))
p2 <- DimPlot(curr_sobj, reduction = "umap", group.by = 'pat_type', cols=c(colors$green, colors$yellow, colors$red)) + theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7)) + labs(title='Macrophages')

i <- 6
curr_type <- colnames(meta)[i]

# Specifying data/fig folders
sample_dir <- paste0(cell_type_dir, "/", curr_type)
sample_fig_dir <- paste0(sample_dir, '/figs/Seurat')
sample_data_dir <- paste0(sample_dir, '/data/Seurat')

load(file=paste0(sample_data_dir, "/seurat_object.RData"))
p3 <- DimPlot(curr_sobj, reduction = "umap", group.by = 'pat_type', cols=c(colors$green, colors$yellow, colors$red)) + theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7)) + labs(title='NK cells')

i <- 7
curr_type <- colnames(meta)[i]

# Specifying data/fig folders
sample_dir <- paste0(cell_type_dir, "/", curr_type)
sample_fig_dir <- paste0(sample_dir, '/figs/Seurat')
sample_data_dir <- paste0(sample_dir, '/data/Seurat')

load(file=paste0(sample_data_dir, "/seurat_object.RData"))
p4 <- DimPlot(curr_sobj, reduction = "umap", group.by = 'pat_type', cols=c(colors$green, colors$yellow, colors$red)) + theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7)) + labs(title='Monocytes')

p <- grid.arrange(p1, p2, p3, p4, nrow=2, ncol=2)

ggsave(p, filename='tmp/gene_expr_umaps.pdf', width=8, height=6)

### VIPER UMAP visualization of the scores

library(gridExtra)
i <- 1
curr_type <- colnames(meta)[i]

# Specifying data/fig folders
sample_dir <- paste0(cell_type_dir, "/", curr_type)
sample_fig_dir <- paste0(sample_dir, '/figs/Seurat')
sample_data_dir <- paste0(sample_dir, '/data/Seurat')

load(file=paste0(sample_data_dir, sprintf("/%s_seurat_object.RData", opt$regulon)))

DefaultAssay(object = curr_sobj) <- opt$regulon
p1 <- DimPlot(curr_sobj, reduction = "umap", group.by = 'pat_type', cols=c(colors$green, colors$yellow, colors$red)) + theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7)) + labs(title='T cells')

i <- 2
curr_type <- colnames(meta)[i]

# Specifying data/fig folders
sample_dir <- paste0(cell_type_dir, "/", curr_type)
sample_fig_dir <- paste0(sample_dir, '/figs/Seurat')
sample_data_dir <- paste0(sample_dir, '/data/Seurat')

load(file=paste0(sample_data_dir, sprintf("/%s_seurat_object.RData", opt$regulon)))

DefaultAssay(object = curr_sobj) <- opt$regulon
p2 <- DimPlot(curr_sobj, reduction = "umap", group.by = 'pat_type', cols=c(colors$green, colors$yellow, colors$red)) + theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7)) + labs(title='Macrophages')

i <- 6
curr_type <- colnames(meta)[i]

# Specifying data/fig folders
sample_dir <- paste0(cell_type_dir, "/", curr_type)
sample_fig_dir <- paste0(sample_dir, '/figs/Seurat')
sample_data_dir <- paste0(sample_dir, '/data/Seurat')

load(file=paste0(sample_data_dir, sprintf("/%s_seurat_object.RData", opt$regulon)))

DefaultAssay(object = curr_sobj) <- opt$regulon
p3 <- DimPlot(curr_sobj, reduction = "umap", group.by = 'pat_type', cols=c(colors$green, colors$yellow, colors$red)) + theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7)) + labs(title='NK cells')

i <- 7
curr_type <- colnames(meta)[i]

# Specifying data/fig folders
sample_dir <- paste0(cell_type_dir, "/", curr_type)
sample_fig_dir <- paste0(sample_dir, '/figs/Seurat')
sample_data_dir <- paste0(sample_dir, '/data/Seurat')

load(file=paste0(sample_data_dir, sprintf("/%s_seurat_object.RData", opt$regulon)))

DefaultAssay(object = curr_sobj) <- opt$regulon
p4 <- DimPlot(curr_sobj, reduction = "umap", group.by = 'pat_type', cols=c(colors$green, colors$yellow, colors$red)) + theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7)) + labs(title='Monocytes')

p <- grid.arrange(p1, p2, p3, p4, nrow=2, ncol=2)

ggsave(p, filename='tmp/tf_umaps.pdf', width=8, height=6)

### VIPER Differential activity analysis (viper_processing_cell_type.R script)

# T_cells
i <- 1
curr_type <- colnames(meta)[i]

# Specifying data/fig folders
sample_dir <- paste0(cell_type_dir, "/", curr_type)
sample_fig_dir <- paste0(sample_dir, '/figs/Seurat')
sample_data_dir <- paste0(sample_dir, '/data/Seurat')

load(file=paste0(sample_data_dir, sprintf("/%s_seurat_object.RData", opt$regulon)))

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

out <- DoHeatmap(subset(curr_sobj, downsample = 1000), features = top10$gene, group.by='pat_type', group.colors=c(colors$green, colors$yellow, colors$red), angle=0, label=F) + NoLegend()
ggsave(out, filename=sprintf('tmp/tf_activity_pat_type_%s.pdf', curr_type), width=8, height=5)

sobj.markers <- FindAllMarkers(curr_sobj, min.pct = 0, logfc.threshold = 0, return.thresh=1.01, verbose = F)
write.table(sobj.markers, file=sprintf('tmp/tf_markers_df_%s.tsv', curr_type), sep='\t')

sobj.markers <- FindMarkers(curr_sobj, ident.1 = 'M', ident.2 = 'C', min.pct = 0, 
                            logfc.threshold = 0, return.thresh=1.01, verbose = F)
write.table(sobj.markers, file=sprintf('tmp/tf_markers_df_M_C_%s.tsv', curr_type), sep='\t')

sobj.markers <- FindMarkers(curr_sobj, ident.1 = 'S', ident.2 = 'C', min.pct = 0, 
                            logfc.threshold = 0, return.thresh=1.01, verbose = F)
write.table(sobj.markers, file=sprintf('tmp/tf_markers_df_S_C_%s.tsv', curr_type), sep='\t')

sobj.markers <- FindMarkers(curr_sobj, ident.1 = 'S', ident.2 = 'M', min.pct = 0, 
                            logfc.threshold = 0, return.thresh=1.01, verbose = F)
write.table(sobj.markers, file=sprintf('tmp/tf_markers_df_S_M_%s.tsv', curr_type), sep='\t')

# source('src/func.R')
# markers_df <- read.table(file=sprintf('tmp/tf_markers_df_%s.tsv', curr_type), sep='\t')
# out_ora <- run_ora(markers_df, p_value_cutoff=1)
# out_gsea <- run_gsea(markers_df)

# Macrophage
i <- 2
curr_type <- colnames(meta)[i]

# Specifying data/fig folders
sample_dir <- paste0(cell_type_dir, "/", curr_type)
sample_fig_dir <- paste0(sample_dir, '/figs/Seurat')
sample_data_dir <- paste0(sample_dir, '/data/Seurat')

load(file=paste0(sample_data_dir, sprintf("/%s_seurat_object.RData", opt$regulon)))

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

out <- DoHeatmap(subset(curr_sobj, downsample = 1000), features = top10$gene, group.by='pat_type', group.colors=c(colors$green, colors$yellow, colors$red), angle=0, label=F) + NoLegend()
ggsave(out, filename=sprintf('tmp/tf_activity_pat_type_%s.pdf', curr_type), width=8, height=5)

sobj.markers <- FindAllMarkers(curr_sobj, min.pct = 0, logfc.threshold = 0, return.thresh=1.01, verbose = F)
write.table(sobj.markers, file=sprintf('tmp/tf_markers_df_%s.tsv', curr_type), sep='\t')

sobj.markers <- FindMarkers(curr_sobj, ident.1 = 'M', ident.2 = 'C', min.pct = 0, 
                            logfc.threshold = 0, return.thresh=1.01, verbose = F)
write.table(sobj.markers, file=sprintf('tmp/tf_markers_df_M_C_%s.tsv', curr_type), sep='\t')

sobj.markers <- FindMarkers(curr_sobj, ident.1 = 'S', ident.2 = 'C', min.pct = 0, 
                            logfc.threshold = 0, return.thresh=1.01, verbose = F)
write.table(sobj.markers, file=sprintf('tmp/tf_markers_df_S_C_%s.tsv', curr_type), sep='\t')

sobj.markers <- FindMarkers(curr_sobj, ident.1 = 'S', ident.2 = 'M', min.pct = 0, 
                            logfc.threshold = 0, return.thresh=1.01, verbose = F)
write.table(sobj.markers, file=sprintf('tmp/tf_markers_df_S_M_%s.tsv', curr_type), sep='\t')

# source('src/func.R')
# markers_df <- read.table(file=sprintf('tmp/tf_markers_df_%s.tsv', curr_type), sep='\t')
# out_ora <- run_ora(markers_df, p_value_cutoff=1)
# out_gsea <- run_gsea(markers_df)

# NK_cell
i <- 6
curr_type <- colnames(meta)[i]

# Specifying data/fig folders
sample_dir <- paste0(cell_type_dir, "/", curr_type)
sample_fig_dir <- paste0(sample_dir, '/figs/Seurat')
sample_data_dir <- paste0(sample_dir, '/data/Seurat')

load(file=paste0(sample_data_dir, sprintf("/%s_seurat_object.RData", opt$regulon)))

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

out <- DoHeatmap(subset(curr_sobj, downsample = 1000), features = top10$gene, group.by='pat_type', group.colors=c(colors$green, colors$yellow, colors$red), angle=0, label=F) + NoLegend()
ggsave(out, filename=sprintf('tmp/tf_activity_pat_type_%s.pdf', curr_type), width=8, height=5)

sobj.markers <- FindAllMarkers(curr_sobj, min.pct = 0, logfc.threshold = 0, return.thresh=1.01, verbose = F)
write.table(sobj.markers, file=sprintf('tmp/tf_markers_df_%s.tsv', curr_type), sep='\t')

sobj.markers <- FindMarkers(curr_sobj, ident.1 = 'M', ident.2 = 'C', min.pct = 0, 
                            logfc.threshold = 0, return.thresh=1.01, verbose = F)
write.table(sobj.markers, file=sprintf('tmp/tf_markers_df_M_C_%s.tsv', curr_type), sep='\t')

sobj.markers <- FindMarkers(curr_sobj, ident.1 = 'S', ident.2 = 'C', min.pct = 0, 
                            logfc.threshold = 0, return.thresh=1.01, verbose = F)
write.table(sobj.markers, file=sprintf('tmp/tf_markers_df_S_C_%s.tsv', curr_type), sep='\t')

sobj.markers <- FindMarkers(curr_sobj, ident.1 = 'S', ident.2 = 'M', min.pct = 0, 
                            logfc.threshold = 0, return.thresh=1.01, verbose = F)
write.table(sobj.markers, file=sprintf('tmp/tf_markers_df_S_M_%s.tsv', curr_type), sep='\t')

# source('src/func.R')
# markers_df <- read.table(file=sprintf('tmp/tf_markers_df_%s.tsv', curr_type), sep='\t')
# out_ora <- run_ora(markers_df, p_value_cutoff=1)
out_gsea <- run_gsea(markers_df)

# T_cells
i <- 7
curr_type <- colnames(meta)[i]

# Specifying data/fig folders
sample_dir <- paste0(cell_type_dir, "/", curr_type)
sample_fig_dir <- paste0(sample_dir, '/figs/Seurat')
sample_data_dir <- paste0(sample_dir, '/data/Seurat')

load(file=paste0(sample_data_dir, sprintf("/%s_seurat_object.RData", opt$regulon)))

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

out <- DoHeatmap(subset(curr_sobj, downsample = 1000), features = top10$gene, group.by='pat_type', group.colors=c(colors$green, colors$yellow, colors$red), angle=0, label=F) + NoLegend()
ggsave(out, filename=sprintf('tmp/tf_activity_pat_type_%s.pdf', curr_type), width=8, height=5)

sobj.markers <- FindAllMarkers(curr_sobj, min.pct = 0, logfc.threshold = 0, return.thresh=1.01, verbose = F)
write.table(sobj.markers, file=sprintf('tmp/tf_markers_df_%s.tsv', curr_type), sep='\t')

sobj.markers <- FindMarkers(curr_sobj, ident.1 = 'M', ident.2 = 'C', min.pct = 0, 
                            logfc.threshold = 0, return.thresh=1.01, verbose = F)
write.table(sobj.markers, file=sprintf('tmp/tf_markers_df_M_C_%s.tsv', curr_type), sep='\t')

sobj.markers <- FindMarkers(curr_sobj, ident.1 = 'S', ident.2 = 'C', min.pct = 0, 
                            logfc.threshold = 0, return.thresh=1.01, verbose = F)
write.table(sobj.markers, file=sprintf('tmp/tf_markers_df_S_C_%s.tsv', curr_type), sep='\t')

sobj.markers <- FindMarkers(curr_sobj, ident.1 = 'S', ident.2 = 'M', min.pct = 0, 
                            logfc.threshold = 0, return.thresh=1.01, verbose = F)
write.table(sobj.markers, file=sprintf('tmp/tf_markers_df_S_M_%s.tsv', curr_type), sep='\t')

# source('src/func.R')
# markers_df <- read.table(file=sprintf('tmp/tf_markers_df_%s.tsv', curr_type), sep='\t')
# out_ora <- run_ora(markers_df, p_value_cutoff=1)
# out_gsea <- run_gsea(markers_df)
