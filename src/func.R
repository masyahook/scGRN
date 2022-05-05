suppressPackageStartupMessages({
  #if (!require("pacman")) install.packages("pacman")
  #list.of.packages <- c("BiocManager","dplyr","SingleR","Matrix","Seurat","future","pheatmap","ggplot2","optparse","hdf5r")
  #BiocManager::install("SingleR")
  #BiocManager::install('limma')
  #BiocManager::install('SingleCellExperiment')
  #pacman::p_load(list.of.packages, character.only = TRUE)
  library(SingleR)
  library(Matrix)
  library(Seurat)
  library(dplyr)
  library(future)
  library(pheatmap)
  library(ggplot2)
  library(optparse)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

enrich_list <- function(x_list){
  gene.df <- bitr(names(x_list), fromType = "SYMBOL",
                  toType = c("ENTREZID"), 
                  OrgDb = org.Hs.eg.db)
  entrez_x <- x_list
  names(entrez_x) <- gene.df$ENTREZID
  
  Sys.setenv("/scratch/tmp")
  
  gse <- gseGO(geneList     = entrez_x, 
               OrgDb        = org.Hs.eg.db, 
               ont          = "BP", 
               minGSSize    = 100, 
               maxGSSize    = 500, 
               pvalueCutoff = 0.05, 
               verbose      = FALSE)
  
  return(gse)
}

a <- read.table('tmp.tsv', sep='\t')
x <- c()

for (i in 2:dim(a)[1]){
  key <- a[i, 'V1']
  val <- a[i, 'V2']
  x[key] = val
}

res <- enrich_list(x)
