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
  library(ReactomePA)
  library(org.Hs.eg.db)
})

compute_ranking <- function(df, rank_type='FC'){
  if (rank_type == 'FC'){
     out <- df['avg_log2FC']
  } else if (rank_type == 'signed_p'){
    out <- df['avg_log2FC'] * (-log10(df['p_val_adj']))
    index <- is.infinite(out[['avg_log2FC']])
    out[index, 'avg_log2FC'] <- max(out[['avg_log2FC']][is.finite(out[['avg_log2FC']])], na.rm=T)
  } else {
    stop(sprintf("rank_type '%s' is incorrect, please change..", rank_type), call.=FALSE)
  }
  
  return(out)
}

gsea_cP <- function(ranked_list, db='ALL', from_mat=T, rank_type='signed_p', 
                        p_value_cutoff=1, ont='BP'){
  
  # Obtaining the ranked list from DE results matrix
  if (from_mat == T){
    ranked_list <- compute_ranking(ranked_list, rank_type=rank_type)
  }
  
  # Transforming to vector format
  if (class(ranked_list) == 'data.frame'){
    x <- unlist(ranked_list)
    names(x) <- rownames(ranked_list)
  } else {
    x <- ranked_list
  }
  
  # Sorting the rankings
  x <- sort(x, decreasing=T)
  
  # Converting gene symbol to ENTREZ ID
  gene.df <- bitr(names(x), fromType = "SYMBOL",
                  toType = c("ENTREZID"), 
                  OrgDb = org.Hs.eg.db)
  entrez_x <- x
  names(entrez_x) <- gene.df$ENTREZID
  
  # Changing temporary directory in RStudio
  Sys.setenv(TMPDIR="/scratch/tmp")
  
  # Running GSEA
  if (db == 'GO'){
    
    # GO annotation
    out <- gseGO(geneList     = entrez_x, 
                 OrgDb        = org.Hs.eg.db, 
                 ont          = ont, 
                 minGSSize    = 100, 
                 maxGSSize    = 500, 
                 pvalueCutoff = p_value_cutoff, 
                 verbose      = FALSE)
    
  } else if (db == 'KEGG'){
    
    # KEGG annotation
    out <- gseKEGG(geneList     = entrez_x,
                   organism     = 'hsa',
                   minGSSize    = 120,
                   pvalueCutoff = p_value_cutoff,
                   verbose      = FALSE)
    
  } else if (db == 'WP'){
    
    # WikiPathways annotation
    out <- gseWP(entrez_x, organism = "Homo sapiens")
    
  } else if (db == 'Reactome'){
    
    # Reactome annotation
    out <- gsePathway(entrez_x, 
                      pvalueCutoff = p_value_cutoff,
                      pAdjustMethod = "BH", 
                      verbose = FALSE)
  } else if (db == 'ALL'){
    
    # Using all annotations
    out <- list()
    tryCatch({
      cat(paste0("Running GSEA on GO data..", '\n'))
      out['GO'] <- gseGO(geneList     = entrez_x, 
                         OrgDb        = org.Hs.eg.db, 
                         ont          = ont, 
                         minGSSize    = 100, 
                         maxGSSize    = 500, 
                         pvalueCutoff = p_value_cutoff, 
                         verbose      = FALSE)
      cat(paste0("    Success!", '\n'))
    }, 
    
    error = function(e) {
      cat(paste0("    Encountered error: '", e, "' when getting GO annotations", '\n'))
    }
    )
    
    tryCatch({
      cat(paste0("Running GSEA on KEGG data..", '\n'))
      out['KEGG'] <- gseKEGG(geneList     = entrez_x,
                           organism     = 'hsa',
                           minGSSize    = 120,
                           pvalueCutoff = p_value_cutoff,
                           verbose      = FALSE)
      cat(paste0("    Success!", '\n'))
    }, 
    
    error = function(e) {
      cat(paste0("    Encountered error: '", e, "' when getting KEGG annotations", '\n'))
    }
    )
    
    tryCatch({
      cat(paste0("Running GSEA on WikiPathways data..", '\n'))
      out['WP'] <- gseWP(entrez_x, organism = "Homo sapiens")
      cat(paste0("    Success!", '\n'))
    }, 
    
    error = function(e) {
      cat(paste0("    Encountered error: '", e, "' when getting WikiPathways annotations", '\n'))
    }
    )
    
    tryCatch({
      cat(paste0("Running GSEA on Reactome data..", '\n'))
      out['Reactome'] <- gsePathway(entrez_x, 
                                    pvalueCutoff = p_value_cutoff,
                                    pAdjustMethod = "BH", 
                                    verbose = FALSE)
      cat(paste0("    Success!", '\n'))
    }, 
    
    error = function(e) {
      cat(paste0("    Encountered error: '", e, "' when getting Reactome annotations", '\n'))
    }
    )
    
  } else {
    stop(sprintf("db '%s' is incorrect, please change..", rank_type), call.=FALSE)
  }
  
  return(out)
}

enrich_DE_results <- function(markers, top_n=10, db_run='ALL', from_mat=T, 
                              rank_type='signed_p', p_value_cutoff=1, ont='BP'){
  
  # Setting some params
  selected_cols <- c("ID", "Description", "p.adjust")
  dbs <- c('GO', 'KEGG', 'Reactome')
  groups <- unique(markers$cluster)
  
  # Running GSEA on each group and saving all results
  all_res <- list()
  cat('\n\nProcessing..\n\n\n')
  for (group in groups){
    cat(paste0(sprintf("Finding enrichment terms for group: '%s'", group), '\n'))
    group_markers <- markers[markers$cluster == group,]
    all_res[[group]] <- gsea_cP(group_markers, db=db_run, from_mat=from_mat, 
                              rank_type=rank_type,
                              p_value_cutoff=p_value_cutoff, ont=ont)
    cat('\n\n')
  }
  
  # # Saving shortened version of results
  short_res <- list()
  cat('\n\nGetting shortened results..\n\n')
  for (db in dbs){

    curr_short_res <- data.frame(row.names=1:top_n)
    for (group in groups){

      # Getting only important info
      curr <- head(all_res[[group]][[db]]@result[,selected_cols], top_n)
      print(curr)
      colnames(curr) <- lapply(colnames(curr),
                               function(x) sprintf('%s_%s', x, group))

      # Filling NA the rows if results for current group are insufficient
      if (nrow(curr) < top_n){
        tmp_n <- nrow(curr) + 1
        curr[tmp_n:top_n,] <- NA
      }
      rownames(curr) <- 1:top_n

      curr_short_res <- cbind(curr_short_res, curr)
    }

    short_res[[db]] <- curr_short_res
  }

  cat('Done!\n')

  out <- list(short=short_res, all=all_res)
  
  return(out)
}

# a <- read.table('tmp.tsv', sep='\t')
# x <- c()
# 
# for (i in 2:dim(a)[1]){
#   key <- a[i, 'V1']
#   val <- a[i, 'V2']
#   x[key] = val
# }
# 
# res <- enrich_list(x)
