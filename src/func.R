suppressPackageStartupMessages({
  #if (!require("pacman")) install.packages("pacman")
  #list.of.packages <- c("BiocManager","dplyr","SingleR","Matrix","Seurat","future","pheatmap","ggplot2","optparse","hdf5r")
  #BiocManager::install("SingleR")
  #BiocManager::install('limma')
  #BiocManager::install('SingleCellExperiment')
  #pacman::p_load(list.of.packages, character.only = TRUE)
  library(ggplot2)
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
    rownames(ranked_list) <- ranked_list$gene
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
  entrez_x <- entrez_x[!is.na(names(entrez_x))]
  
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
    out <- gseWP(entrez_x, organism = "Homo sapiens", 
                 pvalueCutoff = p_value_cutoff,
                 verbose = FALSE)
    
  } else if (db == 'Reactome'){
    
    # Reactome annotation
    out <- gsePathway(entrez_x, minGSSize = 100,
                      pvalueCutoff = p_value_cutoff,
                      pAdjustMethod = "BH", 
                      verbose = FALSE)
  } else if (db == 'ALL'){
    
    # Using all annotations
    out <- list()
    tryCatch({
      cat(paste0("Running GSEA on GO data..", '\n'))
      out[['GO']] <- gseGO(geneList     = entrez_x, 
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
      out[['KEGG']] <- gseKEGG(geneList     = entrez_x,
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
      out[['WP']] <- gseWP(entrez_x, organism = "Homo sapiens",
                           minGSSize = 100,
                           pvalueCutoff = p_value_cutoff,
                           verbose = FALSE)
      cat(paste0("    Success!", '\n'))
    }, 
    
    error = function(e) {
      cat(paste0("    Encountered error: '", e, "' when getting WikiPathways annotations", '\n'))
    }
    )
    
    tryCatch({
      cat(paste0("Running GSEA on Reactome data..", '\n'))
      out[['Reactome']] <- gsePathway(entrez_x, minGSSize = 100,
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

gsea_from_DE_results <- function(markers, top_n=10, db_run='ALL', from_mat=T, 
                                 rank_type='signed_p', p_value_cutoff=1, ont='BP'){
  
  # Setting some params
  selected_cols <- c("ID", "Description", "p.adjust")
  dbs <- c('GO', 'KEGG', 'WP', 'Reactome')
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
  
  # Saving shortened version of results
  short_res <- list()
  cat('\nGetting shortened results..\n\n')
  for (db in dbs){

    curr_short_res <- data.frame(row.names=1:top_n)
    for (group in groups){

      # Getting only important info
      curr <- head(all_res[[group]][[db]]@result[,selected_cols], top_n)
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
  saveRDS(out, 'tmp/tmp_gsea_from_DE_results.rds')
  
  return(out)
}

ora_from_DE_results <- function(markers_df, top_n=10, top_logFC = 1, 
                                top_p_val = 0.05, p_value_cutoff=0.05, 
                                ont='BP'){
  
  # Setting some params
  pat_type_levels <- c('C', 'M', 'S')
  pat_levels <- c('C51', 'C52', 'C100', 'C141', 'C142', 'C144', 'C143', 'C145', 'C146', 'C148', 'C149', 'C152')
  groups <- ifelse((intersect(pat_type_levels, unique(markers_df$cluster)) > 0), 
                   pat_type_levels, pat_levels)
  selected_cols <- c("ID", "Description", "p.adjust")
  dbs <- c('GO', 'KEGG', 'WP', 'Reactome')
  final_markers_df <- markers_df[(markers_df$p_val_adj < top_p_val) & 
                             (markers_df$avg_log2FC > top_logFC),]
  
  tmp_universe <- unique(markers_df$gene)
  universe <- bitr(tmp_universe, fromType = 'SYMBOL', toType = c('ENTREZID'),
                  OrgDb = org.Hs.eg.db)$ENTREZID
  
  # Getting gene marker lists
  markers <- list()
  for (group in groups){
    tmp_markers <- final_markers_df[final_markers_df$cluster == group,]
    tmp_df <- bitr(tmp_markers$gene, fromType = "SYMBOL",
                    toType = c("ENTREZID"), 
                    OrgDb = org.Hs.eg.db)
    markers[[group]] <- tmp_df$ENTREZID
  }
  
  cat('\n\nProcessing..\n\n\n')
  out <- list()
  for (db in dbs){
    cat(paste0(sprintf("Finding enrichment terms for db: '%s'", db), '\n'))
    
    if (db == 'GO'){
      cmd <- sprintf("compareCluster(geneClusters = markers, fun = enrichGO, 
                             universe = universe, OrgDb = org.Hs.eg.db, 
                             ont = ont, pvalueCutoff = p_value_cutoff,
                             qvalueCutoff = p_value_cutoff
                             )")
    } else if (db == 'KEGG'){
      cmd <- sprintf("compareCluster(geneClusters = markers, fun = enrichKEGG,
                             organism     = 'hsa', universe = universe,
                             pvalueCutoff = p_value_cutoff)")
    } else if (db == 'WP'){
      cmd <- sprintf("compareCluster(geneClusters = markers, fun = enrichWP, 
                             organism = 'Homo sapiens', universe = universe,
                             pvalueCutoff = p_value_cutoff)")
    } else if (db == 'Reactome'){
      cmd <- sprintf("compareCluster(geneClusters = markers, 
                             fun = enrichPathway, readable = F,
                             universe = universe, 
                             pvalueCutoff = p_value_cutoff)")
    }
    
    tryCatch({
      
      # Running
      ck <- eval(parse(text=cmd))
      cat('    Success!\n')
      
      # Running the rest of the analysis
      ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
      
      cat(paste0("    Saving figures..", '\n'))
      plot_1 <- dotplot(ck) + theme(axis.text.y = element_text(size = 9)) 
      ggsave(plot_1, filename = sprintf('tmp/dotplot_%s.pdf', db), width=7, 
             height=9)
      
      plot_2 <- cnetplot(ck)
      ggsave(plot_2, filename = sprintf('tmp/cnetplot_%s.pdf', db))
      
      
      out[[db]] <- ck
    }, 
    
    error = function(e) {
      cat(sprintf("    Encountered error: '%s' when using %s annotations\n", e, db))
    }
    )
  }
  
  cat(paste0("\nSaving results..", '\n'))
  saveRDS(out, 'tmp/tmp_ora_from_DE_results.rds')
  cat('Done!\n')
  
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
