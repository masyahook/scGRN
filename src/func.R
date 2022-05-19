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
  library(rlang)
  if (version$minor == '2.0'){
    library(ggnewscale)
  }
  library(stringr)
  library(reticulate)
})

# Define colors
colors <- list(green='#39B600', yellow='#D89000', red='#F8766D', blue='#00B0F6', 
               purple='#9590FF', cyan='#00BFC4', pink='E76BF3', light_pink='#FF62BC',
               saturated_green='#00BF7D')

# Define some useful functions
get_cl_from_tick <- function (x) substr(x, 1, as.integer(gregexpr('\n', x)) - 1)
cl_name_to_int <- function (x) substr(x, as.integer(gregexpr('_', x)) + 1, nchar(x))

compute_ranking <- function(df, rank_type='FC'){
  if (rank_type == 'FC'){
     out <- df['avg_log2FC']
  } else if (rank_type == 'signed_p'){
    out <- df['avg_log2FC'] * (-log10(df['p_val_adj']))
    index <- is.infinite(out[['avg_log2FC']])
    out[index, 'avg_log2FC'] <- max(out[['avg_log2FC']][is.finite(out[['avg_log2FC']])], na.rm=T)
  } else if (rank_type == 'centrality'){
    out <- df['centrality']
  } else {
    stop(sprintf("rank_type '%s' is incorrect, please change..", rank_type), call.=FALSE)
  }
  
  return(out)
}

remove_geom <- function(ggplot2_object, geom_type) {
  # Delete layers that match the requested type.
  layers <- lapply(ggplot2_object$layers, function(x) {
    if (class(x$geom)[1] == geom_type) {
      NULL
    } else {
      x
    }
  })
  # Delete the unwanted layers.
  layers <- layers[!sapply(layers, is.null)]
  ggplot2_object$layers <- layers
  ggplot2_object
}

remove_legend_title <- function(ggplot2_object) {
  # Delete the unwanted layers.
  layers <- ggplot2_object$layers[-7]
  ggplot2_object$layers <- layers
  ggplot2_object
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

run_gsea <- function(markers, is_clusters=F, top_n=10, db_run='ALL', 
                                 from_mat=T, rank_type='signed_p', 
                                 p_value_cutoff=1, ont='BP'){
  
  # Setting some params
  selected_cols <- c("ID", "Description", "p.adjust")
  dbs <- c('GO', 'KEGG', 'WP', 'Reactome')
  groups <- unique(markers$cluster)
  
  if (is_clusters == T){
    rank_type <- 'centrality'
  }
  
  # Running GSEA on each group and saving all results
  all_res <- list()
  cat('\n\nProcessing..\n\n\n')
  for (group in groups){
    cat(paste0(sprintf("Finding enrichment terms for group: '%s'", group), '\n'))
    group_markers <- markers[markers$cluster == group,]
    all_res[[group]] <- gsea_cP(group_markers, db=db_run, from_mat=from_mat, 
                              rank_type=rank_type, ont=ont,
                              p_value_cutoff=p_value_cutoff)
    cat('\n\n')
  }
  
  # Saving shortened version of results
  short_res <- list()
  cat('\nGetting shortened results..\n\n')
  for (db in dbs){

    curr_short_res <- data.frame(row.names=1:top_n)
    for (group in groups){
      
      if (group %in% names(all_res)){
        if (db %in% names(all_res[[group]])){
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
      }

    }

    short_res[[db]] <- curr_short_res
  }

  cat('Done!\n')

  out <- list(short=short_res, all=all_res)
  cat("\nSaving results to 'tmp/tmp_gsea.rds'..\n")
  saveRDS(out, 'tmp/tmp_gsea.rds')
  cat('Done!\n')
  
  return(out)
}

run_ora <- function(markers_df, is_clusters=F, top_n_dotplot = 10, 
                    top_n_cnetplot = 5, top_logFC = 1, top_p_val = 0.05, 
                    cell_type_for_community_ana = NaN, ytick_dotplot_size = 8,
                    xtick_dotplot_size = 8, p_value_cutoff=0.05, ont='BP', 
                    suffix=''){
  
  dbs <- c('GO', 'KEGG', 'WP', 'Reactome')
  pd <- import("pandas")
  
  # Setting some params
  if (is_clusters == T){
    groups <- unique(markers_df$cluster)
    universe = unique(sort(as.data.frame(org.Hs.egGO)$gene_id))
    
    final_markers_df <- markers_df
  } else {
    pat_type_levels <- c('C', 'M', 'S')
    dis_type_levels <- c('H', 'C')
    new_type_levels <- c('NS', 'S')
    pat_levels <- c('C51', 'C52', 'C100', 'C141', 'C142', 'C144', 'C143', 'C145', 'C146', 'C148', 'C149', 'C152')
    if (length(intersect(pat_levels, unique(markers_df$cluster))) == 12){
      groups <- pat_levels
    } else if (length(intersect(pat_type_levels, unique(markers_df$cluster))) == 3){
      groups <- pat_type_levels
    } else if (length(intersect(dis_type_levels, unique(markers_df$cluster))) == 2){
      groups <- dis_type_levels
    } else {
      groups <- new_type_levels
    }
    
    # Other
    final_markers_df <- markers_df[(markers_df$p_val_adj < top_p_val) & 
                                     (markers_df$avg_log2FC > top_logFC),]
    tmp_universe <- unique(markers_df$gene)
    universe <- bitr(tmp_universe, fromType = 'SYMBOL', toType = c('ENTREZID'),
                     OrgDb = org.Hs.eg.db)$ENTREZID
  }
  
  # Getting gene marker lists
  markers <- list()
  symbol_markers <- list()
  for (group in groups){
    tryCatch({
      tmp_markers <- final_markers_df[final_markers_df$cluster == group,]
      tmp_df <- bitr(tmp_markers$gene, fromType = "SYMBOL",
                     toType = c("ENTREZID"), 
                     OrgDb = org.Hs.eg.db)
      markers[[group]] <- tmp_df$ENTREZID
      symbol_markers[[group]] <- tmp_markers$gene
    }, 
      error = function(e) {
        cat(sprintf("    Encountered error: '%s' when getting ENTREZ IDs for %s\n", e, group))
        markers[[group]] <- c()
        symbol_markers[[group]] <- tmp_markers$gene
    })
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
      
      if (is_clusters == T){
        x_title <- 'Cluster'
      } else (
        x_title <- 'Patient type'
      )
      
      if (is_clusters == T){
        cell_type = str_split(cell_type_for_community_ana, pattern=' ')[[1]][1]
        data_type <- str_split(cell_type_for_community_ana, pattern=' ')[[1]][2]
        if (data_type == 'all'){
          df_name <- 'raw_data_communities_info.pickle'
        } else {
          df_name <- sprintf('raw_data_%s_type_communities_info.pickle', data_type)
        }
        path_to_df <- sprintf(
          paste0('/gpfs/projects/bsc08/bsc08890/res/covid_19/cell_types/%s', 
                '/data/grnboost2/leiden_communities/%s'),
          cell_type, df_name
        )
        df <- pd$read_pickle(path_to_df)
        num_nodes_per_cluster <- unlist(df$num_nodes)
        names(num_nodes_per_cluster) <- lapply(
          0:(length(num_nodes_per_cluster) - 1), 
          function(x) sprintf('cluster_%d', x))
      }
      
      if (db == 'Reactome'){
        ck@compareClusterResult$Description <- lapply(
          ck@compareClusterResult$Description,
          function (x) str_replace(x, 'Home Sapiens\r: ', '')
        )
      }
      
      cat(paste0("    Saving figures..", '\n'))
      
      plot_1 <- dotplot(ck, showCategory=top_n_dotplot) + 
        theme(axis.text.y = element_text(size = 9)) + 
        labs(y = 'Associated pathways', x = x_title)
      if (is_clusters == T){
        new_xticks <- lapply(ggplot_build(plot_1)$layout$panel_params[[1]]$x.sec$scale$range$range, 
                            function (x) paste0(cl_name_to_int(get_cl_from_tick(x)), '\n', sprintf('(%s)', num_nodes_per_cluster[[get_cl_from_tick(x)]])))
        plot_1 <- plot_1 + scale_x_discrete(labels=new_xticks) +
          theme(axis.text.y = element_text(size = ytick_dotplot_size), 
                axis.text.x = element_text(size = xtick_dotplot_size))
      }
      ggsave(plot_1, filename = sprintf('tmp/dotplot_%s_%s.pdf', suffix, db), 
             width=7, 
             height=9)
      
      if (length(groups) == 3){
        plot_2 <- cnetplot(ck, showCategory=top_n_cnetplot) + 
          scale_fill_manual(values=c(colors$green, colors$yellow, colors$red))
        plot_2 <- remove_legend_title(plot_2)
      } else if (length(groups) == 2){
        plot_2 <- cnetplot(ck, showCategory=top_n_cnetplot) + 
          scale_fill_manual(values=c(colors$green, colors$yellow))
        plot_2 <- remove_legend_title(plot_2)
      } else {
        plot_2 <- cnetplot(ck, cex_label_category=0.8, cex_label_gene=0.6, showCategory=top_n_cnetplot)
        new_labels <- lapply(ggplot_build(plot_1)$layout$panel_params[[1]]$x.sec$scale$range$range, 
                             function (x) paste0(cl_name_to_int(get_cl_from_tick(x)), ' ', sprintf('(%s)', num_nodes_per_cluster[[get_cl_from_tick(x)]])))
        plot_2 <- plot_2 + scale_fill_discrete(labels=new_labels)
        plot_2 <- remove_legend_title(plot_2)
      }
      
      ggsave(plot_2, filename = sprintf('tmp/cnetplot_%s_%s.pdf', suffix, db),
             width=7, height=7)
      
      
      out[[db]] <- ck
    }, 
    
    error = function(e) {
      cat(sprintf("    Encountered error: '%s' when using %s annotations\n", e, db))
    }
    )
  }
  
  cat(sprintf("\nSaving results to 'tmp/tmp_ora_%s.rds'..\n", suffix))
  saveRDS(out, sprintf('tmp/tmp_ora_%s.rds', suffix)) 
  cat('Done!\n')
  
  return(out)
  
}

ora_for_wordcloud <- function(markers_df, db='Reactome', p_value_cutoff=0.05,
                              ont='BP'){
  
  groups <- unique(markers_df$cluster)
  universe = unique(sort(as.data.frame(org.Hs.egGO)$gene_id))
  
  # Getting gene marker lists
  markers <- list()
  symbol_markers <- list()
  for (group in groups){
    tryCatch({
      tmp_markers <- markers_df[markers_df$cluster == group,]
      tmp_df <- bitr(tmp_markers$gene, fromType = "SYMBOL",
                     toType = c("ENTREZID"), 
                     OrgDb = org.Hs.eg.db)
      markers[[group]] <- tmp_df$ENTREZID
      symbol_markers[[group]] <- tmp_markers$gene
    }, 
    error = function(e) {
      cat(sprintf("    Encountered error: '%s' when getting ENTREZ IDs for %s\n", e, group))
      markers[[group]] <- c()
      symbol_markers[[group]] <- tmp_markers$gene
    })
  }
  
  cat('\n\nProcessing..\n\n\n')
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
    
    # Getting the output data
    ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    
    if (db == 'Reactome'){
      ck@compareClusterResult$Description <- lapply(
        ck@compareClusterResult$Description, 
        function (x) str_replace(x, 'Home Sapiens\r: ', '')
      )
    }
    
    out <- as.data.frame(ck)
    out <- data.frame(
      Cluster=unlist(out$Cluster), Description=unlist(out$Description)
    )
  }, 
  
  error = function(e) {
    cat(sprintf("    Encountered error: '%s' when using %s annotations\n", e, db))
  }
  )
  return(out)
}

DoMultiBarHeatmap <- function (object, 
                               features = NULL, 
                               cells = NULL, 
                               group.by = "ident", 
                               additional.group.by = NULL, 
                               group.bar = TRUE, 
                               disp.min = -2.5, 
                               disp.max = NULL, 
                               slot = "scale.data", 
                               assay = NULL, 
                               label = TRUE, 
                               size = 5.5, 
                               hjust = 0, 
                               angle = 45, 
                               raster = TRUE, 
                               draw.lines = TRUE, 
                               lines.width = NULL, 
                               group.bar.height = 0.02, 
                               combine = TRUE) 
{
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  features <- features %||% VariableFeatures(object = object)
  ## Why reverse???
  features <- rev(x = unique(x = features))
  disp.max <- disp.max %||% ifelse(test = slot == "scale.data", 
                                   yes = 2.5, no = 6)
  possible.features <- rownames(x = GetAssayData(object = object, 
                                                 slot = slot))
  if (any(!features %in% possible.features)) {
    bad.features <- features[!features %in% possible.features]
    features <- features[features %in% possible.features]
    if (length(x = features) == 0) {
      stop("No requested features found in the ", slot, 
           " slot for the ", assay, " assay.")
    }
    warning("The following features were omitted as they were not found in the ", 
            slot, " slot for the ", assay, " assay: ", paste(bad.features, 
                                                             collapse = ", "))
  }
  data <- as.data.frame(x = as.matrix(x = t(x = GetAssayData(object = object, 
                                                             slot = slot)[features, cells, drop = FALSE])))
  
  object <- suppressMessages(expr = StashIdent(object = object, 
                                               save.name = "ident"))
  group.by <- group.by %||% "ident"
  groups.use <- object[[c(group.by, additional.group.by)]][cells, , drop = FALSE]
  plots <- list()
  for (i in group.by) {
    data.group <- data
    group.use <- groups.use[, c(i, additional.group.by), drop = FALSE]
    
    for(colname in colnames(group.use)){
      if (!is.factor(x = group.use[[colname]])) {
        group.use[[colname]] <- factor(x = group.use[[colname]])
      }  
    }
    
    if (draw.lines) {
      lines.width <- lines.width %||% ceiling(x = nrow(x = data.group) * 
                                                0.0025)
      placeholder.cells <- sapply(X = 1:(length(x = levels(x = group.use[[i]])) * 
                                           lines.width), FUN = function(x) {
                                             return(Seurat:::RandomName(length = 20))
                                           })
      placeholder.groups <- data.frame(foo=rep(x = levels(x = group.use[[i]]), times = lines.width))
      placeholder.groups[additional.group.by] = NA
      colnames(placeholder.groups) <- colnames(group.use)
      rownames(placeholder.groups) <- placeholder.cells
      
      group.levels <- levels(x = group.use[[i]])
      
      group.use <- sapply(group.use, as.vector)
      rownames(x = group.use) <- cells
      
      group.use <- rbind(group.use, placeholder.groups)
      
      na.data.group <- matrix(data = NA, nrow = length(x = placeholder.cells), 
                              ncol = ncol(x = data.group), dimnames = list(placeholder.cells, 
                                                                           colnames(x = data.group)))
      data.group <- rbind(data.group, na.data.group)
    }
    
    #group.use = group.use[order(group.use[[i]]), , drop=F]
    group.use <- group.use[with(group.use, eval(parse(text=paste('order(', paste(c(i, additional.group.by), collapse=', '), ')', sep='')))), , drop=F]
    
    plot <- Seurat:::SingleRasterMap(data = data.group, raster = raster, 
                                     disp.min = disp.min, disp.max = disp.max, feature.order = features, 
                                     cell.order = rownames(x = group.use), group.by = group.use[[i]])
    
    if (group.bar) {
      pbuild <- ggplot_build(plot = plot)
      group.use2 <- group.use
      cols <- list()
      na.group <- Seurat:::RandomName(length = 20)
      for (colname in colnames(group.use2)){
        if (colname == group.by){
          colid = paste0('Identity (', colname, ')')
        } else {
          colid = colname
        }
        #Create labels
        levels(x = group.use2[[colname]]) <- unique(group.use2[[colname]])
        if (draw.lines) {
          levels(x = group.use2[[colname]]) <- c(levels(x = group.use2[[colname]]), na.group)  
          group.use2[placeholder.cells, colname] <- na.group
          cols[[colname]] <- c(scales::hue_pal()(length(x = levels(x = group.use2[[colname]]))-1), "#FFFFFF")
        } else {
          cols[[colname]] <- c(scales::hue_pal()(length(x = levels(x = group.use2[[colname]]))))
        }
        names(x = cols[[colname]]) <- levels(x = group.use2[[colname]])
        
        
        y.range <- diff(x = pbuild$layout$panel_params[[1]]$y.range)
        y.pos <- max(pbuild$layout$panel_params[[1]]$y.range) + y.range * 0.015
        y.max <- y.pos + group.bar.height * y.range
        pbuild$layout$panel_params[[1]]$y.range <- c(pbuild$layout$panel_params[[1]]$y.range[1], y.max)
        
        plot <- suppressMessages(plot + 
                                   annotation_raster(raster = t(x = cols[[colname]][group.use2[[colname]]]), xmin = -Inf, xmax = Inf, ymin = y.pos, ymax = y.max) + 
                                   annotation_custom(grob = grid::textGrob(label = colid, hjust = 0, gp = grid::gpar(cex = 0.75)), ymin = mean(c(y.pos, y.max)), ymax = mean(c(y.pos, y.max)), xmin = Inf, xmax = Inf) +
                                   coord_cartesian(ylim = c(0, y.max), clip = "off"))
        
        #temp <- as.data.frame(cols[[colname]][levels(group.use[[colname]])])
        #colnames(temp) <- 'color'
        #temp$x <- temp$y <- 1
        #temp[['name']] <- as.factor(rownames(temp))
        
        #temp <- ggplot(temp, aes(x=x, y=y, fill=name)) + geom_point(shape=21, size=5) + labs(fill=colname) + theme(legend.position = "bottom")
        #legend <- get_legend(temp)
        #multiplot(plot, legend, heights=3,1)
        
        if ((colname != group.by) && label) {
          x.max <- max(pbuild$layout$panel_params[[1]]$x.range)
          x.divs <- pbuild$layout$panel_params[[1]]$x$break_positions()
          group.use$x <- x.divs
          label.unique <- paste(group.use[[colname]], group.use[[group.by]], sep="+=$")
          label.x.pos <- tapply(X = group.use$x, INDEX = label.unique,
                                FUN = median) * x.max
          label.x.pos <- data.frame(group = names(x = label.x.pos), 
                                    label.x.pos)
          label.x.pos <- label.x.pos[!grepl("NA", label.x.pos$group),]
          label.x.pos$group <- label.x.pos$group %>% lapply( function(x) gsub("\\+\\=\\$.*","", x))
          plot <- plot + geom_text(stat = "identity", 
                                   data = label.x.pos, aes_string(label = "group", 
                                                                  x = "label.x.pos"), y = y.max + y.max * 
                                     0.03 * 0.5, angle = angle, hjust = hjust, 
                                   size = size)
          
          plot <- suppressMessages(plot + coord_cartesian(ylim = c(0, 
                                                                   y.max + y.max * 0.002 * max(nchar(x = unique(group.use[[colname]])), na.rm=TRUE) * 
                                                                     size), clip = "off"))
        }
      }
    }
    plot <- plot + theme(line = element_blank())
    plots[[i]] <- plot
  }
  if (combine) {
    plots <- CombinePlots(plots = plots)
  }
  return(plots)
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
