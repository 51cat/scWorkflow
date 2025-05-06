library(CellChat)
library(Seurat)
library(magrittr)
library(tidyverse)
library(SeuratObject)
library(uwot)
library(argparse)


dirnameScript <- function(){
  # get full directory path of current script located
  cmd = commandArgs(trailingOnly = FALSE)
  scriptName = sub('--file=', "", cmd[grep('^--file=', cmd)])
  if (length(scriptName) > 0) {
    path = normalizePath(scriptName)
    dirname = dirname(path)
  } else {
    print('Warning: not a runtime environment, using current directory instead.')
    dirname = getwd()
  }
  return(dirname)
}

##### 读入配色库 #####

source(str_glue("{dirnameScript()}/COLORS/load_color.r"))


save_gg <- function(p, filename, width, height, format = c('pdf', 'png')) {
  for (d in format) {
    ggplot2::ggsave(filename = str_glue("{filename}.{d}"), plot = p, width = width, height = height, device = d)
  }
}

mk.outdir <- function(dir) {
  if (!dir.exists(dir)){
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  }
  return(dir)
}



mk_input <- function(seurat.rds.path, group_col, celltype_col, control, treatment) {
  
  obj <- readRDS(seurat.rds.path)
  
  if (control != 'None' & treatment != 'None') {
    cells <- rownames(obj@meta.data[obj@meta.data[[group_col]] %in% c(control, treatment), ])
    obj <- subset(obj, cells = cells)
    ncell <- length(unique(obj@meta.data[[celltype_col]]))
    print(ncell)
  }
  
  objlist <- SplitObject(obj,split.by = group_col)
  print(objlist)
  
  # 筛选共有细胞
  all_celltypes <- lapply(objlist, function(seurat_obj) {
    unique(seurat_obj@meta.data[[celltype_col]]) 
  })
  
  common_celltypes <- Reduce(intersect, all_celltypes)
  
  print(all_celltypes)
  print(common_celltypes)
  
  seurat_list_filtered <- lapply(objlist, function(seurat_obj) {
    cell_ids <- rownames(seurat_obj@meta.data)
    
    # 找到符合共有细胞类型的细胞
    cells_to_keep <- cell_ids[seurat_obj@meta.data[[celltype_col]] %in% common_celltypes]
    
    subset(seurat_obj, cells = cells_to_keep)
  })
  
  # 输出细胞数统计信息
  print(lapply(seurat_list_filtered, function(x) {
    table(x@meta.data[[group_col]], x@meta.data[[celltype_col]])
  }))
  
  return(seurat_list_filtered)
}

process_data <- function(objlist, celltype_col, spec = 'human') {
  
  cellchatlist <- purrr::map(objlist,function(obj){
    data <- GetAssayData(object = obj, slot = 'data')
    meta <- obj@meta.data
    cellchat <- createCellChat(object = obj,
                               meta = meta,
                               group.by = celltype_col,
                               assay = 'RNA')
    if (spec == 'human'){
      cellchat@DB <- CellChatDB.human #human
    }else {
      cellchat@DB <- CellChatDB.mouse #mouse
    }
    
    cellchat <- subsetData(cellchat) #save time and memory
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    options(future.globals.maxSize = 2000 * 1024^2)
    cellchat <- computeCommunProb(cellchat)
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    #pathway level
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    gc()
    cellchat
  })
  
  cellchat <- mergeCellChat(cellchatlist, add.names = names(cellchatlist))
  cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
  #Compute signaling network similarity for datasets 1 2
  cellchat <- netEmbedding(cellchat, umap.method = 'uwot', type = "functional")
  #Manifold learning of the signaling networks for datasets 1 2
  cellchat <- netClustering(cellchat, type = "functional",do.parallel = F)
  
  return(list(
    merge = cellchat,
    obj_list = cellchatlist
  ))
}

extract_res_dataframe <- function(cellchat, outdir, prfx) {
  df.net<- subsetCommunication(cellchat)
  write.csv(df.net, str_glue("{outdir}/{prfx}_net_lr.csv"))
  
  df.netp<- subsetCommunication(cellchat,slot.name="netP")
  write.csv(df.netp,str_glue("{outdir}/{prfx}_net_pathway.csv"))
}

#### analysis function

show.Similarity <- function(cellchat_input, outdir, treatment, control) {
  cellchat_obj <- cellchat_input[['merge']]
  out.p <- mk.outdir(str_glue("{outdir}/{treatment}_vs_{control}/"))
  #Visualization in 2D-space
  p1 <- netVisual_embeddingPairwise(cellchat_obj, type = "functional", label.size = 2)
  #2D visualization of signaling networks from datasets 1 2
  p2 <- rankSimilarity(cellchat_obj, type  = "functional")
  #Compute the distance of signaling networks between datasets 1 2
  
  h <- min(0.11*(nrow(cellchat_obj@netP$similarity$functional$dr[['1-2']]) / 2), 10)
  w <- 3.65
  
  save_gg(p1, str_glue("{out.p}/functional_scatter.pdf"), 5, 5)
  save_gg(p2, str_glue("{out.p}/functional_rank_bar.pdf"), w, h)
  
  
}

plot_ccc_bar <- function(cellchat_obj, gr_level, gcpal) {
  df.count <- as.data.frame(sapply(cellchat_obj@net, function(x) sum(x$count)))
  colnames(df.count) <-  "val"
  df.count <- df.count %>% mutate(type='Number of inferred interactions ', group = rownames(.)) %>% as_tibble() 
  
  df.strength <- as.data.frame(sapply(cellchat_obj@net, function(x) sum(x$weight)))
  colnames(df.strength) <-  "val"
  df.strength <- df.strength %>% mutate(type='Interaction strength', group = rownames(.)) %>% as_tibble() 
  
  df <- rbind(df.count, df.strength)
  df$group <- factor(df$group, levels = gr_level)
  p <- ggplot(df, aes(group, val, fill = group))  + geom_bar(stat = 'identity', width = 0.75) + 
    scale_fill_manual(values = gcpal) +facet_wrap(~type, scale='free') +theme_classic()
  
  return(list(gg = p, table = df))
}

###### color  ######
show.compareInteractions <- function(cellchat_input, outdir, treatment, control, celltype_col, only_bar = FALSE, group_cpal = NULL, cell_cpal = NULL) {
  cellchat_obj <- cellchat_input[['merge']]
  cellchatlist <- cellchat_input[['obj_list']]
  ncell <- cellchat_obj@meta[[celltype_col]] %>% unique() %>% length
  
  res <- plot_ccc_bar(cellchat_obj, names(cellchatlist),  group_cpal)
  k1 <- res[['gg']]
  count.df <- res[['table']]
  
  # bar
  print(names(cellchatlist))
  
  out.p <- mk.outdir(str_glue("{outdir}/{treatment}_vs_{control}/compare/"))
  #gg1 <- compareInteractions(cellchat_obj, show.legend = F, group = names(cellchatlist), color.use = group_cpal) 
  #gg1 <- compareInteractions(cellchat_obj, show.legend = F , group = names(cellchatlist)) 
  #gg2 <- compareInteractions(cellchat_obj, show.legend = F, group = names(cellchatlist), measure = "weight", color.use = group_cpal)
  #gg2 <- compareInteractions(cellchat_obj, show.legend = F, measure = "weight", group =names(cellchatlist))
  #k1 <- gg1 + gg2
  save_gg(k1, str_glue("{out.p}/Interactions_compare_bar"), 8.8, 4)
  count.df  %>% write_csv(str_glue("{out.p}/Interactions_count.csv"))
  
  if (only_bar) {
    return(NULL)
  }
  
  # circle: 两个样本
  # https://stackoverflow.com/questions/49602032/how-to-convert-plot-to-ggplot-in-r
  
  k <- 0.9*ncell/10
  
  if (k < 1.5) {
    k <- 1.5
  }
  
  p1 <- ggplotify::as.ggplot(function() netVisual_diffInteraction(cellchat_obj, weight.scale = T, measure = "weight", vertex.label.cex = 0.85, color.use = cell_cpal))
  p2 <- ggplotify::as.ggplot(function() netVisual_diffInteraction(cellchat_obj, weight.scale = T, vertex.label.cex = 0.85, color.use = cell_cpal))
  save_gg(p1+p2, str_glue("{out.p}/Interactions_compare_circle_{treatment}_vs_{control}.pdf"), 10*k, 5*k)
  
  # heatmap: 两个样本
  k <- 0.6
  
  gg1 <-  ggplotify::as.ggplot(netVisual_heatmap(cellchat_obj, width = ncell*k, height = ncell*k, color.use = cell_cpal))
  gg2 <-  ggplotify::as.ggplot(netVisual_heatmap(cellchat_obj, measure = "weight", width= ncell*k, height = ncell*k, color.use = cell_cpal))
  k3<-gg1 + gg2
  save_gg(k3, str_glue("{out.p}/Interactions_compare_heatmap_{treatment}_vs_{control}.pdf"), ncell*k*2,ncell*k)
  
}

show.eachsample_interaction <- function(cellchat_input, celltype_col, outdir, cell_cpal = NULL) {
  cellchatlist <- cellchat_input[['obj_list']]
  cellchat_obj <- cellchat_input[['merge']]
  ncell <- cellchat_obj@meta[[celltype_col]] %>% unique() %>% length
  
  out.p <- mk.outdir(str_glue("{outdir}/sample_interaction/"))
  # 展示每个样本的互作强度和数量
  weight.max.number <- getMaxWeight(cellchatlist, attribute = c("idents","count"))
  weight.max.weight <- getMaxWeight(cellchatlist, attribute = c("idents","weight"))
  
  
  k <- 0.9*ncell/10
  
  if (k < 1.5) {
    k <- 1.5
  }
  
  p.list.number <- list()
  p.list.strength <- list()
  
  for (inx in 1:length(cellchatlist)) {
    p1 <- ggplotify::as.ggplot(function() netVisual_circle(cellchatlist[[inx]]@net$count, 
                                                           weight.scale = T, 
                                                           label.edge= F, 
                                                           color.use = cell_cpal,
                                                           vertex.label.cex = 0.85,
                                                           edge.weight.max = weight.max.number[2], 
                                                           title.name = paste0("Number of interactions - ", names(cellchatlist)[inx]))
    )
    
    p2 <- ggplotify::as.ggplot(function() netVisual_circle(cellchatlist[[inx]]@net$weight, 
                                                           weight.scale = T, 
                                                           label.edge= F, 
                                                           vertex.label.cex = 0.85,
                                                           color.use = cell_cpal,
                                                           edge.weight.max = weight.max.weight[2], 
                                                           title.name = paste0("Interaction weights/strength - ", names(cellchatlist)[inx]))
    )
    
    p.list.number[[inx]] <- p1
    p.list.strength[[inx]] <- p2
  }
  p.number <- ggpubr::ggarrange(plotlist = p.list.number, ncol = length(p.list.number))
  p.strength <- ggpubr::ggarrange(plotlist = p.list.strength, ncol =  length(p.list.strength))
  save_gg(p.number, str_glue("{out.p}/Interactions_number.pdf"), min(5*k*length(cellchatlist), 32), 5*k)
  save_gg(p.strength, str_glue("{out.p}/Interactions_strength.pdf"), min(5*k*length(cellchatlist), 32), 5*k)
}

###### color  ######
show.ranknet <- function(cellchat_input, outdir, group_cpal = 'npg') {
  out.p <- mk.outdir(str_glue("{outdir}"))
  cellchat <- cellchat_input[['merge']]
  n.dataset <- cellchat@meta$datasets %>% unique() %>% length
  dataset <- cellchat@meta$datasets %>% unique()
  
  total_pathway <- lapply(dataset, function(x) cellchat@netP[[x]]$pathways) %>% unlist() %>% unique %>% length
  print(n.dataset)
  
  if (total_pathway != 0) {
    gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, comparison = c(1:n.dataset), color.use = group_cpal) 
    gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, comparison = c(1:n.dataset), color.use = group_cpal )
    
    save_gg(gg1, str_glue("{out.p}/rank_Net_stacked"), 6, total_pathway*0.15)
    save_gg(gg2, str_glue("{out.p}/rank_Net"), 6, total_pathway*0.15)
  }
}


show.incoming_outcoming <- function(cellchat_input, outdir, celltype_col, cell_cpal = cell_cpal) {
  out.p <- mk.outdir(str_glue("{outdir}"))
  cellchatlist <- cellchat_input[['obj_list']]
  pathway.union  <- lapply(cellchatlist, function(x) x@netP$pathways) %>% unlist %>% unique
  
  for (target in c('outgoing', 'incoming', 'all')) {
    if (target == 'all') {color.heatmap = "OrRd"} else {color.heatmap = "GnBu"}
    plot_list <- list()
    for (inx in seq_along(cellchatlist)) {
      ncell <- cellchatlist[[inx]]@meta[[celltype_col]] %>% unique() %>% length
      
      cellchatlist[[inx]] <- netAnalysis_computeCentrality(cellchatlist[[inx]])
      ht <- netAnalysis_signalingRole_heatmap(
        cellchatlist[[inx]], pattern = target,
        color.use = cell_cpal, 
        signaling = pathway.union, title = names(cellchatlist)[inx],
        width = ncell*1.5, 
        height = length(pathway.union)*0.35, color.heatmap=color.heatmap)
      plot_list[[inx]] <- ggplotify::as.ggplot(ht)
    }
    combined_plot <- reduce(plot_list, `+`)
    save_gg(combined_plot, str_glue("{out.p}/Interactions_{target}"), min(16*length(combined_plot), 50), min(length(pathway.union)*0.5,50))
  }
}

show.pathwayGene_bubble <- function(cellchat_input, outdir, celltype_col, pathway_list = NULL) {
  out.p <- mk.outdir(str_glue("{outdir}/bubble/pathwaygene/"))
  cellchat <- cellchat_obj[['merge']]
  dataset <- cellchat@meta$datasets %>% unique()
  ncell <- cellchat@meta[[celltype_col]] %>% unique() %>% length()
  pathway.union <- lapply(dataset, function(x) cellchat@netP[[x]]$pathways) %>% unlist() %>% unique
  
  # 使用交集?
  # 
  pathway.union <- Reduce(intersect, lapply(dataset, function(x) cellchat@netP[[x]]$pathways)) %>% unique
  for (i in 1:length(pathway.union)){
    print(pathway.union[i])
    sub <- subsetCommunication(cellchat, signaling = pathway.union[i], slot.name = 'net')
    
    #ncell_pair <- lapply(dataset,  function(x)  str_c(sub[[x]]$source, sub[[x]]$target)) %>% unlist %>% unique %>% length
    #ngene <- lapply(dataset,  function(x)  sub[[x]]$interaction_name) %>% unlist %>% unique %>% length
    #print(ncell_pair)
    #print(ngene)
    res <- netVisual_bubble(cellchat, sources.use =c(1:ncell) , targets.use = c(1:ncell),  comparison = c(1, 2), signaling = pathway.union[i], remove.isolate = TRUE, return.data = TRUE)
    p <- res[['gg.obj']]
    df <- res[['communication']]
    ncell_pair <- df$source.target %>% unique() %>% length
    ngene <- df$interaction_name %>% unique() %>% length
    
    if (ncell_pair <= 10) {
      w <- 5
    } else {
      w <- min(20, ncell_pair*0.5)
    }
    
    if (ngene <= 10) {
      h <- 5
    } else {
      h <- min(20, ngene*0.3)
    }
    
    save_gg(p, str_glue("{out.p}/{pathway.union[i]}_gene_bubble.pdf"), w, h)
    df %>% write_csv(str_glue("{out.p}/{pathway.union[i]}_gene_bubble.csv"))
  }
  
}

show.diffpathway <- function(cellchat_input, outdir, treatment,  control, celltype_col, thresh.pc=0.1, thresh.fc = 0.1, thresh.p = 0.05) {
  out.p <- mk.outdir(str_glue("{outdir}/bubble/pathwayDEG/"))
  cellchat <- cellchat_input[['merge']]
  ncell <- cellchat@meta[[celltype_col]] %>% unique() %>% length()
  cellchat <- identifyOverExpressedGenes(
    cellchat, group.dataset = "datasets", 
    pos.dataset = treatment, features.name = treatment, only.pos = FALSE, thresh.pc = thresh.pc, thresh.fc = thresh.fc, thresh.p = thresh.p)
  net <- netMappingDEG(cellchat, features.name = treatment)
  
  net.up <- subsetCommunication(cellchat, net = net, datasets = treatment, ligand.logFC = thresh.fc, receptor.logFC =thresh.fc)
  net.down <- subsetCommunication(cellchat, net = net, datasets = control, ligand.logFC = -thresh.fc, receptor.logFC = -thresh.fc)
  
  gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
  gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
  
  pairLR.use.up = net.up[, "interaction_name", drop = F]
  pairLR.use.down = net.down[, "interaction_name", drop = F]
  
  if (nrow(pairLR.use.up) != 0) {
    print("up:")
    print(nrow(pairLR.use.up))
    
    res.up <- netVisual_bubble(
      cellchat, pairLR.use = pairLR.use.up, 
      sources.use = c(1:ncell), 
      targets.use = c(1:ncell), comparison = c(1, 2),  
      angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", treatment),
      return.data = TRUE)
    p.up <- res.up[['gg.obj']]
    df.up <- res.up[['communication']]
    ncell_pair <- df.up$source.target %>% unique() %>% length
    
    if (ncell_pair <= 10) {
      w <- 8
    } else {
      w <- min(20, ncell_pair*0.5)
    }
    
    if (nrow(pairLR.use.up) <= 10) {
      h <- 8
    } else {
      h <- min(20, nrow(pairLR.use.up)*0.3)
    }
    
    save_gg(p.up, str_glue("{out.p}/up_bubble.pdf"), w, h)
    df.up %>% write_csv(str_glue("{out.p}/up_gene_bubble.csv"))
    
    
  } 
  
  if (nrow(pairLR.use.down) != 0){
    print("down:")
    print(nrow(pairLR.use.down))
    res.down <- netVisual_bubble(
      cellchat, pairLR.use = pairLR.use.down, 
      sources.use =c(1:ncell), targets.use =c(1:ncell), comparison = c(1, 2),  
      angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", treatment),
      return.data = TRUE)
    
    p.down <- res.down[['gg.obj']]
    df.down <- res.down[['communication']]
    ncell_pair <- df.down$df.down %>% unique() %>% length
    
    if (ncell_pair <= 10) {
      w <- 5
    } else {
      w <- min(20, ncell_pair*0.5)
    }
    
    if (nrow(pairLR.use.down) <= 10) {
      h <- 5
    } else {
      h <- min(20, nrow(pairLR.use.down)*0.3)
    }
    
    save_gg(p.down, str_glue("{out.p}/down_bubble.pdf"), w, h)
    df.down %>% write_csv(str_glue("{out.p}/down_gene_bubble.csv"))
  }
}


get_pathway_cellchatlist <- function(pathway, cellchatlist) {
  new.list <- list()
  for (name in names(cellchatlist)) {
    target <- cellchatlist[[name]]
    if (pathway %in% target@netP$pathways) {
      new.list[[name]] <- target
    }
  }
  return(new.list)
}


show.pathway <- function(cellchat_input, outdir, celltype_col, group_col, cell_cpal = NULL){
  cellchatlist <- cellchat_input[['obj_list']]
  cellchat_obj <- cellchat_input[['merge']]
  
  ncell <- cellchat_obj@meta[[celltype_col]] %>% unique() %>% length
  
  out.p <- mk.outdir(str_glue("{outdir}/pathway_show/"))
  
  pathway.union <- lapply(cellchatlist, function(x) x@netP$pathways)%>% unlist %>% unique
  
  for (inx in seq_along(pathway.union)) {
    
    filter.cellchatlist <- get_pathway_cellchatlist(pathway.union[inx], cellchatlist)
    
    weight.max <- getMaxWeight(filter.cellchatlist, slot.name = c("netP"), attribute = pathway.union[inx])
    print(pathway.union[inx])
    print(weight.max[1])
    
    # circle
    p.circle.list <- lapply(filter.cellchatlist, function(x) {
      ggplotify::as.ggplot(function() netVisual_aggregate(
        x, 
        signaling = pathway.union[inx], 
        layout = "circle", 
        edge.weight.max = weight.max[1], 
        color.use = cell_cpal,
        edge.width.max = 10, 
        vertex.weight = NULL,
        vertex.label.cex = 0.85,
        signaling.name = str_c(pathway.union[inx], ":", unique(x@meta[[group_col]]))
      ))
    })
    
    
    k <- 0.9*ncell/10
    
    if (k < 1.5) {
      k <- 1.5
    }
    
    p.circle <- ggpubr::ggarrange(plotlist = p.circle.list,  ncol = length(p.circle.list))
    save_gg(p.circle, str_glue("{out.p}/{pathway.union[inx]}_circle"), min(5*k*length(filter.cellchatlist), 32), 5*k)
    
    
    # heatmap 
    ncell <- filter.cellchatlist[[1]]@meta[[celltype_col]] %>% unique() %>% length
    k <- 0.5
    p.heatmap.list <- lapply(filter.cellchatlist, function(x) ggplotify::as.ggplot(
      netVisual_heatmap(x, signaling = pathway.union[inx], color.use = cell_cpal, color.heatmap = "Reds",title.name = str_c(pathway.union[inx], ":" ,unique(x@meta[[group_col]])))
    ))
    p.hetmap <- ggpubr::ggarrange(plotlist = p.heatmap.list, ncol = length(p.heatmap.list))
    save_gg(p.hetmap, str_glue("{out.p}/{pathway.union[inx]}_heatmap"),ncell*k*2,ncell*k)
    
    
    # gene expression
    
    #levels <- cellchat@meta[[group_col]] %>% unique %>% as.character()
    #cellchat@meta[[group_col]] = factor(cellchat@meta[[group_col]], levels = levels)
    
    #p <- plotGeneExpression(cellchat, signaling = pathway.union[inx], split.by = group_col) 
    genes <- extractEnrichedLR(cellchat_obj, signaling = pathway.union[inx], geneLR.return = TRUE, enriched.only = TRUE)$geneLR
    genes.df <- tibble(gene = genes, pathway =  pathway.union[inx])
    genes.df %>% write_csv(str_glue("{out.p}/{pathway.union[inx]}_gene.csv"))
    #sufx <- stringr::str_c(levels, collapse = "_")
    #save_gg(p, str_glue("{out.p}/{pathway.union[inx]}_violin_gene_{sufx}"),ncell*0.8, length(genes)*1.1)
  }
}

show.netanalysis <- function(cellchat_input, outdir, cell_cpal = NULL) {
  out.p <- mk.outdir(str_glue("{outdir}"))
  cellchatlist <- cellchat_input[['obj_list']]
  num.link <- sapply(cellchatlist, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
  weight.MinMax <- c(min(num.link), max(num.link)) 
  
  gg <- list()
  for (i in 1:length(cellchatlist)) {
    gg[[i]] <- netAnalysis_computeCentrality(cellchatlist[[i]])
    gg[[i]] <- netAnalysis_signalingRole_scatter(gg[[i]], title = names(cellchatlist)[i], weight.MinMax = weight.MinMax, color.use = cell_cpal)
  }
  
  #Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
  p <- patchwork::wrap_plots(plots = gg, ncol = length(gg))
  # output
  n.plots <- length(gg)
  save_gg(p, str_glue("{out.p}/Interactions_compare_scatter"), 4*n.plots,4)
}

parse_comparisons <- function(input_str) {
  comparison_strs <- strsplit(input_str, ",")[[1]]
  parsed_list <- lapply(comparison_strs, function(x) {
    elements <- strsplit(trimws(x), "vs")[[1]]
    return(elements)
  })
  return(parsed_list)
}



parser <- ArgumentParser(description = "")

parser$add_argument("--rds", help = "")
parser$add_argument("--outdir", help = "")

parser$add_argument("--compare_str", help = "", default = 'None')
parser$add_argument("--compare_method", help = "", default = 'pairs')


parser$add_argument("--celltype_col", help = "")
parser$add_argument("--group_col", help = "")
parser$add_argument("--spec", help = "", default = 'human')

parser$add_argument("--group_sort", help = "", default = "None")
parser$add_argument("--gcpal", help = "", default = "npg")
parser$add_argument("--ccpal", help = "", default = "Paired")

parser$add_argument("--cellchat_rds_dir", help = "", default = "None")

parser$add_argument("--only_tables", help = "", default = "False")
parser$add_argument("--py_exec", help = "", default = "/usr/bin/py39")
args <- parser$parse_args()


Sys.setenv(RETICULATE_PYTHON=args$py_exec)

## 两两比较
seurat.rds.path <- args$rds
group_col <- args$group_col
celltype_col <- args$celltype_col
compare_str <- args$compare_str
outdir <- args$outdir
py_exec <- args$py_exec
only_tables <- args$only_tables


cellchat_rds_dir <- args$cellchat_rds_dir

compare_str <- args$compare_str
compare_method <- args$compare_method

group_sort <- args$group_sort
gcpal <- args$gcpal
ccpal <- args$ccpal
spec <- args$spec






# 设定颜色
rds.tmp <- readRDS(seurat.rds.path)
gr.use <- rds.tmp@meta.data[[group_col]] %>% unique()
cell.use <- rds.tmp@meta.data[[celltype_col]] %>% unique()

group_cpal <- get_color(gr.use , length(gr.use ), palette = gcpal)
cell_cpal <- get_color(cell.use , length(cell.use ), palette = ccpal)

print(group_cpal)
print(cell_cpal)


if (compare_method == 'pairs') {
  
  
  comps <- parse_comparisons(compare_str)
  for (comp in comps) {
    
    treatment <- comp[[1]]
    control <- comp[[2]]

    print(str_glue("treatment: {treatment}"))
    print(str_glue("control: {control}"))
    # save result
    outdir.comp <- str_glue("{outdir}/{treatment}_vs_{control}")
      
    out.p.df <- mk.outdir(str_glue("{outdir.comp}/table/"))
    res.out <- mk.outdir(str_glue("{outdir}/innput_obj/"))
    
    if (cellchat_rds_dir == 'None') {
      
      
      cellchat_input <- mk_input(
        seurat.rds.path, 
        group_col = group_col, 
        celltype_col=celltype_col, 
        control=control, 
        treatment=treatment)
      
      cellchat_obj <- process_data(
        cellchat_input,  
        celltype_col=celltype_col,  
        spec = spec )
      cellchat_obj[['run_mode']] <- 'two_sample'
      
      saveRDS(cellchat_obj, str_glue("{res.out}/{treatment}_vs_{control}.rds"))
      # 保存表
      for (obj in cellchat_obj[['obj_list']]) {
        extract_res_dataframe(obj, out.p.df, unique(obj@meta[[group_col]]))
      }
    }else {
      cellchat_obj <- readRDS(str_glue("{cellchat_rds_dir}/{treatment}_vs_{control}.rds"))
      print("LOAD : CELLCHAT RDS FINISH!")
    }
    
    
    # analysis
    show.Similarity(cellchat_obj, str_glue("{outdir.comp}"), treatment = treatment, control = control)
    show.compareInteractions(cellchat_obj,  str_glue("{outdir.comp}"), treatment = treatment, control = control, celltype_col= celltype_col, group_cpal = group_cpal, cell_cpal = cell_cpal)
    show.eachsample_interaction(cellchat_obj, str_glue("{outdir.comp}"), celltype_col = celltype_col, cell_cpal = cell_cpal)
    show.netanalysis(cellchat_obj,str_glue("{outdir.comp}"),cell_cpal = cell_cpal)
    show.ranknet(cellchat_obj, str_glue("{outdir.comp}"), group_cpal = group_cpal)
    
    show.incoming_outcoming(cellchat_obj, str_glue("{outdir.comp}"), celltype_col = celltype_col, cell_cpal = cell_cpal)
    show.pathwayGene_bubble(cellchat_obj, str_glue("{outdir.comp}"), celltype_col = celltype_col)
    show.diffpathway(cellchat_obj, str_glue("{outdir.comp}"), treatment=treatment,  control=control, celltype_col=celltype_col, thresh.pc=0.1, thresh.fc = 0.1, thresh.p = 0.05)
    show.pathway(cellchat_obj, str_glue("{outdir.comp}"),celltype_col = celltype_col, group_col = group_col, cell_cpal = cell_cpal)
    print('Finish!')
  }
}



if (compare_method == 'multi') {
  
  print('run: multi!')
  print(seurat.rds.path)
  outdir.comp <- str_glue("{outdir}/multi")
  res.out <- mk.outdir(str_glue("{outdir}/innput_obj/"))
  out.p.df <- mk.outdir(str_glue("{outdir.comp}/table/"))
  
  if (cellchat_rds_dir == 'None') {
    
    cellchat_input <- mk_input(seurat.rds.path, group_col = group_col, celltype_col=celltype_col, control='None', treatment='None')
    cellchat_obj <- process_data(cellchat_input,  celltype_col=celltype_col)
    cellchat_obj[['run_mode']] <- 'multi'
    
    if (group_sort != 'None') {
      group_sort <- str_split(group_sort, ',')[[1]]
      cellchat_obj$obj_list <- cellchat_obj$obj_list[group_sort]
      cellchat_obj$merge@meta$datasets <- factor(cellchat_obj$merge@meta$datasets, levels = group_sort)
    }
    
    saveRDS(cellchat_obj, str_glue("{res.out}/multi.rds"))
    
    for (obj in cellchat_obj[['obj_list']]) {
      extract_res_dataframe(obj, out.p.df, unique(obj@meta[[group_col]]))
    }
  }else {
    cellchat_obj <- readRDS(str_glue("{cellchat_rds_dir}/multi.rds"))
    print("LOAD : CELLCHAT RDS FINISH! - multi")
  }
  
  
  
  show.compareInteractions(cellchat_obj, outdir.comp, treatment = 'None', control= 'None', celltype_col = celltype_col, only_bar = TRUE,group_cpal = group_cpal)
  show.ranknet(cellchat_obj, str_glue("{outdir.comp}"), group_cpal = group_cpal)
  show.eachsample_interaction(cellchat_obj, str_glue("{outdir.comp}"), celltype_col = celltype_col, cell_cpal = cell_cpal)
  show.pathway(cellchat_obj, str_glue("{outdir.comp}"),celltype_col = celltype_col, group_col = group_col,cell_cpal = cell_cpal)
  show.netanalysis(cellchat_obj,str_glue("{outdir.comp}"), cell_cpal = cell_cpal)
  
}
print('Finish!')