library(monocle3)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(argparse)
library(ggVennDiagram)
library(colorRamp2)
library(RColorBrewer)
library(ComplexHeatmap)


set.seed(12306)

mk.outdir <- function(dir) {
  if (!dir.exists(dir)){
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  }
  return(dir)
}

save_gg <- function(p, filename, width, height, format = c('pdf', 'png')) {
  for (d in format) {
    ggplot2::ggsave(filename = str_glue("{filename}.{d}"), plot = p, width = width, height = height, device = d)
  }
}

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

##### 读入配色库/富集模块 #####

source(str_glue("{dirnameScript()}/COLORS/load_color.r"))
source(str_glue("{dirnameScript()}/enrich_func.r"))


# auto find root
get_earliest_principal_node <- function(cds, time_bin=root_celltype, cluster_col = cluster_col_use){
  
  cell_ids <- which(colData(cds)[, cluster_col] == time_bin)
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}


process_data <- function(scRNAsub, root_celltype, cluster_col, use_seurat_umap = FALSE, reduction_name = "umap") {
  # 更改sample_name列名, 如果不改会和monocle3冲突报错
  if ("sample_name" %in% colnames(scRNAsub@meta.data)) {
    print("sample_name columns in meta.data, auto change: sample_name --> sample_name_sufx")
    scRNAsub@meta.data <- scRNAsub@meta.data %>% rename(sample_name_sufx=sample_name)
  }
  
  data <- Seurat::GetAssayData(scRNAsub, assay = 'RNA', slot = 'counts')
  cell_metadata <- scRNAsub@meta.data
  gene_annotation <- data.frame(gene_short_name = rownames(data))
  rownames(gene_annotation) <- rownames(data)
  cds <- new_cell_data_set(data,
                           cell_metadata = cell_metadata,
                           gene_metadata = gene_annotation)
  
  #preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
  
  cds <- preprocess_cds(cds, num_dim = 150)
  #cds <- align_cds(cds, alignment_group = "orig.ident")
  
  cds <- cds[,Matrix::colSums(exprs(cds)) != 0]
  cds <- reduce_dimension(cds, preprocess_method = "PCA")
  
  if (use_seurat_umap) {
    print("use seurat umap..2")
    cds.embed <- cds@int_colData$reducedDims$UMAP
    int.embed <- Seurat::Embeddings(scRNAsub, reduction = reduction_name)
    int.embed <- int.embed[rownames(cds.embed),]
    cds@int_colData$reducedDims$UMAP <- int.embed
    print("use seurat umap..Finish")
  }
  
  cds <- cluster_cells(cds)
  cds <- learn_graph(cds)
  cds <- order_cells(cds, 
                     root_pr_nodes = get_earliest_principal_node(cds, root_celltype, cluster_col))
  
  return(cds)
}

show.tracre <- function(cds, group_col, celltype_col, outdir, group_cpal = 'npg', cell_cpal = 'Paired') {
  out.p <- mk.outdir(str_glue("{outdir}/tarce_show/"))
  
  # 设置颜色
  groups <- cds@colData[[group_col]] %>% unique
  cells <- cds@colData[[celltype_col]] %>% unique
  npartitions <- partitions(cds) %>% unique()
  
  group_cpal <- get_color(groups, length(groups), palette = group_cpal)
  cell_cpal <- get_color(cells, length(cells), palette = cell_cpal)
  #partition_cpal <- get_color(npartitions, length(npartitions), palette = "qPBI")
  
  # 绘制
  pp.func <- function(color_cells_by,...) {
    plot_cells(
      cds, 
      color_cells_by = color_cells_by, 
      label_cell_groups = FALSE, 
      label_leaves = FALSE,  
      label_branch_points = FALSE, 
      cell_size =0.75, 
      alpha  = 0.6) 
  }
  
  # partition
  p.partition <- pp.func(color_cells_by = "partition")
  p.partition.split <- p.partition + facet_wrap( as.formula(paste("~",group_col )),nrow = 1)
  
  # cell
  p.cell <- pp.func(color_cells_by = celltype_col) + 
    scale_color_manual(values = cell_cpal) 
  
  p.cell.split <- p.cell + facet_wrap( as.formula(paste("~",group_col )), nrow = 1)
  
  # group
  p.group <- pp.func(color_cells_by = group_col) + 
    scale_color_manual(values = group_cpal) 
  
  p.group.split <- p.group + facet_wrap( as.formula(paste("~",group_col )), nrow = 1)
  
  
  # pse
  p.pse <- pp.func(color_cells_by = 'pseudotime')  
  
  p.pse.split <- p.pse + facet_wrap( as.formula(paste("~",group_col )))
  
  
  save_gg(p.partition, str_glue("{out.p}/trace_partition"), 4.5, 4)
  save_gg(p.partition.split, str_glue("{out.p}/trace_partition_split"), 4*length(groups), 4)
  
  save_gg(p.cell, str_glue("{out.p}/trace_celltype"), 4.5, 4)
  save_gg(p.cell.split, str_glue("{out.p}/trace_celltype_split"), 4*length(groups), 4)
  
  save_gg(p.group, str_glue("{out.p}/trace_group"), 4.5, 4)
  save_gg(p.group.split, str_glue("{out.p}/trace_group_split"), 4*length(groups), 4)
  
  save_gg(p.pse, str_glue("{out.p}/trace_pseudotime"), 4.5, 4)
  save_gg(p.pse.split, str_glue("{out.p}/trace_pseudotime_split"), 4*length(groups), 4)
  
}


find_track_genes <- function(
    cds, 
    outdir,
    skip_process = FALSE,
    q_value_th = 1e-3,
    morans_I_th = 0,
    top = 10) {
  out.p <- mk.outdir(str_glue("{outdir}/gene_analysis/"))
  if (skip_process == FALSE) {
    Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=16)
    Track_genes %>% write_csv(str_glue('{out.p}/Track_genes.csv'))
  }else {
    Track_genes <- read_csv(str_glue('{out.p}/Track_genes.csv'))
  }
  
  
  
  Track_genes_top <- Track_genes %>% 
    filter(q_value <= q_value_th) %>% 
    filter(morans_I>morans_I_th)%>%
    top_n(n=top, morans_I) %>%
    pull(gene_short_name) %>% 
    as.character()
  
  Track_genes_all <- Track_genes %>% 
    filter(q_value <= q_value_th) %>% 
    filter(morans_I>morans_I_th)%>%
    pull(gene_short_name) %>% 
    as.character()
  
  
  return(
    list(all = Track_genes_all, top = Track_genes_top)
  )
}


get.val <- function(start, target = 'larger') {
  if (target == 'larger') {
    return(median(seq(start,1, length.out = 3)))
  } 
  
  if (target == 'less') {
    return(median(seq(0, start, length.out = 3)))
  }
}


find_target_modules_reso_val <- function(cds, genes, target_modules_number, init_reso, max_iter = 100) {
  print(str_glue("Trying resolution: {init_reso}"))
  
  gene_module_df <- find_gene_modules(cds[genes,], resolution = init_reso, cores = 32)
  nmodules <- gene_module_df %>% pull(module) %>% unique() %>% length()
  
  target_modules_number_max <- target_modules_number + 3
  target_modules_number_min <- target_modules_number - 1
  
  print(str_glue("Trying resolution: {init_reso}, found {nmodules} modules"))
  
  if (nmodules <= target_modules_number_max && nmodules >= target_modules_number_min) {
    return(gene_module_df)
  }
  
  if (max_iter <= 0) {
    warning("Maximum iterations reached, returning last result.")
    return(gene_module_df)
  }
  
  if (nmodules > target_modules_number_max) {
    init_reso <- init_reso - init_reso * 0.35
    #nit_reso <- get.val(init_reso, 'less')
    return(find_target_modules_reso_val(cds, genes, target_modules_number, init_reso, max_iter - 1))
  }
  
  if (nmodules < target_modules_number_min) {
    init_reso <- init_reso + init_reso * 0.35
    #init_reso <- get.val(init_reso, 'larger')
    return(find_target_modules_reso_val(cds, genes, target_modules_number, init_reso, max_iter - 1))
  }
}


gene.module.analysis <- function(
    cds,
    group_col,
    celltype_col,
    module_genes,
    show_genes,
    outdir,
    group_cpal = 'npg', 
    cell_cpal = 'Paired',
    resolution = 0.01
){
  out.p <- mk.outdir(str_glue("{outdir}/gene_analysis/"))
  out.p.show <- mk.outdir(str_glue("{outdir}/gene_analysis/gene_show_top10"))
  
  # 设置颜色
  groups <- cds@colData[[group_col]] %>% unique
  cells <- cds@colData[[celltype_col]] %>% unique
  npartitions <- partitions(cds) %>% unique()
  
  group_cpal <- get_color(groups, length(groups), palette = group_cpal)
  cell_cpal <- get_color(cells, length(cells), palette = cell_cpal)
  partition_cpal <- get_color(npartitions, length(npartitions), palette = "qPBI")
  
  
  pp.func <- function(color_cells_by, gene, ...) {
    plot_genes_in_pseudotime(cds[gene,], color_cells_by=color_cells_by, 
                             min_expr=0.5, ncol = 1, cell_size = 0.75)
  }
  
  for (gene in show_genes) {
    p1 <- pp.func(group_col, gene) + scale_color_manual(values = group_cpal)
    save_gg(p1, str_glue("{out.p.show}/{gene}_group"), 4, 2)
  }
  
  for (gene in show_genes) {
    
    p1 <- pp.func(celltype_col, gene) + scale_color_manual(values = cell_cpal)
    save_gg(p1, str_glue("{out.p.show}/{gene}_celltype"), 4, 2)
  }
  
  for (gene in show_genes) {
    p1 <- pp.func('pseudotime', gene)
    save_gg(p1, str_glue("{out.p.show}/{gene}_pseudotime"), 4, 2)
  }
  
  for (gene in show_genes) {
    p1 <- plot_cells(cds,
                     genes=gene,label_cell_groups=FALSE,
                     show_trajectory_graph=FALSE)
  }   
  
  # 模块分析
  print(resolution)
  print("++++++++")
  print(module_genes %>% head)
  #gene_module_df <- find_gene_modules(cds[module_genes,], resolution=resolution, cores = 32)
  
  gene_module_df <- find_target_modules_reso_val(
    cds,module_genes, length(cells), resolution, max_iter = 50)
  
  cell_group_df <- tibble::tibble(
    cell=row.names(colData(cds)),
    cell_group=clusters(cds)[colnames(cds)])
  
  ncell_group <- cell_group_df$cell_group %>% unique %>% length
  nmodules <- gene_module_df[['module']] %>% unique %>% length
  
  gene_module_df %>% write_csv(str_glue('{out.p}/gene_module.csv'))
  
  
  # 根据不同分组画热图
  ngroup <- colData(cds)[[group_col]] %>% unique %>% length
  
  if (ngroup <= 1) {
    g.use.all <- c(celltype_col, "all")
  }else  {
    g.use.all <- c(celltype_col, group_col, "all")
  }
  
  
  for (g.use in  g.use.all) {
    
    if (g.use != 'all'){
      cell_group_df <- tibble::tibble(
        cell=row.names(colData(cds)),
        cell_group=colData(cds)[[g.use]])
    }else {
      cell_group_df <- tibble::tibble(
        cell=row.names(colData(cds)),
        cell_group=str_glue("{colData(cds)[[celltype_col]]}:{colData(cds)[[group_col]]}"))
    }
    
    
    ncell_group <- cell_group_df$cell_group %>% unique %>% length
    nmodules <- gene_module_df[['module']] %>% unique %>% length
    
    agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
    row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
    colors <- get_color(c(1:1000), n = 1000, palette = "gsea")
    
    pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,color =  colors,fontsize_row = 6,
                       scale="column", clustering_method="ward.D2",treeheight_row  = 0,treeheight_col = 0, border_color = NA,
                       width = min(ncell_group*0.6, 12), height = min(nmodules*0.5, 30), filename = str_glue('{out.p}/module_{g.use}.pdf'))
    
    pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,color =  colors,fontsize_row = 6,
                       scale="column", clustering_method="ward.D2",treeheight_row  = 0,treeheight_col = 0, border_color = NA,
                       width = min(ncell_group*0.6, 12), height = min(nmodules*0.5, 30), filename = str_glue('{out.p}/module_{g.use}.png'))
    
    
  }
  
  # umap
  nmodules <- gene_module_df[['module']] %>% unique %>% length
  w1 <- min(sqrt(nmodules)*1.25, 20) %>% round(2)
  h1 <-  min(sqrt(nmodules)*1.25, 20) %>% round(2)
  
  p <- plot_cells(cds,
                  genes=gene_module_df,label_cell_groups=FALSE,
                  show_trajectory_graph=FALSE)
  
  save_gg(p, str_glue("{out.p}/module_UMAP"), w1, h1)
  
  # group
  gr.all <- colData(cds)[[group_col]] %>% unique()
  for (gr in gr.all) {
    cds_subset <- cds[, colData(cds)[[group_col]] == gr]
    p <- plot_cells(cds_subset,
                    genes=gene_module_df,label_cell_groups=FALSE,
                    show_trajectory_graph=FALSE)
    
    save_gg(p, str_glue("{out.p}/module_UMAP_{gr}"), w1, h1)
  }
  
}


#### 基因随拟时间变化分析 #### 


genes_pse_analysis2 <- function(cds, genes, k, show_genes = NULL) {
  #pseudotime_values <- pseudotime(cds)
  #valid_cells <- is.finite(pseudotime_values)
  #cds_filtered <- cds[, valid_cells]
  
  pt.matrix <- as.matrix(cds@assays@data$counts[match(genes,rowData(cds)[,1]),order(pseudotime(cds))])
  pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
  pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
  
  
  rownames(pt.matrix) <- genes
  
  
  ht <- Heatmap(
    pt.matrix,                     
    name                         = "z-score",
    km = k,
    col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
    show_row_names               = T,
    show_column_names            = FALSE,
    row_names_gp                 = gpar(fontsize = 6),
    row_title_rot                = 0, 
    cluster_rows                 = TRUE,
    cluster_row_slices           = FALSE,
    cluster_columns              = FALSE
  )
  
  ht_drawn <- draw(ht)
  row_order_list <- row_order(ht_drawn) 
  gene_cluster <- list()
  for (name in names(row_order_list)) {
    inx <- row_order_list[[name]]
    
    genes <- pt.matrix[inx,] %>% rownames()
    gene_cluster[[name]] <- genes
  }
  df <- enframe(gene_cluster, name = "cluster", value = "gene") %>% unnest_longer(gene) 
  sp <- df[['cluster']] %>% as.numeric()
  names(sp) <- df[['gene']]
  
  
  ht2 <- Heatmap(
    pt.matrix[names(sp),],                     
    name                         = "z-score",
    row_split = sp,
    col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
    show_row_names               = FALSE,
    show_column_names            = FALSE,
    row_names_gp                 = gpar(fontsize = 6),
    row_title_rot                = 0, 
    cluster_rows                 = FALSE,
    cluster_row_slices           = FALSE,
    cluster_columns              = FALSE
  )
  
  
  return(list(p = ht2, sp = sp, order= df))
}



get.top <- function(folder_path, keyword, top = 5) {
  pattern <- paste0("_", keyword, "_Results\\.csv$") 
  
  file_list <- list.files(path = folder_path, pattern = pattern, full.names = TRUE)
  
  all_data <- file_list %>%
    set_names() %>%  
    map_dfr(read_csv, .id = "source_file")
  
  all_data <- all_data %>%
    mutate(GeneRatio = map_dbl(GeneRatio, ~ {
      parts <- str_split(.x, "/", simplify = TRUE)
      as.numeric(parts[1]) / as.numeric(parts[2])
    }))
  all_data <- all_data %>% filter(p.adjust <= 0.05)
  
  top5_data <- all_data %>%
    group_by(group) %>%
    slice_max(order_by = GeneRatio, n = top, with_ties = FALSE) %>%
    ungroup() 
  top5.list <- split(top5_data$Description, top5_data$group)
  
  return(top5.list)
}


mk_ann <- function(sp, ann.lst) {
  
  ha <- rowAnnotation(textbox = anno_textbox(sp, ann.lst, word_wrap = TRUE, 
                                             add_new_line = TRUE))
  return(ha)
}



write_results_file <- function(cds, outdir) {
  out.p <- mk.outdir(str_glue("{outdir}"))
  save_monocle_objects(cds, str_glue("{outdir}/monocle3_result_cds"))
  #pseudotime <- pseudotime(cds, reduction_method = 'UMAP') %>% as.dataframe()
  
  write.csv(data.frame(barcode = rownames(cds@colData), time = pseudotime(cds)),
            file = str_glue("{outdir}/pseudotime_values.csv"), row.names = FALSE)
}


parser <- ArgumentParser(description = "")

parser$add_argument("--rds", help = "")
parser$add_argument("--cds", help = "", default = 'None')
parser$add_argument("--outdir", help = "")
parser$add_argument("--celltype_col", help = "")
parser$add_argument("--group_col", help = "")
parser$add_argument("--root_celltype", help = "")
parser$add_argument("--target_genes", help = "", default = NULL)

parser$add_argument("--use_seurat_umap", help = "", default = "True")
parser$add_argument("--reduction_name", help = "", default = 'umap')

parser$add_argument("--gcpal", help = "", default = "npg")
parser$add_argument("--ccpal", help = "", default = 'Paired')

parser$add_argument("--spec", help = "", default = 'human')

args <- parser$parse_args()

rds_path <- args$rds
outdir <- args$outdir
cluster_col <- args$celltype_col
group_col <- args$group_col
root_celltype <- args$root_celltype

use_seurat_umap <- args$use_seurat_umap
reduction_name <- args$reduction_name
target_genes <- args$target_genes
group_cpal <- args$gcpal
cell_cpal <- args$ccpal
spec <- args$spec



if (!is.null(target_genes)) {
  target_genes <- read_csv(target_genes)$gene %>% unique
  print(target_genes)
}

print("start...")

# 前处理
if (args$cds == 'None') {
  rds <- readRDS(rds_path)
  
  if (use_seurat_umap %in% c("True", "T", "TRUE")) {
    print("use seurat umap..1")
    cds <- process_data(rds,  root_celltype, cluster_col, use_seurat_umap = TRUE, reduction_name = reduction_name)
  }else {
    cds <- process_data(rds,  root_celltype, cluster_col, use_seurat_umap = FALSE)
  }
  
  write_results_file(cds, outdir)
}else {
  print('load CDS..')
  cds <- load_monocle_objects(args$cds)
  print('Finish load CDS!')
}

##########################


print('start show.tracre...')
show.tracre(cds, group_col, cluster_col, outdir,group_cpal = group_cpal, cell_cpal = cell_cpal)
print('Finish load show.tracre!')


# 基因分析

res <- find_track_genes(
  cds, 
  outdir,
  skip_process = FALSE,
  q_value_th = 1e-3,
  morans_I_th = 0,
  top = 10
)

# 获取top 轨迹基因和全部轨迹基因

top_track_gene <- res[['top']]
all_track_gene <- res[['all']]

print(top_track_gene)

# 基因模块分析

print('start gene.module.analysis...')
gene.module.analysis(
  cds,
  group_col,
  cluster_col,
  module_genes = all_track_gene,
  show_genes = top_track_gene,
  outdir,
  group_cpal = group_cpal, 
  cell_cpal = cell_cpal,
  resolution = 0.01)


if (!is.null(target_genes)) {
  
  out.p <- mk.outdir(str_glue("{outdir}/target_gene/"))
  
  res <- find_track_genes(
    cds, 
    outdir,
    skip_process = TRUE,
    q_value_th = 1e-3,
    morans_I_th = 0,
    top = 150
  )
  
  top_track_gene <- res[['top']]
  all_track_gene <- res[['all']]
  
  gene_list <- list(
    `Track` = top_track_gene,
    `Target` = target_genes
  )
  save_gg(ggVennDiagram(gene_list, set_size = 2.5), str_glue("{out.p}/veen_top"), 8, 4)
  
  gene_list <- list(
    `Track` = all_track_gene,
    `Target` = target_genes
  )
  save_gg(ggVennDiagram(gene_list, set_size = 2.5), str_glue("{out.p}/veen_all"), 8, 4)
  
  # 基因模块分析
  target_genes.show <- intersect(target_genes, top_track_gene)
  target_genes.module <- intersect(target_genes, all_track_gene)
  
  gene.module.analysis(
    cds,
    group_col,
    cluster_col,
    module_genes = target_genes.module,
    show_genes = target_genes.show,
    out.p,
    group_cpal = group_cpal, 
    cell_cpal = cell_cpal,
    resolution = 0.1)
  
}


### 富集分析 ###

if (spec == "human") {
  orgDb <- 'org.Hs.eg.db'
  organism <- 'hsa'
  ref <- 'human'
  
}else {
  orgDb <- 'org.Mm.eg.db'
  organism <- 'mmu'
  ref <- 'mouse'
}


### 基因module 文件
file.module <- str_glue('{outdir}/gene_analysis/gene_module.csv')
print(file.module)
out.p <- mk.outdir(str_glue('{outdir}/gene_analysis/module_enrich_results/'))

if (file.exists(file.module)) {
  module_df <- read_csv(file.module) %>% mutate(module = str_c('Module_', module))
  module_genes <- split(module_df$id, module_df$module)
  print(names(module_genes))
  
  
  for (name in names(module_genes)) {
    genes <- module_genes[[name]]
    print(name)
    # 富集分析
    analysis_GO(genes, orgDb = orgDb, output_prefix = str_glue("{out.p}/{name}_"), ont = 'BP')
    analysis_GO(genes, orgDb = orgDb, output_prefix = str_glue("{out.p}/{name}_"), ont = 'CC')
    analysis_GO(genes, orgDb = orgDb, output_prefix = str_glue("{out.p}/{name}_"), ont = 'MF')
    analysis_KEGG(genes, orgDb = orgDb, organism=organism , output_prefix = str_glue("{out.p}/{name}_"))
    analysis_Reactome(genes, orgDb = orgDb,organism = ref, output_prefix = str_glue("{out.p}/{name}_"))    
  }
  
}else {
  print("Not Found GENE MODULE file!")
}



file.module <- str_glue('{outdir}/target_gene/gene_analysis/gene_module.csv')
out.p <- mk.outdir(str_glue('{outdir}/target_gene/gene_analysis/module_enrich_results/'))

if (file.exists(file.module)) {
  module_df <- read_csv(file.module) %>% mutate(module = str_c('Module_', module))
  module_genes <- split(module_df$id, module_df$module)
  print(names(module_genes))
  
  for (name in names(module_genes)) {
    genes <- module_genes[[name]]
    # 富集分析
    analysis_GO(genes, orgDb = orgDb, output_prefix = str_glue("{out.p}/{name}_"), ont = 'BP')
    analysis_GO(genes, orgDb = orgDb, output_prefix = str_glue("{out.p}/{name}_"), ont = 'CC')
    analysis_GO(genes, orgDb = orgDb, output_prefix = str_glue("{out.p}/{name}_"), ont = 'MF')
    analysis_KEGG(genes, orgDb = orgDb, organism=organism , output_prefix = str_glue("{out.p}/{name}_"))
    analysis_Reactome(genes, orgDb = orgDb,organism = ref, output_prefix = str_glue("{out.p}/{name}_"))    
  }
  
}else {
  print("Not Found GENE MODULE file!")
}

###################################


#### 基因随拟时间变化分析 #### 

## top 500 的轨迹基因分析 
file.track <- str_glue('{outdir}/gene_analysis/Track_genes.csv')
out.p <- mk.outdir(str_glue('{outdir}/gene_analysis/gene_pse_analysis/'))
out.p.enrich <- mk.outdir(str_glue('{outdir}/gene_analysis/gene_pse_analysis/enrich/'))

if (file.exists(file.track)) {
  Track_genes <- read_csv(file.track)
  
  Track_genes_top <- Track_genes %>% 
    filter(q_value <= 1e-3) %>% 
    filter(morans_I>0)%>%
    top_n(n=500, morans_I) %>%
    pull(gene_short_name) %>% 
    as.character()
  
  ncell <- colData(cds)[[cluster_col]] %>% unique() %>% length()
  
  res <- genes_pse_analysis2(cds, Track_genes_top, ncell)
  
  pdf(str_glue('{out.p}/top500_gene_pse_exp.pdf'),height = 9, width = 6,onefile = F)
  print(res[['p']])
  res[['order']] %>% write_csv(str_glue('{out.p}/gene_cluster.csv'))
  dev.off()
  
  # 富集分析
  module_genes <- split(res[['order']]$gene, res[['order']]$cluster)
  
  for (name in names(module_genes)) {
    print(name)
    genes <- module_genes[[name]]
    # 富集分析
    analysis_GO(genes, orgDb = orgDb, output_prefix = str_glue("{out.p.enrich}/{name}_"), ont = 'BP', gr_name = name)
    analysis_GO(genes, orgDb = orgDb, output_prefix = str_glue("{out.p.enrich}/{name}_"), ont = 'CC', gr_name = name)
    analysis_GO(genes, orgDb = orgDb, output_prefix = str_glue("{out.p.enrich}/{name}_"), ont = 'MF', gr_name = name)
    analysis_KEGG(genes, orgDb = orgDb, organism=organism , output_prefix = str_glue("{out.p.enrich}/{name}_"), gr_name = name)
    analysis_Reactome(genes, orgDb = orgDb,organism = ref, output_prefix = str_glue("{out.p.enrich}/{name}_"), gr_name = name)    
  }
  
  
}else {
  print('NOT found track gene file')
}

### 标记top item
# GO
BP.list <- get.top(out.p.enrich, 'BP')
CC.list <- get.top(out.p.enrich, 'CC')
MF.list <- get.top(out.p.enrich, 'MF')
kegg.list <- get.top(out.p.enrich, 'kegg')
reactome.list <- get.top(out.p.enrich, 'reactome')

sp <- res[['sp']]

print(sp)
print(BP.list)
print(mk_ann(sp, BP.list))


p.all <- list(
  BP = res[['p']] + mk_ann(sp, BP.list),
  CC = res[['p']] + mk_ann(sp, CC.list),
  MF = res[['p']] +mk_ann(sp, MF.list),
  kegg = res[['p']] + mk_ann(sp, kegg.list),
  reactome = res[['p']] +mk_ann(sp, reactome.list)
)

### 热图

for (p.name in names(p.all)) {
  pdf(str_glue('{out.p}/top500_gene_pse_exp_{p.name}.pdf'),height = 9, width = 10,onefile = F)
  print(p.all[[p.name]])
  dev.off()
}


###  针对所有 common target gene的分析 ### 
file.track <- str_glue('{outdir}/target_gene/gene_analysis/gene_module.csv')
out.p <- mk.outdir(str_glue('{outdir}/target_gene/gene_analysis/gene_pse_analysis/'))
out.p.enrich <- mk.outdir(str_glue('{outdir}/target_gene/gene_analysis/gene_pse_analysis/enrich/'))

if (file.exists(file.track)) {
  Track_genes <- read_csv(file.track) %>% pull(id) %>% unique()
  ncell <- colData(cds)[[cluster_col]] %>% unique() %>% length()
  
  res <- genes_pse_analysis2(cds, Track_genes_top, ncell)
  
  pdf(str_glue('{out.p}/target_gene_pse_exp.pdf'),height = 9, width = 6,onefile = F)
  print(res[['p']])
  res[['order']] %>% write_csv(str_glue('{out.p}/target_gene_cluster.csv'))
  dev.off()
  
  # 富集分析
  module_genes <- split(res[['order']]$gene, res[['order']]$cluster)
  
  for (name in names(module_genes)) {
    print(name)
    genes <- module_genes[[name]]
    # 富集分析
    analysis_GO(genes, orgDb = orgDb, output_prefix = str_glue("{out.p.enrich}/{name}_"), ont = 'BP', gr_name = name)
    analysis_GO(genes, orgDb = orgDb, output_prefix = str_glue("{out.p.enrich}/{name}_"), ont = 'CC', gr_name = name)
    analysis_GO(genes, orgDb = orgDb, output_prefix = str_glue("{out.p.enrich}/{name}_"), ont = 'MF', gr_name = name)
    analysis_KEGG(genes, orgDb = orgDb, organism=organism , output_prefix = str_glue("{out.p.enrich}/{name}_"), gr_name = name)
    analysis_Reactome(genes, orgDb = orgDb,organism = ref, output_prefix = str_glue("{out.p.enrich}/{name}_"), gr_name = name)    
  }
  ### 标记top item
  # GO
  BP.list <- get.top(out.p.enrich, 'BP')
  CC.list <- get.top(out.p.enrich, 'CC')
  MF.list <- get.top(out.p.enrich, 'MF')
  kegg.list <- get.top(out.p.enrich, 'kegg')
  reactome.list <- get.top(out.p.enrich, 'reactome')
  
  sp <- res[['sp']]
  
  p.all <- list(
    BP = res[['p']] + mk_ann(sp, BP.list),
    CC = res[['p']] +mk_ann(sp, CC.list),
    MF = res[['p']] +mk_ann(sp, MF.list),
    kegg = res[['p']] +mk_ann(sp, kegg.list),
    reactome = res[['p']] +mk_ann(sp, reactome.list)
  )
  
  ### 热图
  
  for (p.name in names(p.all)) {
    pdf(str_glue('{out.p}/target_gene_pse_exp_{p.name}.pdf'),height = 9, width = 10,onefile = F)
    print(p.all[[p.name]])
    dev.off()
  }
  
  
}else {
  print('NOT found track gene file')
}

