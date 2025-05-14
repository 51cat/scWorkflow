library(SCP)
library(scCustomize)
library(ggplot2)
library(tidyverse)
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

source(str_glue("{dirnameScript()}/COLORS/load_color.r"))



fetch.seurat <- function(rds, group_col, target) {
  cell.use <- rownames(rds@meta.data[rds@meta.data[[group_col]] %in% target,])
  rds@meta.data[[group_col]] <- as.character(rds@meta.data[[group_col]])
  return(subset(rds, cells = cell.use))
} 


mk.outdir <- function(dir) {
  if (!dir.exists(dir)){
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  }
  return(dir)
}

save_gg <- function(p, filename, width, height, format = c('pdf', 'png')) {
  for (d in format) {
    ggplot2::ggsave(filename = str_glue("{filename}.{d}"), plot = p, width = width, height = height, device = d, limitsize = FALSE)
  }
}
parse_comparisons <- function(input_str) {
  comparison_strs <- strsplit(input_str, ",")[[1]]
  parsed_list <- lapply(comparison_strs, function(x) {
    elements <- strsplit(trimws(x), "vs")[[1]]
    return(elements)
  })
  return(parsed_list)
}

addGeneExpressionStatusToMetadata <- function(seurat_obj, gene_name, expression_threshold = 0) {
  if (!(gene_name %in% rownames(seurat_obj@assays$RNA@data))) {
    stop("指定的基因不在 Seurat 对象的基因列表中！")
  }
  
  gene_expression <- seurat_obj@assays$RNA@data[gene_name, ]
  
  gene_name_use <- str_replace_all(gene_name, "-", "_")
  expression_status <- ifelse(gene_expression != expression_threshold, 
                              paste0(gene_name_use, "+"), 
                              paste0(gene_name_use, "-"))
  
  seurat_obj@meta.data[[paste0(gene_name_use, "_status")]] <- expression_status
  seurat_obj@meta.data[[paste0(gene_name_use, "_status")]] <- as_factor(seurat_obj@meta.data[[paste0(gene_name_use, "_status")]])
  
  return(list(obj = seurat_obj, gene_name = gene_name_use))
}

plot_gene_group_1 <- function(rds, genes, outdir, celltype_col, group_col, celltype_use = NULL, group_use = NULL, compare_str = NULL, set_w = NULL, set_h = NULL, cpal = 'Paired', gcpal = 'npg') {
  
  mk.outdir(outdir)
  
  # 点图 
  out.dot <- str_glue("{outdir}/gene_dot_exp")
  # 热图
  out.heatmap <- str_glue("{outdir}/gene_heatmap_exp")
  
  if (!is.null(celltype_use)){
    rds <- fetch.seurat(rds, celltype_col, celltype_use)
    rds@meta.data[[celltype_col]] <- factor(rds@meta.data[[celltype_col]], levels = celltype_use)

    if (length(celltype_use) == 1) {
      exp_method <- 'raw'
      slot <- 'data'
    }else {
      exp_method <- 'zscore'
      slot <- 'counts'
    }
  }else {
    exp_method <- 'zscore'
    slot <- 'counts'
  }
  
  if (!is.null(group_use)){
    rds <- fetch.seurat(rds, group_col, group_use)
  }
  
  all.cell <- rds@meta.data[[celltype_col]] %>% unique
  all.group <- rds@meta.data[[group_col]] %>% unique

  cpal <- get_color(all.cell, length(all.cell), palette = cpal)
  gcpal <- get_color(all.group, length(all.group), palette = gcpal)
  cpal.list <- list(cpal= cpal)
  
  ncell <- all.cell %>% length
  ngroup <-  all.group %>% length

  
  

  ngene <- length(genes)
  ht8 <- GroupHeatmap(rds,exp_method = exp_method, slot = slot,
                      features = genes, group.by = celltype_col, split.by = group_col, 
                      cluster_rows = FALSE, cluster_columns = FALSE, cluster_row_slices = FALSE, cluster_column_slices = FALSE,
                      add_dot = TRUE, add_reticle = TRUE, heatmap_palette = "viridis",cell_split_palcolor = gcpal, group_palcolor =  cpal.list, #group_palette = "Paired",# group_palcolor = cpal, #,group_palette = "Paired",
                      nlabel = 0, show_row_names = TRUE,show_column_names=F,
                      ht_params = list(row_gap = unit(0, "mm"), row_names_gp = gpar(fontsize = 10))) 
  
  if (ncell > 3) {
    w <- min(30, ncell*1.5)
    h <- min(30, ngene*0.9)
  }else {
    w <- 8
    h <- min(30, ngene*1.2)
  }

  if (!is.null(set_h)) {
    h <- as.double(set_h)
  } 
  if (!is.null(set_w)) {
    w <- as.double(set_w)
  } 
  
  print(w)
  print(h)
  
  save_gg(ht8$plot, out.dot, w, h)
  
  # 热图(1)
  ht1 <- GroupHeatmap(rds,
                      features = genes,
                      group.by = celltype_col,
                      exp_method = exp_method, slot = slot,
                      split.by = group_col,
                      show_row_names = TRUE,show_column_names=FALSE, cell_split_palcolor  = gcpal,group_palcolor =  cpal.list,cell_annotation_palcolor  = cpal,#group_palette = "Paired", ##group_palette = "Paired"# group_palcolor  = cpal,
  )
  print(cpal)
  save_gg(ht1$plot, out.heatmap , w, h*1.15)
  
  all_cells <- rds@meta.data[[celltype_col]] %>% unique
  # 小提琴图
  if (!is.null(compare_str)) {
    for (cell in all_cells) {
      rds.use <- fetch.seurat(rds, celltype_col, cell)
      w <- min(30, ngroup*1.5)
      for (gene in genes) {
        out.vio <- str_glue("{outdir}/{gene}_vio_exp")
        p <- FeatureStatPlot(
          rds.use , stat.by = gene, group.by = group_col, add_box = TRUE, stack = TRUE,add_trend = TRUE, palcolor = gcpal,
          comparisons = parse_comparisons(compare_str))
        save_gg(p, str_glue('{out.vio}_{cell}'), w, 5)
      }
    }
  }
  # gene-0-1
  if (!is.null(compare_str)) {
    for (cell in all_cells) {
      rds.use <- fetch.seurat(rds, celltype_col, cell)
      w <- min(30, ngroup*1.4)
      for (gene in genes) {
        out.status <- str_glue("{outdir}/{gene}_status_{cell}")
        rds.2.list <- addGeneExpressionStatusToMetadata(rds.use, gene, 0)
        rds.2 <- rds.2.list[['obj']]
        gene_name <- rds.2.list[['gene_name']]
        color_vec <- c('red', 'gray')
        names(color_vec) <- c(str_glue('{gene_name}+'), str_glue('{gene_name}-'))
        print(color_vec)
        p <- CellStatPlot(rds.2 , stat.by = str_glue(gene_name, '_status'), group.by = group_col, label = TRUE) + 
          scale_fill_manual(values = color_vec)
        save_gg(p, out.status , w, 3.5)
      }
    }
  }
}


parser <- ArgumentParser(description = "")

parser$add_argument("--rds", type = "character", required = TRUE, help = "")
parser$add_argument("--gene_list", type = "character", required = TRUE, help = "文件一列, 列名是gene")
parser$add_argument("--outdir", type = "character", required = TRUE, help = "")
parser$add_argument("--celltype_col", type = "character", required = TRUE, help = "")
parser$add_argument("--group_col", type = "character", required = TRUE, help = "")
parser$add_argument("--cell_use", type = "character", required = FALSE, help = "细胞类型逗号分隔, 细胞类型中有空格的用::代替空格", default = 'all')
parser$add_argument("--compare_str", type = "character", required = FALSE, help = "AvsB,AvsC,...", default = NULL)
parser$add_argument("--set_w", type = "character", required = FALSE, help = "None", default = NULL)
parser$add_argument("--set_h", type = "character", required = FALSE, help = "None", default = NULL)

parser$add_argument("--ccpal", type = "character", required = FALSE, help = "None", default = 'Paired')
parser$add_argument("--gcpal", type = "character", required = FALSE, help = "None", default = 'npg')

args <- parser$parse_args()


rds <- readRDS(args$rds)
rds_V4 <- Convert_Assay(seurat_object = rds, convert_to = "V3")
genes <- read_tsv(args$gene_list)$gene %>% unique()

if (args$cell_use != "all") {
  cell_filter <- str_split(args$cell_use, ",")[[1]]
  cell_filter <- sub("::", " ", cell_filter)
}else {
  cell_filter <- NULL
}

print(genes)
print(args$compare_str)

plot_gene_group_1(
  rds_V4, 
  genes = genes,
  outdir = args$outdir,
  celltype_col = args$celltype_col,
  group_col = args$group_col,
  celltype_use = cell_filter,
  compare_str = args$compare_str,
  set_h = args$set_h,
  set_w = args$set_w,
  cpal = args$ccpal, gcpal = args$gcpal
)


#rds <- readRDS("/public/home/zhliu/singlecell/project/PBMC_V2/02.celltype/v1/pbmc_v1_2.rds")

#rds_V4 <- Convert_Assay(seurat_object = rds, convert_to = "V3")

#cells <- c(
#'CD4_Treg_FOXP3'
#'CD4_TEM_ANXA1',
#'CD8_TEM_GNLY',
#'CD8_TEM_GZMK',
#'CD8_TEM_GNLY'
#)


#plot_gene_group_1(
#    rds_V4, 
#    genes = c("FOXP3","CTLA4","IL2RA","IL10RA","IKZF4","LAG3"),
#    outdir = "test_out/",
#    celltype_col = 'subcelltype',
#    group_col = 'diabetes',
#    celltype_use = cells,
#    compare_str = "GDMvsCtrl,T2DMvsCtrl,GDMvsT2DM"

#)
