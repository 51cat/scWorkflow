library(SCP)
library(Seurat)
library(tidyverse)
library(argparse)
library(ggpubr)
library(qs)


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


get.markers <- function(obj, method, ident.use, ident.1, ident.2 = NULL) {
  Idents(obj) <- obj@meta.data[[ident.use]]
  df <- FindMarkers(
    object = obj,
    assay = method,
    slot = "data",
    ident.1 = ident.1,
    ident.2 = ident.2,
    test.use = "wilcox",
    min.pct = -Inf,
    logfc.threshold = 0,
    min.cells.group = 0,
    min.diff.pct = -Inf,
    verbose = T,
    min.cells.feature = 0) 
  
  df <- df %>% 
    tibble::rownames_to_column(var = "gene") %>%
    dplyr::mutate(direction = dplyr::if_else(avg_log2FC >0, "up", "down")) %>%
    dplyr::select(-c("pct.1", "pct.2"))
  return(df)
}

compare_within <- function(obj, outdir) {
  prfx <- obj[['prfx']]
  out.p <- mk.outdir(str_glue("{outdir}/{prfx}/compare_within/"))
  seurat_obj <- obj[['obj']]
  compare_col <- obj[['compare_col']]
  group_col <- obj[['group_col']]
  methods <- obj[['methods_use']]
  
  all.group <- seurat_obj@meta.data[[group_col]] %>% unique
  all.celltype <- seurat_obj@meta.data[[compare_col]] %>% unique
  inx <- 1
  diff.df <- list()
  
  for (method.use in methods) {
    DefaultAssay(seurat_obj) <- method.use
    print(str_glue('start compare {method.use}'))
    
    for (gr.use in all.group) {
      cells <- rownames(seurat_obj@meta.data[seurat_obj@meta.data[[group_col]] == gr.use, ])
      rds.use <- subset(seurat_obj, cells = cells)
      
      # 对每个cell分析
      for (cell.use in all.celltype) {
        marker_df <- get.markers(rds.use, method.use, compare_col, ident.1 = cell.use, ident.2 = NULL)
        marker_df <- marker_df %>% mutate(cluster = cell.use, compare_method = 'within', method = method.use, group = gr.use)
        diff.df[[inx]] <- marker_df
        inx <- inx + 1
      }
    }
    print(str_glue('finish! {method.use}'))
  }
  
  df.out <- do.call('rbind',  diff.df)
  df.out %>% write_csv(str_glue("{out.p}/treatment_within_compare.csv"))
}

compare_between <- function(obj, outdir) {
  out.p <- mk.outdir(str_glue("{outdir}/compare_between/"))
  seurat_obj <- obj[['obj']]
  compare_col <- obj[['compare_col']]
  use_cell <- obj[['use_cell']]
  compare_list  <- obj[['compare_list']]
  group_col <- obj[['group_col']]
  methods <- obj[['methods_use']]
  inx <- 1
  diff.df <- list()
  
  for (method.use in methods) {
    print(str_glue('start compare {method.use} - {use_cell}'))
    for (s in compare_list) {
      marker_df <- get.markers(seurat_obj, method.use, compare_col, ident.1 = s[1], ident.2 = s[2])
      marker_df <- marker_df %>% mutate(cluster = use_cell, compare_method = 'between', method = method.use, group = str_glue("{s[1]}_vs_{s[2]}"))
      diff.df[[inx]] <- marker_df
      inx <- inx + 1
    }
    
  }
  df.out <- do.call('rbind',  diff.df)
  df.out %>% write_csv(str_glue("{out.p}/{use_cell}_treatment_between_compare.csv"))
  print(str_glue('finish compare {method.use} - {use_cell} !!'))
}

load_input_data <- function(dir, group_sort = NULL) {
  files <- list.files(path = str_glue("{dir}/score_obj/"), pattern = ".*_score\\.qs$", full.names = TRUE)
  print(files)
  input_data <- list()
  if (is.null(group_sort)) {
    for (inx in seq_along(files)) {
      input_data[[inx]] <-  qread(files[[inx]])
    }
  }else {
    for (inx in seq_along(files)) {
      
      print('change factor!')
      rds.input <-  qread(files[[inx]])
      obj <- rds.input[['obj']]
      group_col <- rds.input[['group_col']]
      obj@meta.data[[group_col]] <- factor(obj@meta.data[[group_col]], levels = group_sort)
      rds.input[['obj']] <- obj
      input_data[[inx]] <- rds.input
    }
    
  }
  print('load data success!')
  return(input_data)
}

plot_compare_within <- function(within_obj, outdir, target_pathway_show, ccpal = 'Paired') {
  #out.p <- mk.outdir(str_glue("{outdir}/compare_within/"))
  seurat_obj <- within_obj[['obj']]
  compare_col <- within_obj[['compare_col']]
  group_col <- within_obj[['group_col']]
  methods <- within_obj[['methods_use']]
  
  all.group <- seurat_obj@meta.data[[group_col]] %>% unique
  all.celltype <- seurat_obj@meta.data[[compare_col]] %>% unique
  
  cell_cpal <- get_color(all.celltype, length(all.celltype), palette = ccpal)

  ngroup <- length(all.group)
  ncell <- length(all.celltype)


  inx <- 1
  diff.df <- list()
  
  for (method.use in methods) {
    DefaultAssay(seurat_obj) <- method.use
    for (gr.use in all.group) {
      out.p <- mk.outdir(str_glue("{outdir}/compare_within/{method.use}/{gr.use}/"))
      cells <- rownames(seurat_obj@meta.data[seurat_obj@meta.data[[group_col]] == gr.use, ])
      rds.use <- subset(seurat_obj, cells = cells)
      DefaultAssay(rds.use) <- method.use
      
      target_pathway_input <- intersect(target_pathway_show, rownames(rds.use))
      
      for (pathway_show in target_pathway_input) {
        p1 <- FeatureStatPlot(
          rds.use, 
          palcolor = cell_cpal,
          stat.by = c(pathway_show), 
          group.by =compare_col, 
          title = gr.use, add_box = T, add_trend= T
        ) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
        # FeatureStatPlot(pancreas_sub, stat.by = c("G2M_score", "Fev"), group.by = "SubCellType", plot_type = "box")
        p2 <- FeatureStatPlot(
          rds.use, 
          stat.by = c(pathway_show), 
          group.by =compare_col, 
          palcolor = cell_cpal,
          plot_type = "box",
          title = gr.use, add_trend= T
        ) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
        save_gg(p1, str_glue('{out.p}/{pathway_show}_{method.use}_{gr.use}_violin'), ncell*0.6, 4)
        save_gg(p2, str_glue('{out.p}/{pathway_show}_{method.use}_{gr.use}_box'), ncell*0.6, 4)
      }
    }
  }
}

plot_compare_between <- function(between_obj_list, target_pathway, method, gcpal = 'npg') {
  p.list.vio <- list()
  p.list.box <- list()
  
  inx <- 1
  
  for (between_obj in between_obj_list) {
    seurat_obj <- between_obj[['obj']]
    use_cell <- between_obj[['use_cell']]
    compare_list  <- between_obj[['compare_list']]
    group_col <- between_obj[['group_col']]
    methods <- between_obj[['methods_use']]

    all.group <- seurat_obj@meta.data[[group_col]] %>% unique
    ngroup <- length(all.group)

    gr_cpal <- get_color(all.group, length(all.group), palette = gcpal)

    w <- 0.5*ngroup
    h <- 4

    if (w <= 3) {
      w <- 3
    }

    
    DefaultAssay(seurat_obj) <- method
    
    target_pathway_input <- intersect(target_pathway, rownames(seurat_obj))
    
    p.vionlin.1 <- FeatureStatPlot(
      seurat_obj, 
      palcolor = gr_cpal,
      stat.by = target_pathway, 
      group.by = group_col, comparisons = compare_list, 
      title = use_cell, add_box = T, add_trend=T ) + theme(legend.position = 'none')
    
    p.box.1 <- FeatureStatPlot(
      seurat_obj, 
      stat.by = target_pathway, 
      palcolor = gr_cpal,
      group.by = group_col, comparisons = compare_list, 
      title = use_cell, plot_type = "box",add_trend=T ) + theme(legend.position = 'none')
    
    p.list.vio[[use_cell]] <- p.vionlin.1
    p.list.box[[use_cell]] <- p.box.1
  }
  
  
  return(list(box=p.list.box, vio = p.list.vio, w = w, h = h))
}

run.diff.anlysis <- function(input_data, outdir) {
  for (input_use in input_data) {
    if (input_use[['compare_method']] == 'within') {
      compare_within(input_use, outdir)
      
    }
    if (input_use[['compare_method']] == 'between') {
      compare_between(input_use, outdir)
    }
  }
}


parser <- ArgumentParser(description = "R script with command line arguments for irGSEA")

parser$add_argument("--score_step_outdir", type = "character", default = "",
                    help = "Path to the score outdir")

parser$add_argument("--outdir", type = "character", default = "",
                    help = "")
parser$add_argument("--group_sort", type = "character", default = "None",
                    help = "group sort")

parser$add_argument("--gcpal", type = "character", default = "npg",
                    help = "group sort")

parser$add_argument("--ccpal", type = "character", default = "Paired",
                    help = "group sort")

args <- parser$parse_args()

score_step_outdir <- args$score_step_outdir
outdir <- args$outdir

gcpal <- args$gcpal
ccpal <- args$ccpal

print(score_step_outdir)
print(outdir)

# 分析差异通路
if (args$group_sort != 'None') {
  g_sort <- str_split(args$group_sort, ',')[[1]]
  input_data <- load_input_data(score_step_outdir, group_sort = g_sort)
}else {
  input_data <- load_input_data(score_step_outdir)
}

run.diff.anlysis(input_data, outdir)

# 绘制通路(总通路数大于50, 绘制每组的top5, 否则绘制全部)

# within

obj_within <- Filter(function(x) x$compare_method == "within", input_data)

if (length(obj_within) != 0) {
  
  for (obj_within_use in obj_within) {
    
    #obj_within <- obj_within[[1]]
    # 获取通路数量
    print(obj_within_use)
    methos_test <- obj_within_use[['methods_use']][1]
    total_pathway <- rownames(obj_within_use[['obj']][[methos_test]])
    
    if (length(total_pathway) <= 50) {
      target_pathway_show <- total_pathway
      print(target_pathway_show)
      prfx <- obj_within_use[['prfx']]
      outdir.within <- str_glue("{outdir}/{prfx}/")
      plot_compare_within(obj_within_use, outdir.within, target_pathway_show, ccpal =ccpal)
    }else {
      prfx <- obj_within_use[['prfx']]
      outdir.within <- str_glue("{outdir}/{prfx}/")
      diff.df <- read_csv(str_glue("{outdir.within}/compare_within/treatment_within_compare.csv"))
      
      target_pathway_show <- diff.df %>% filter(p_val <= 0.05) %>%
        group_by(cluster, method, group) %>%  
        arrange(p_val) %>%  
        slice_head(n = 5) %>% 
        pull(gene) %>% unique()
      print(target_pathway_show)
      plot_compare_within(obj_within_use, outdir.within, target_pathway_show, ccpal =ccpal)
    }
  }
}

# between
obj_between <- Filter(function(x) x$compare_method != "within", input_data)

if (length(obj_between) != 0) {
  methods <- lapply(obj_between, function(x) x[['methods_use']]) %>% unlist() %>% unique()
  cells <- lapply(obj_between, function(x) x[['use_cell']]) %>% unlist() %>% unique()
  
  # 获取通路数量
  methos_test <- obj_between[[1]][['methods_use']][1]
  total_pathway <- rownames(obj_between[[1]][['obj']][[methos_test]])
  
  if (length(total_pathway) <= 50) {
    for (method.use in methods) {
      out.p <- mk.outdir(str_glue("{outdir}/compare_between/{method.use}/"))
      
      for (target_pathway in total_pathway) {
        p.out <- plot_compare_between(obj_between, target_pathway, method.use, gcpal = gcpal)
        # save 
        p.box.all <- p.out[['box']]
        p.vio.all <- p.out[['vio']]
        
        for (cell in names(p.box.all)) {
          save_gg(p.box.all[[cell]], str_glue('{out.p}/{target_pathway}_{method.use}_{cell}_box'), p.out[['w']], p.out[['h']])
          save_gg(p.vio.all[[cell]], str_glue('{out.p}/{target_pathway}_{method.use}_{cell}_vio'), p.out[['w']], p.out[['h']])          
        }
      }
    }
  }
  
  else {
    files <- list.files(path = str_glue("{outdir}/compare_between/"), pattern = "*_between_compare.csv", full.names = TRUE)
    data_list <- lapply(files, read_csv)
    diff.df <- bind_rows(data_list)
    
    total_pathway <- diff.df %>% filter(p_val <= 0.05) %>%
      group_by(cluster, method, group) %>%  
      arrange(p_val) %>%  
      slice_head(n = 5) %>% 
      pull(gene) %>% unique()
    print(total_pathway)
    
    for (method.use in methods) {
      out.p <- mk.outdir(str_glue("{outdir}/compare_between/{method.use}/"))
      for (target_pathway in total_pathway) {
        p.out <- plot_compare_between(obj_between, target_pathway, method.use, gcpal = gcpal)
        # save 
        p.box.all <- p.out[['box']]
        p.vio.all <- p.out[['vio']]

        for (cell in names(p.box.all)) {
          save_gg(p.box.all[[cell]], str_glue('{out.p}/{target_pathway}_{method.use}_{cell}_box'), p.out[['w']], p.out[['h']])
          save_gg(p.vio.all[[cell]], str_glue('{out.p}/{target_pathway}_{method.use}_{cell}_vio'), p.out[['w']], p.out[['h']])          
        }
      }       
    }
  }
}
print('Finish!')