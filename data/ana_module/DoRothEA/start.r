## We load the required packages
library(Seurat)
library(decoupleR)
library(argparse)
# Only needed for data handling and plotting
library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(ggplot2)
library(pheatmap)
library(tidyverse)
library(progress)
library(SCP)
library(rlang)


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

DB <- str_glue("{dirnameScript()}/TFnet/")
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


load.db <- function(spec = 'human') {
    if (spec == 'human') {
        net <- read_tsv(str_glue('{DB}/human_TF_decoupleR.tsv'))
        return(net)
    }
    if (spec == 'mouse') {
        net <- read_tsv(str_glue('{DB}/mouse_TF_decoupleR.tsv'))
        return(net)
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

analysis <- function(seurat_obj, celltype_col) {
    Seurat::Idents(seurat_obj) <- seurat_obj@meta.data[[celltype_col]]
    mat <- as.matrix(seurat_obj@assays$RNA$data)
    acts <- decoupleR::run_ulm(mat = mat, 
                           net = net, 
                           .source = 'source', 
                           .target = 'target',
                           .mor='mor', 
                           minsize = 5)
    
    # Extract ulm and store it in tfsulm in pbmc
    seurat_obj[['tfsulm']] <- acts %>%
                        tidyr::pivot_wider(id_cols = 'source', 
                                        names_from = 'condition',
                                        values_from = 'score') %>%
                        tibble::column_to_rownames('source') %>%
                        Seurat::CreateAssayObject(.)
    # Change assay
    DefaultAssay(object = seurat_obj) <- "tfsulm"
    seurat_obj <- Seurat::ScaleData(seurat_obj)
    seurat_obj@assays$tfsulm$data <- seurat_obj@assays$tfsulm$scale.data

    return(list(
        rds = seurat_obj,
        tf_table = acts
    ))
}


write_out_table <- function(seurat_obj, tf_table, db, outdir) {

    out.p <- mk.outdir(str_glue("{outdir}/tf_table/"))

    acts2 <- tf_table %>% rename(barcode = condition) 
    meta <- seurat_obj@meta.data %>% mutate(barcode = rownames(.))
    
    act.out <- acts2 %>% left_join(meta)
    reg <- net %>% filter(source %in% act.out$source) 
    # write
    #act.out %>% write_csv(str_glue("{out.p}/acts.csv"))
    reg %>% write_csv(str_glue("{out.p}/reg_gene.csv"))

}




# 计算平均表达量
cal_mean_exp <- function(seurat_obj, group_col, celltype_col) {
    mat <- AverageExpression(seurat_obj, group.by = c(celltype_col, group_col))$tfsulm
    
    #ann.df <- tibble(
    #        cluster = seurat_obj@meta.data[[celltype_col]], 
    #        group = seurat_obj@meta.data[[group_col]], 
    #        colname = str_c(cluster, "_", group)) %>% 
    #        distinct() %>% 
    #        column_to_rownames(var = 'colname')
    
    

    mat.long <- mat %>% 
            t() %>% 
            as_tibble() #%>% 
            #gather(key = "TF", value = "act", -c("cluster", "group"))
    
    return(list(
        mean_exp = mat,
        mean_exp_out = mat.long
    ))
}

cal_mean_exp2 <- function(seurat_obj, by) {
    mat <- AverageExpression(seurat_obj, group.by = by)$tfsulm
    mat.long <-  mat %>% 
        as.data.frame() %>% 
        mutate(gene = rownames(.)) %>% 
        as_tibble() %>%
        gather(key = "group", value = "act", -gene)
    

    return(list(
        mean_exp = mat,
        mean_exp_out = mat.long
    ))
}


# 分析差异转录因子
get_markers_all <- function(rds, compare_str, group.by, cluster_col, not_compare_cell = FALSE) {
  
  if (not_compare_cell == FALSE) {
    cells <- rds@meta.data[[cluster_col]] %>% unique()
  }else {
    cluster_col <- 'cell_col'
    rds@meta.data[['cell_col']] <- 'all'
    cells <- rds@meta.data[[cluster_col]] %>% unique()
  }


  cells <- rds@meta.data[[cluster_col]] %>% unique()
  print(cells)
  
  comparisons <- parse_comparisons(compare_str)
  inx <- 1
  df.list <- list()

  for (cell in cells) {
     for (s in comparisons) {
        cell.use <- rownames(rds@meta.data[rds@meta.data[[cluster_col]] %in% cell,])
        rds.sub <- subset(rds, cells = cell.use)

        # 检查分组是否存在, 细胞数是否小于3
        flag1 <- all(s %in% rds.sub@meta.data[[group.by]])
        flag2 <- all(
            rds.sub@meta.data[rds.sub@meta.data[[group.by]] %in% s,] %>% 
            group_by(!!sym(group.by)) %>% 
            summarise(total = n()) %>% pull(total) >= 3) 
        
        if (flag1 & flag2) {
            df <- FindMarkers(rds.sub, 
                    ident.1 = s[1], 
                    ident.2 = s[2], 
                    group.by = group.by,
                    only.pos =  F, 
                    verbose = T, 
                    min.pct= 0.25,
                    logfc.threshold = 0.25)%>% 
                    mutate(gene = rownames(.)) %>%
        mutate(cluster = cell, group = str_glue("{s[1]} vs {s[2]}"))
        df.list[[inx]] <- df
        inx <- inx + 1 
        print(str_glue("finish {cell} {s[1]}  {s[2]}"))
        }else {
            print(str_glue('skip {cell} {s[1]}  {s[2]}!'))
        }
     }
  }
  df.all <- do.call("rbind", df.list) 
  df.all <- df.all %>% filter(p_val_adj <= 0.05)
  return(df.all)
}


plot.all.TF <- function(rds, group_col, cluster_col, outdir) {
    colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
    colors.use <- grDevices::colorRampPalette(colors = colors)(100)
    
    mat <- AverageExpression(rds, group.by = c(cluster_col, group_col))$tfsulm
    ann.df <- tibble(
        cluster = rds@meta.data[[cluster_col]], 
        group = rds@meta.data[[group_col]], 
        colname = str_c(cluster, "_", group)) %>% 
        distinct() %>% 
        column_to_rownames(var = 'colname')
    
    # heatmap
    # Choose color palette
    colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
    colors.use <- grDevices::colorRampPalette(colors = colors)(100)
    
    my_breaks <- c(seq(-2, 0, length.out = ceiling(100 / 2) + 1),
                   seq(0.05, 2, length.out = floor(100 / 2)))
    if (cluster_col == 'seurat_clusters') {
        rownames(ann.df) <- str_c('g', rownames(ann.df))
    }

    print(ann.df)
    
    # Plot
    pheatmap::pheatmap(mat = mat %>% t() ,
                       color = colors.use,
                       border_color = NA,
                       breaks = my_breaks,
                       cluster_cols = T,
                       cluster_rows = F,
                       treeheight_row = 20,
                       treeheight_col = 20, main = 'TF_heatmap', annotation_row = ann.df,
                       filename = str_glue("{outdir}/TF_heatmap_all.pdf"), width = 40, fontsize_col = 4)
}

compare.all.TF.single <- function(seurat_obj, group_col, cluster_col, outdir, compare_str, group_cpal ='npg') {
    # box + violin by SCP
    DefaultAssay(object = seurat_obj) <- "tfsulm"

    total_TF <- nrow(seurat_obj)
    print(str_glue("Total TF: {total_TF}"))
    inx <- 0
    #setTxtProgressBar(pb, inx)
    cells.all <- seurat_obj@meta.data[[cluster_col]] %>% unique
    print(cells.all)
    for (cell.use in cells.all) {
        out.p <- mk.outdir(str_glue("{outdir}/all_TF_compare/{cell.use}/"))

        cells <- rownames(seurat_obj@meta.data[seurat_obj@meta.data[[cluster_col]] == cell.use,])

        rds.sub <- subset(seurat_obj, cells = cells)
        gene_sums <- Matrix::rowSums( rds.sub@assays$tfsulm@counts)
        target_TF_plot <- names(gene_sums[gene_sums > 0])

        for (TF in target_TF_plot) {
            p <- FeatureStatPlot(rds.sub, stat.by = TF, group.by = group_col, add_box = TRUE, stack = TRUE,add_trend = TRUE, comparisons = parse_comparisons(compare_str),  palette = group_cpal)
            save_gg(p, str_glue("{out.p}/{cell.use}_{TF}_compare"), 4, 4)
        }
    }
    
}




parser <- ArgumentParser(description = "")

parser$add_argument("--rds", help = "")
parser$add_argument("--tf_rds", help = "", default = 'None')
parser$add_argument("--outdir", help = "")
parser$add_argument("--celltype_col", help = "")
parser$add_argument("--group_col", help = "")
parser$add_argument("--compare_str", help = "")
parser$add_argument("--spec", help = "", default = "human")
parser$add_argument("--ntop_TF", help = "", default = '5')
parser$add_argument("--drop_celltype", help = "", default = 'None')
#parser$add_argument("--target_TF", help = "", default = 'None')
parser$add_argument("--show_all_TF", help = "", default = 'FALSE')
parser$add_argument("--compare_method", help = "", default = '1')

parser$add_argument("--gcpal", help = "", default = 'npg')
args <- parser$parse_args()

rds_path <- args$rds
tf_rds_path <- args$tf_rds
outdir <- args$outdir
celltype_col <- args$celltype_col
group_col <- args$group_col
spec <- args$spec
compare_str <- args$compare_str
target_TF <- args$target_TF
ntop_TF <- as.integer(args$ntop_TF)
compare_all_TF <- args$show_all_TF
compare_method <- str_split(args$compare_method, ",")[[1]]

group_cpal <- args$gcpal

if (tf_rds_path == 'None') {
    net <- load.db(spec)
    seurat_obj <- readRDS(rds_path)

    if (args$drop_celltype != 'None') {
        print("drop celltype!")
        print(args$drop_celltype)
        drop_celltype <- str_split(args$drop_celltype, ",")[[1]]
        cell.use <- rownames(seurat_obj@meta.data[!seurat_obj@meta.data[[celltype_col]] %in% drop_celltype,])
        rds.sub <- subset(seurat_obj, cells = cell.use)
        
        seurat_obj <- rds.sub
        print(seurat_obj@meta.data[[celltype_col]] %>% unique)
    }

    # 转录因子分析
    results <- analysis(seurat_obj, celltype_col)
    acts.df <- results[['tf_table']]
    seurat_obj_tf <- results[['rds']]

    # 输出表格
    write_out_table(seurat_obj_tf, acts.df, net, outdir)

    # 绘制所有转录因子
    saveRDS(seurat_obj_tf, str_glue("{outdir}/tf.rds"))
} else {
    seurat_obj_tf <- readRDS(tf_rds_path)
    print(DefaultAssay(seurat_obj_tf))
    out.p <- mk.outdir(str_glue("{outdir}"))
}

#plot.all.TF(seurat_obj_tf, group_col, celltype_col, outdir)




if ('1' %in% compare_method) {
    #  差异分析(组间+celltype)
    diff.TF.df <- get_markers_all(seurat_obj_tf, compare_str, group_col, celltype_col)
    diff.TF.df %>% write_csv(str_glue("{outdir}/TF_diff.csv"))

    # 平均表达量计算
    TF.mean <- cal_mean_exp2(seurat_obj_tf, c(group_col, celltype_col))
    mean.mat <- TF.mean$mean_exp
    TF.mean$mean_exp_out %>% write_csv(str_glue("{outdir}/mean_exp_all_TF.csv"))

    # top  热图
    #diff.TF.df$fc = diff.TF.df$pct.1 - diff.TF.df$pct.2
    top.TF <- diff.TF.df %>% 
        filter(p_val_adj <= 0.05) %>%
        group_by(cluster, group) %>% 
        top_n(ntop_TF, avg_log2FC) %>% pull(gene) %>% unique()


    mat.top <- mean.mat[top.TF, ]

    w <- min(dim(mat.top)[1]*0.1, 20)
    h <- min(dim(mat.top)[2]*0.60, 20)

    colors <- get_color(c(1:1000), n = 1000, palette = "gsea")
    pheatmap(mat.top %>% as.matrix, scale = "row", 
                color = colors, #colorRampPalette(c('blue','white','red'))(1000),
                cluster_cols = F, cluster_rows = T, border_color = NA,
                filename = str_glue("{outdir}/TF_heatmap_top_{ntop_TF}.pdf"), main = str_glue('Top {ntop_TF} TF'), height = h, width = w, treeheight_row = 0)
    
    pheatmap(mat.top %>% as.matrix, scale = "row", 
                color = colors, #colorRampPalette(c('blue','white','red'))(1000),
                cluster_cols = F, cluster_rows = T, border_color = NA,
                filename = str_glue("{outdir}/TF_heatmap_top_{ntop_TF}.png"), main = str_glue('Top {ntop_TF} TF'), height = h, width = w, treeheight_row = 0)
}




if ('2' %in% compare_method) {

    # 差异分析 (不分celltype)
    #diff.TF.df <- get_markers_all(seurat_obj_tf, compare_str, group_col, celltype_col, not_compare_cell = TRUE)
    #diff.TF.df %>% write_csv(str_glue("{outdir}/TF_diff_group.csv"))

    # 平均表达量计算
    TF.mean.group <- cal_mean_exp2(seurat_obj_tf, c(group_col))
    TF.mean.cell <- cal_mean_exp2(seurat_obj_tf, c(celltype_col))

    mean.mat.group <- TF.mean.group$mean_exp
    mean.mat.cell <- TF.mean.cell$mean_exp

    TF.mean.group$mean_exp_out %>% write_csv(str_glue("{outdir}/mean_exp_all_TF_group.csv"))
    TF.mean.cell$mean_exp_out %>% write_csv(str_glue("{outdir}/mean_exp_all_TF_cell.csv"))

    #TF.mean$mean_exp_out %>% write_csv(str_glue("{outdir}/mean_exp_all_TF_group.csv"))

    # top  热图
    #diff.TF.df$fc = diff.TF.df$pct.1 - diff.TF.df$pct.2
    #top.TF <- diff.TF.df %>% group_by(group) %>% top_n(ntop_TF, fc) %>% pull(gene) %>% unique()

    # gruop
    ## FindAllmarkers
    
    Idents(seurat_obj_tf) <- seurat_obj_tf@meta.data[[group_col]]
    
    diff.TF.df <- FindAllMarkers(seurat_obj_tf, only.pos =  F, 
                    verbose = T, 
                    min.pct= 0.25,
                    logfc.threshold = 0.25)


    top.TF.group <- diff.TF.df %>%
        filter(p_val_adj <= 0.05) %>% 
        group_by(cluster) %>% 
        top_n(ntop_TF, avg_log2FC) %>% pull(gene) %>% unique()

    print(top.TF.group)
    print(mean.mat.group)

    mat.top <- mean.mat.group[top.TF.group, ]
    diff.TF.df %>% write_csv(str_glue("{outdir}/TF_diff_group.csv"))



    colors <- get_color(c(1:1000), n = 1000, palette = "gsea")
    pheatmap(mat.top %>% as.matrix, scale = "row", 
                color =colors,
                cluster_cols = F, cluster_rows = T, border_color = NA,
                filename = str_glue("{outdir}/TF_heatmap_top_{ntop_TF}_group.pdf"), main = str_glue('Top {ntop_TF} TF'),height = 8, width = 3, treeheight_row = 0)

    pheatmap(mat.top %>% as.matrix, scale = "row", 
                color = colors, #colorRampPalette(c('blue','white','red'))(1000),
                cluster_cols = F, cluster_rows = T, border_color = NA,
                filename = str_glue("{outdir}/TF_heatmap_top_{ntop_TF}_group.png"), main = str_glue('Top {ntop_TF} TF'), height = 8, width = 3, treeheight_row = 0)
    
    # cell
    ## FindAllmarkers

    Idents(seurat_obj_tf) <- seurat_obj_tf@meta.data[[celltype_col]]
    
    diff.TF.df <- FindAllMarkers(seurat_obj_tf, only.pos =  F, 
                    verbose = T, 
                    min.pct= 0.25,
                    logfc.threshold = 0.25)
    diff.TF.df %>% write_csv(str_glue("{outdir}/TF_diff_cell.csv"))


    top.TF.cell <- diff.TF.df %>% 
        filter(p_val_adj <= 0.05) %>% 
        group_by(cluster) %>% 
        top_n(ntop_TF, avg_log2FC) %>% pull(gene) %>% unique()

    mat.top <- mean.mat.cell[top.TF.cell, ]
    

    mat.top <- mean.mat.cell[top.TF.cell, ]


    w <- min(dim(mat.top)[1]*0.15, 20)
    h <- min(dim(mat.top)[2]*1.5, 20)

    colors <- get_color(c(1:1000), n = 1000, palette = "gsea")
    pheatmap(mat.top %>% as.matrix, scale = "row", 
                color =colors,
                cluster_cols = F, cluster_rows = T, border_color = NA,
                filename = str_glue("{outdir}/TF_heatmap_top_{ntop_TF}_cell.pdf"), main = str_glue('Top {ntop_TF} TF'),height = h, width = w, treeheight_row = 0)

    pheatmap(mat.top %>% as.matrix, scale = "row", 
                color = colors, #colorRampPalette(c('blue','white','red'))(1000),
                cluster_cols = F, cluster_rows = T, border_color = NA,
                filename = str_glue("{outdir}/TF_heatmap_top_{ntop_TF}_cell.png"), main = str_glue('Top {ntop_TF} TF'), height = h, width = w, treeheight_row = 0)

}





# 展示所有转录因子的分析结果
if (compare_all_TF == 'True') {
    compare.all.TF.single(seurat_obj_tf, group_col, celltype_col, outdir, compare_str, group_cpal)
}

