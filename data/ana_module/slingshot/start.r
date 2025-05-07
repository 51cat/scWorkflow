library(tidyverse)
library(Seurat)
library(SCP)
library(rlang)
library(argparse)

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


process_data <- function(rds, celltype_col, reduction_name = 'umap') {
    rds <- RunSlingshot(srt = rds,group.by = celltype_col, reduction = reduction_name, show_plot = FALSE)
    return(rds)
}


show_traj <- function(rds, celltype_col, group_col, outdir, reduction_name = 'umap') {
    out.p <- mk.outdir(str_glue('{outdir}/trace_show/'))
    all_Lineage <- grep("^Lineage", colnames(rds@meta.data), value = TRUE)
    ngroup <- rds@meta.data[[group_col]] %>% unique %>% length
    print(all_Lineage)
    # save df
    for (L in all_Lineage) {
        rds@meta.data %>% select(!!sym(L)) %>% mutate(barcode = rownames(.)) %>% as_tibble() %>% write_csv(str_glue("{outdir}/{L}_time.csv"))
    }
    
     
    p1 <- CellDimPlot(rds, group.by = celltype_col, reduction = reduction_name, lineages = all_Lineage)
    save_gg(p1, str_glue("{out.p}/{celltype_col}_all_Lineage"), 8, 5)
    p2 <- CellDimPlot(rds, group.by = group_col, reduction = reduction_name, lineages = all_Lineage)
    save_gg(p2, str_glue("{out.p}/{group_col}_all_Lineage"), 8, 5)

    # 单独展示
    for (Lineage in all_Lineage) {
        p <- CellDimPlot(rds, group.by = celltype_col, reduction = reduction_name, lineages = Lineage)
        save_gg(p, str_glue("{out.p}/{celltype_col}_{Lineage}"), 8, 5)
        p <- CellDimPlot(rds, group.by = group_col, reduction = reduction_name, lineages = Lineage)
        save_gg(p, str_glue("{out.p}/{group_col}_{Lineage}"), 8, 5)
    }

    # 分样本
    p3 <- CellDimPlot(rds, group.by = celltype_col, reduction = reduction_name, lineages = all_Lineage, split.by = group_col) + guides(color = guide_legend(ncol = 1))
    
    w <- min(30, ngroup*8)
    save_gg(p2, str_glue("{out.p}/{group_col}_split_all_Lineage"), w, 5)


    for (Lineage in all_Lineage) {
        p <- CellDimPlot(rds, group.by = celltype_col, reduction = reduction_name, lineages = Lineage, split.by = group_col)
        w <- min(30, ngroup*8)
        save_gg(p, str_glue("{out.p}/{group_col}_split_{Lineage}"), w, 5)
    }

    # 时间展示
    for (Lineage in all_Lineage) {
        p1 <- FeatureDimPlot(rds, features = Lineage, reduction = reduction_name, theme_use = "theme_blank")
        save_gg(p1, str_glue("{out.p}/{Lineage}_time"), 8, 5)
        w <- min(30, ngroup*8)
        p2 <- FeatureDimPlot(rds, features = Lineage, reduction = reduction_name, theme_use = "theme_blank", split.by = group_col)
        save_gg(p2, str_glue("{out.p}/{Lineage}_split_time"), w, 5)
    }

}


parser <- ArgumentParser(description = "")

parser$add_argument("--rds", type = "character", required = TRUE, help = "")
parser$add_argument("--outdir", type = "character", required = TRUE, help = "")
parser$add_argument("--reduction_name", type = "character", required = FALSE, help = "", default = 'umap')

parser$add_argument("--group_col", type = "character", required = TRUE, help = "")
parser$add_argument("--celltype_col", type = "character", required = TRUE, help = "")


args <- parser$parse_args()

rds.path <-  args$rds
outdir <-  args$outdir
group_col <-  args$group_col
celltype_col <-  args$celltype_col
reduction_name <-  args$reduction_name

rds<- readRDS(rds.path)
print("load RDS success")

rds <- process_data(rds, celltype_col, reduction_name)
show_traj(rds, celltype_col, group_col, outdir, reduction_name)
print('finish')