library(Seurat)
library(tidyverse)
library(ggalluvial)
library(ggpubr) 
library(ggsci)
library(argparse)
library(SCP)
library(ggsci)

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

save_gg <- function(p, filename, width, height, format = c('pdf', 'png')) {
    for (d in format) {
        ggsave(filename = str_glue("{filename}.{d}"), plot = p, width = width, height = height, device = d)
    }
}

save_df <- function(df, filename, format = 'csv') {
    if (format == 'csv') {
        df %>% write_csv(filename)
    }else {
        df %>% write_tsv(filename)
    }
    
}

mk.outdir <- function(dir) {
    if (!dir.exists(dir)){
        dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    }
    return(dir)
}

show.scatter <- function(rds, group_col, celltype_col, outdir, reduction_name = 'umap', gcpal = 'npg', ccpal = 'Paired') {
    out.p <- mk.outdir(str_glue("{outdir}/{reduction_name}_plot/"))

    p.celltype.out <- str_glue("{out.p}/celltype_{reduction_name}") 
    p.celltype.out1 <- str_glue("{out.p}/celltype_{reduction_name}_no_legend") 
    p.celltype.out2 <- str_glue("{out.p}/celltype_{reduction_name}_2")
    p.celltype.out3 <- str_glue("{out.p}/celltype_{reduction_name}_no_legend_2")
    p.group.out <- str_glue("{out.p}/group_{reduction_name}")
    p.group.split.out <- str_glue("{out.p}/celltype_{reduction_name}_split")
    p.group.split.out2 <- str_glue("{out.p}/group_{reduction_name}_split")


    vec_cell <- rds@meta.data[[celltype_col]] %>% unique
    vec_gr <- rds@meta.data[[group_col]] %>% unique

    cell_cpal <- get_color(vec_cell, length(vec_cell), palette = ccpal)
    group_cpal <- get_color(vec_gr, length(vec_gr), palette = gcpal)



    p.celltype.1 <- CellDimPlot(
        rds, raster = FALSE,
        group.by = celltype_col, 
        palcolor = cell_cpal,
        label = T, 
        reduction = reduction_name, 
        theme_use = ggplot2::theme_classic, 
        theme_args = list(base_size = 10) )+ guides(color = guide_legend(ncol = 1))
    
    p.celltype.1.nolegend <- CellDimPlot(
        rds, raster = FALSE,
        group.by = celltype_col, 
        palcolor = cell_cpal,
        label = T, 
        reduction = reduction_name, 
        theme_use = ggplot2::theme_classic, 
        theme_args = list(base_size = 10),legend.position = "none" )

    p.celltype.2 <- CellDimPlot(
        rds, raster = FALSE,
        group.by = celltype_col, 
        palcolor = cell_cpal,
        reduction =reduction_name,
        label = TRUE, label_insitu = TRUE, label_repel = TRUE, label_segment_color = "red") + guides(color = guide_legend(ncol = 1))
    
    p.celltype.3 <- CellDimPlot(
        rds, raster = FALSE,
        group.by = celltype_col, 
        palcolor = cell_cpal,
        reduction =reduction_name,
        label = TRUE, label_insitu = TRUE, label_repel = TRUE, label_segment_color = "red", legend.position = "none")

    p.group <- CellDimPlot(
        rds, raster = FALSE,
        group.by = group_col, 
        palcolor = group_cpal,
        reduction = reduction_name, 
        theme_use = ggplot2::theme_classic, 
        theme_args = list(base_size = 10))

    ngroup <- rds@meta.data[[group_col]] %>% unique %>% length
    
    p.group.split <- CellDimPlot(
        rds, raster = FALSE,
        group.by = celltype_col, 
        reduction = reduction_name, 
        palcolor = cell_cpal,
        theme_use = ggplot2::theme_classic, 
        theme_args = list(base_size = 10), legend.position = "none", split.by = group_col)
    
    p.group.split2 <- CellDimPlot(
        rds, raster = FALSE,
        group.by = group_col, 
        palcolor = group_cpal,
        reduction = reduction_name, 
        theme_use = ggplot2::theme_classic, 
        theme_args = list(base_size = 10),
        split.by = group_col,legend.position = "none")




    save_gg(p.celltype.1, p.celltype.out, 7, 7.5)
    save_gg(p.celltype.1.nolegend, p.celltype.out1, 7, 7.5)
    save_gg(p.celltype.2, p.celltype.out2, 7, 7.5)
    save_gg(p.celltype.3, p.celltype.out3, 7, 7.5)
    save_gg(p.group, p.group.out, 7, 7.5)

    save_gg(p.group.split, p.group.split.out, min(30, 3.5*ngroup), 4)
    save_gg(p.group.split2, p.group.split.out2, min(30, 3.5*ngroup), 5)

}

get.cell.pct <- function(rds, group_col, celltype_col, outdir) {
    ncell <- rds@meta.data %>% 
        group_by_at(c(group_col, celltype_col)) %>% 
        summarise(ncell = n()) 
    
    total <- rds@meta.data %>% 
        group_by_at(group_col) %>% 
        summarise(total = n()) 
    cell.pct <- ncell %>% left_join(total, by = group_col) %>% mutate(pct = 100*ncell/total) 
    save_df(cell.pct, str_glue("{outdir}/cell_count_{group_col}.csv"))
    return(cell.pct)
}


plot.bar.group <- function(df, x, fill, group_sort = NULL, y = 'pct', gcpal = 'npg') {
    
    if (!is.null(group_sort)) {
        print('Change factor..')
        df[[fill]] <- factor(df[[fill]], levels = group_sort)
    }
    
    vec <- df[[fill]] %>% unique
    #pal3 <- SCP::palette_scp(vec, palette = "Paired", n = length(vec))
    group_cpal <- get_color(vec, length(vec), palette = gcpal)
    print(group_cpal)
    ggplot(df, aes_string(x = x, y = y, fill = fill, stratum=fill, alluvium=fill)) +
        geom_bar(stat = "identity", position = "stack", width = 0.6, color = 'black') +
        geom_flow(width = 0.6)+
        scale_fill_manual(values = group_cpal) +
        theme_bw()+
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
        labs(
            x = "",
            y = "Percentage(%)",
            fill = "")

}


plot_cellpct <- function(rds, pct_df, group_col, celltype_col, outdir, group_sort, gcpal = 'npg', ccpal = 'Paired') {
    out.p <- mk.outdir(str_glue("{outdir}/bar_plot/"))
    #p.celltype.group <- CellStatPlot(rds, stat.by = group_col, group.by = celltype_col, plot_type = "trend")+ guides(fill = guide_legend(ncol = 1))+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    
    df <- pct_df %>%
        group_by_at(celltype_col) %>%
        mutate(total = sum(pct)) %>%
        mutate(pct = pct / total * 100) %>%
        ungroup()
    p.celltype.group <- plot.bar.group(df, celltype_col, group_col, group_sort = group_sort, gcpal= gcpal)+ guides(fill = guide_legend(ncol = 1))

    p.group.celltype <- CellStatPlot(rds, stat.by = celltype_col, group.by = group_col, plot_type = "trend", palette = ccpal)+ guides(fill = guide_legend(ncol = 1))+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

    p.celltype.group.out <- str_glue("{out.p}/x_celltype_y_{group_col}")
    p.group.celltype.out <- str_glue("{out.p}/x_{group_col}_y_celltype")

    ncell <- rds@meta.data[[celltype_col]] %>% unique %>% length()
    ngroup <- rds@meta.data[[group_col]] %>% unique %>% length()

    w <- min(ncell*0.8, 20)
    h <- min(ngroup*1.2, 20)


    save_gg(p.celltype.group, p.celltype.group.out, w,h)


    w <- min(ngroup*1.5, 35)
    h <- min(ncell*0.4, 20)


    save_gg(p.group.celltype, p.group.celltype.out, w, h)


}


parser <- ArgumentParser(description = "")

parser$add_argument("--rds", type = "character", required = TRUE, help = "")
parser$add_argument("--outdir", type = "character", required = TRUE, help = "")
parser$add_argument("--reduction_name", type = "character", required = FALSE, help = "", default = 'umap')

parser$add_argument("--group_col", type = "character", required = TRUE, help = "")
parser$add_argument("--celltype_col", type = "character", required = TRUE, help = "")
parser$add_argument("--show", type = "character", required = FALSE, help = "", default = 'scatter,bar')


parser$add_argument("--gcpal", type = "character", required = FALSE, help = "", default = 'npg')
parser$add_argument("--ccpal", type = "character", required = FALSE, help = "", default = 'Paired')
parser$add_argument("--group_sort", type = "character", required = FALSE, help = "", default = 'None')

args <- parser$parse_args()

rds.path <-  args$rds
outdir <-  args$outdir
group_col <-  args$group_col
cluster_col <-  args$celltype_col
reduction_name <-  args$reduction_name

if (args$group_sort != 'None') {
    group_sort <- str_split(args$group_sort, ",")[[1]]
}else {
    group_sort <- NULL
}

rds<- readRDS(rds.path)
print("load RDS success")

show.type <- str_split(args$show, ",")[[1]]

if ('scatter' %in% show.type) {
    show.scatter(rds, group_col, cluster_col, outdir, reduction_name, ccpal = args$ccpal,  gcpal = args$gcpal)
}

if ('bar' %in% show.type) {
    pct_df <- get.cell.pct(rds, group_col, cluster_col, outdir)
    plot_cellpct(rds, pct_df, group_col, cluster_col, outdir, group_sort = group_sort, ccpal = args$ccpal,  gcpal = args$gcpal)
}