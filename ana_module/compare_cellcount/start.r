library(Seurat)
library(tidyverse)
library(ggpubr)
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


save_gg <- function(p, filename, width, height, format = c('pdf', 'png')) {
    for (d in format) {
        ggsave(filename = str_glue("{filename}.{d}"), plot = p, width = width, height = height, device = d)
    }
}

mk.outdir <- function(dir) {
    if (!dir.exists(dir)){
        dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    }
    return(dir)
}



get.cell.pct <- function(rds, groups, celltype_col) {
    ncell <- rds@meta.data %>% 
        group_by_at(c(groups, celltype_col)) %>% 
        summarise(ncell = n()) 
    
    total <- rds@meta.data %>% 
        group_by_at(groups) %>% 
        summarise(total = n()) 
    cell.pct <- ncell %>% left_join(total, by = groups) %>% mutate(pct = 100*ncell/total) 
    return(cell.pct)
}


parse_comparisons <- function(input_str) {
  comparison_strs <- strsplit(input_str, ",")[[1]]
  parsed_list <- lapply(comparison_strs, function(x) {
    elements <- strsplit(trimws(x), "vs")[[1]]
    return(elements)
  })
  return(parsed_list)
}


compare.cell.pct <- function(rds, compare_str, rep_col, tre_col, celltype_col, outdir, gcpal = 'npg', method = "wilcox.test") {
    out.p <- mk.outdir(outdir)

    combinations  <- parse_comparisons(compare_str)

    print(combinations)

    cells <- rds@meta.data[[celltype_col]] %>% unique

    pct <- get.cell.pct(rds, c(rep_col, tre_col), celltype_col)
    
    groups <- rds@meta.data[[tre_col]] %>% unique

    group_cpal <- get_color(groups, length(groups), palette = gcpal)

    for (cell in cells) {
        use <- pct %>% filter(get(celltype_col) == cell)
        p <- ggplot(use, aes_string(tre_col, 'pct', fill = tre_col)) + 
            geom_boxplot(width = 0.5) + 
            geom_jitter(size = 0.1) + 
            stat_compare_means(comparisons=combinations,
                method=method,label = 'p.signif', size = 6, vjust = 0.5, hide.ns = T)+ scale_fill_manual(values=group_cpal ) +
            ggtitle(cell) + theme_classic()
        save_gg(p, str_glue("{out.p}/{cell}_box"), 3.5, 3.5)
        #ggsave(filename = str_glue("{cell}_box.png"), plot = p, width = 4, height = 4, device = "png")
    }
    pct %>% write_csv(str_glue("{out.p}/cell_pct.csv"))
}

parser <- ArgumentParser(description = "")

parser$add_argument("--rds", help = "")
parser$add_argument("--compare_str", help = "")
parser$add_argument("--outdir", help = "")

parser$add_argument("--tre_col", help = "") # 处理列
parser$add_argument("--rep_col", help = "", default = 'orig.ident') # 重复列

parser$add_argument("--celltype_col", help = "") # 细胞类型列

parser$add_argument("--compare_method", help = "wilcox.test or kruskal.test", default = 'wilcox.test')

parser$add_argument("--gcpal", help = "", default = 'npg')
args <- parser$parse_args()


compare.cell.pct(
    readRDS(args$rds),
    args$compare_str,
    args$rep_col,
    args$tre_col,
    args$celltype_col,
    args$outdir,
    args$gcpal,
    args$compare_method
)