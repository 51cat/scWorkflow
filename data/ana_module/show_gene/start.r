library(scCustomize)
library(SCP)
library(tidyverse)
library(argparse)
library(vroom)

print(sessionInfo())

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

to_v4 <- function(rds) {
    rds_V3 <- Convert_Assay(seurat_object = rds, convert_to = "V3")
    return(rds_V3)
}

plot_gene <- function(rds, gene, group_col, outdir, subdir = NULL, split_col =NULL, reduction_name = 'umap') {
    
    if (is.null(subdir) | subdir == 'None') {
        out.p <- mk.outdir(str_glue('{outdir}/single_gene/'))
    }else {
        out.p <- mk.outdir(str_glue('{outdir}/single_gene/{subdir}/'))
    }
    total.gr <- rds@meta.data[[group_col]] %>% unique %>% length
    p.UMAP <- FeatureDimPlot(
        srt = rds, features = gene, split.by = split_col, theme_use = "theme_blank", reduction = reduction_name
        )
    p.violin <- FeatureStatPlot(rds, stat.by = gene, group.by = group_col, add_box = TRUE)
    # save
    save_gg(p.UMAP, str_glue("{out.p}/{gene}_UMAP"), 4,4)
    save_gg(p.violin, str_glue("{out.p}/{gene}_{group_col}_violin"), max(total.gr*0.45, 7),4)
}

plot.all.gene <- function(rds, genes, group_by, outdir) {
    out.p <- mk.outdir(str_glue('{outdir}/all_gene/'))
    total.gr <- rds@meta.data[[group_by]] %>% unique %>% length
    if (is.list(genes)) {
        p2.dotplot <- DotPlot(rds, genes)+theme(axis.text.x = element_text(angle = 90, hjust = 1))
        out.name <- str_glue("{out.p}/dotplot_{group_by}_1")
    }
    

    save_gg(p2.dotplot, out.name, max(total.gr*0.5, 10),max(5, length(genes)*0.4))
    
}



parser <- ArgumentParser()
parser$add_argument("--rds", required=TRUE, help="Path to seurat RDS file")
parser$add_argument("--gene_list", required=TRUE, help="Path to gene list file, first column name must be genes")
parser$add_argument("--group_col", required=TRUE, help="Column for grouping")
parser$add_argument("--outdir", required=TRUE, help="Output directory")
parser$add_argument("--reduction_name", required=FALSE, default = 'umap')
args <- parser$parse_args()

rds_p <- args$rds
genes_file <- args$genes_file
group_cols <- str_split(args$group_cols, ",")[[1]]
outdir <- args$outdir
reduction_name <- args$reduction_name


rds <- readRDS(rds_p)
dir.create(outdir,recursive = T)

genes_df <- tryCatch({
    vroom(genes_file)
}, error = function(e) {
    read_tsv(genes_file)
})

genes_filter <- intersect(rownames(rds), genes_df$gene)
genes_df <- genes_df %>% filter(genes %in% genes_filter)

if (ncol(genes_df) == 2) {
    names(genes_df) <- c('gene', 'celltype')
}else {
    genes_df <- genes_df %>% mutate(celltype = 'None')
}

print(genes_df)
print(group_cols)

cells.all <- genes_df %>% pull(celltype) %>% unique()
print(cells.all)
rds <- to_v4(rds)

for (group_use in group_cols) {

    p.group <- CellDimPlot(rds, group.by = group_use, reduction = reduction_name,label = TRUE, label_insitu = TRUE)

    
    save_gg( p.group , str_glue('{outdir}/{group_use}_{reduction_name}'), 12,12)

    for (cell_use in cells.all) {
        genes_input <- genes_df %>% filter(celltype == cell_use) %>% pull(gene) %>% unique

        for (g in genes_input) {
            plot_gene(rds, g, group_use,  outdir = outdir, subdir = cell_use, reduction_name=reduction_name)
        }

    }
}

for (group_use in group_cols) {
    plot.all.gene(rds, unique(genes_filter), group_use, outdir)
}

print('Finish')
