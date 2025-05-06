library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(argparse)


parser <- ArgumentParser()
parser$add_argument("--rds", required=TRUE, help="Path to seurat RDS file")
parser$add_argument("--outdir", required=TRUE, help="")
parser$add_argument("--celltype_col", required=TRUE, help="celltype column name")
parser$add_argument("--genes", required=FALSE, help=",分隔")
parser$add_argument("--cells", required=FALSE, help=",分隔")
parser$add_argument("--gene_df", required=FALSE, help="genes celltype", default = 'None')
args <- parser$parse_args()


dir.create(args$outdir, recursive = T)

if (args$gene_df != 'None') {
    df <- read_tsv(args$gene_df)
    genes <- unique(df$genes)
    cells_sort <- unique(df$celltype)
}else {
    genes <- str_split(args$genes, ",")[[1]] %>% unique
    cells_sort <- str_split(args$cells, ",")[[1]] %>% unique
}

rds <- readRDS(args$rds)
Idents(rds) <- rds@meta.data[[args$celltype_col]]



p <- DotPlot(rds, features = genes) +theme(axis.text.x.top = element_text(angle = 90, vjust = 1,hjust=0) )+
  scale_y_discrete(position="right")+scale_x_discrete(position="top")+scale_colour_gradientn(colours = brewer.pal(11, name = "Reds"))+
  theme(axis.line = element_line(linewidth = 0)) + theme(panel.border = element_rect(color = "black", size = 0.5))



p$data <- p$data %>% mutate(id = as.character(id)) %>%  mutate(id = factor(id, levels= cells_sort))

h.base <- p$data$id %>% unique() %>% length()
w.base <- length(genes)
h <- max(h.base*0.3, 5)
w <- max(w.base*0.6, 12)

print(h)
print(w)

ggsave(filename = str_glue("{args$outdir}/{args$celltype_col}_marker.pdf"), plot = p,width = w, height = h)
ggsave(filename = str_glue("{args$outdir}/{args$celltype_col}_marker.png"), plot = p,width = w, height = h)
