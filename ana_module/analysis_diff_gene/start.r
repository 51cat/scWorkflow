library(tidyverse)
library(Seurat)
library(argparse)
library(scRNAtoolVis)
library(ggrepel)

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
}

plot.vol <- function(data, p_value_cutoff = 0.05, log2FC_cutoff = 0.25) {
  data.use <- data %>% mutate(
    stat = case_when(
      p_val_adj > p_value_cutoff | abs(avg_log2FC) < log2FC_cutoff ~ "Not_significant",
      p_val_adj < p_value_cutoff & avg_log2FC >= log2FC_cutoff ~ "up",
      p_val_adj < p_value_cutoff & avg_log2FC <= -log2FC_cutoff ~ "down",
      
    )
  )
  label.df <- data.use %>% 
    filter(stat %in% c('up', 'down')) %>% 
    group_by(stat) %>% 
    top_n(5, wt = abs(avg_log2FC))
  
  # plot
  
  p <- ggplot(data.use, aes(x = avg_log2FC, y = -log10(p_val_adj), color = stat)) +
    geom_point() +
    scale_color_manual(values = c("up" = "red", "down" = "blue", "Not_significant" = "gray")) +
    labs(x = "Average log2 Fold Change", y = "-log10(P-value)") +
    geom_vline(xintercept = log2FC_cutoff, color = 'gray40')+
    geom_vline(xintercept = -log2FC_cutoff, color = 'gray40')+
    geom_hline(yintercept = -log10(p_value_cutoff), color = 'gray40') + 
    theme_bw() +
    theme(legend.title = element_blank())
  p <- p+ geom_text_repel(data = label.df, aes(label = gene), size = 3, color = "black")
  return(p)
}


plot.diff.vol.line <- function(DEGs, prfx, logfc_threshold = 0.25, pval_threshold = 0.05) {
  
  DEGs$difference <- DEGs$pct.1 - DEGs$pct.2
  
  DEGs_sig <- DEGs[DEGs$p_val_adj < pval_threshold & abs(DEGs$avg_log2FC) > logfc_threshold, ]
  DEGs_sig$label <- rownames(DEGs_sig)
  
  p<- ggplot(DEGs, aes(x = difference, y = avg_log2FC)) + 
    geom_point(size = 0.5, color = "grey60") + 
    geom_text_repel(data = DEGs_sig, aes(label = label), 
                    color = "black", fontface = "italic", size = 3) +
    geom_point(data = DEGs[DEGs$p_val_adj < pval_threshold & DEGs$avg_log2FC > logfc_threshold, ],
               aes(x = difference, y = avg_log2FC),
               size = 0.5, color = "red") +  
    geom_point(data = DEGs[DEGs$p_val_adj < pval_threshold & DEGs$avg_log2FC < -logfc_threshold, ],
               aes(x = difference, y = avg_log2FC),
               size = 0.5, color = "blue") +
    labs(x = "Delta Percent", 
         y = "Log-fold Change", 
         title = paste0(prfx, "-case rich vs control poor")) +
    theme_classic() +
    theme(axis.text.x = element_text(colour = 'black', size = 10),
          axis.text.y = element_text(colour = 'black', size = 10),
          axis.title = element_text(colour = 'black', size = 10),
          axis.line = element_line(color = 'black', size = 0.4),
          plot.title = element_text(face = "bold", size = 10)) +
    geom_hline(yintercept = 0, lty = 2, lwd = 0.4) +
    geom_vline(xintercept = 0, lty = 2, lwd = 0.4)
  
  return(p)
  
}

get_markers <- function(rds, 
                        target_cell, 
                        celltype_col,
                        g1, 
                        g2, 
                        group.by, 
                        only.p = F,
                        logfc.threshold = 0,
                        min.pct = 0.25
) {
  cell.use <- rownames(rds@meta.data[rds@meta.data[[celltype_col]] %in% target_cell,])
  rds.sub <- subset(rds, cells = cell.use)
  df <- FindMarkers(rds.sub, 
                    ident.1 = g1, 
                    ident.2 = g2, 
                    group.by = group.by,
                    only.pos = only.p, 
                    verbose = T, 
                    min.pct= min.pct,
                    logfc.threshold = logfc.threshold)%>% 
    mutate(gene = rownames(.)) %>%
    mutate(cluster = target_cell, group = str_glue("{g1} vs {g2}"))
  return(df)
}

plot_ndegs <- function(df, group_use, rds, celltype_col, fc=0.25, p_th=0.05) {
  df <- df %>% 
   # filter(group == group_use) %>% 
    filter(p_val_adj <= p_th) %>% 
    filter(abs(avg_log2FC) >= fc) 
  
  cellpct.df <- get.cell.pct(rds, celltype_col)
  

  
  ndegs <- df %>% group_by(group, cluster) %>% summarise(nGenes = n()) %>% ungroup()
  ndegs <- ndegs %>% expand(group, cluster) %>% full_join(ndegs)
  ndegs[is.na(ndegs)] <- 0
  ndegs <- ndegs %>% filter(group == group_use)
  ndegs$cluster <- factor(ndegs$cluster, levels = unique(ndegs$cluster[order(-ndegs$nGenes)]))
  ndegs <- ndegs %>% left_join(cellpct.df)
  
  p <- ggplot(ndegs) + geom_bar(aes(cluster, nGenes), stat = 'identity') +  
    ggtitle(str_glue('{group_use}: Number of DEGs'), subtitle = str_glue('Foldchange >= {fc}; p.adj <= {p_th}'))+theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
  
  df.label <- ndegs %>% filter(cell_pct >= mean(ndegs$cell_pct)) %>% filter(nGenes >= mean(ndegs$nGenes))
  p2 <- ggplot(ndegs, aes(nGenes, cell_pct)) + 
    geom_point(aes(size = nGenes, fill = cell_pct),shape = 21) + 
    scale_fill_viridis_c() + 
    ggrepel::geom_label_repel(df.label, mapping = aes(nGenes, cell_pct, label = cluster), 
                              max.overlaps = Inf, min.segment.length = 0, size = 3 ) + 
    geom_hline(yintercept = mean(ndegs$cell_pct), color = 'gray40', lty='dashed')+
    geom_vline(xintercept = mean(ndegs$nGenes),color = 'gray40', lty='dashed')+
    theme_classic()
    
  return(list(table =ndegs, plot = p, plot.cellpct = p2 ))
}

get.cell.pct <- function(rds, celltype_col) {
    ncell.df <- rds@meta.data 

    ncell.df <- ncell.df %>% group_by_at(celltype_col) %>% 
        summarise(ncell = n()) 

    total <- nrow(rds@meta.data)
    cell.pct <- ncell.df %>% mutate(pct = 100*ncell/total) 
    cell.pct <- cell.pct %>% select(!!!syms(c(celltype_col, "pct")))
    names(cell.pct) <- c('cluster', 'cell_pct')
    print(cell.pct)
    return(cell.pct)
}


summarize_diff_genes <- function(df.list, rds, celltype_col, outdir = './') {
  
  mk.outdir(outdir)
  
  marker.df <- do.call("rbind", df.list) 
  
  save_df(marker.df, str_glue("{outdir}/deg_all.csv"))
  groups <- marker.df$group %>% unique()
  
  for (gr in groups) { 
    
    p.ndegs.obj <- plot_ndegs(marker.df, gr, rds, celltype_col, fc=0.25, p_th=0.05) 
    ndeg.df <- p.ndegs.obj[['table']]
    save_df(ndeg.df, str_glue("{outdir}/{gr}_ndegs.csv"))
    save_gg(p.ndegs.obj[['plot']], str_glue("{outdir}/deg_bar_{gr}"), nrow(ndeg.df)*0.2, 4)
    save_gg(p.ndegs.obj[['plot.cellpct']], str_glue("{outdir}/deg_bar_cellpct_{gr}"), 4, 3.5)
    
    
  }
  
  groups <- marker.df$group %>% unique()
  
  # 组间曼哈顿图
  
  for (gr in groups) {
    df.use <- marker.df %>% filter(group == gr) 
    
    p <- jjVolcano(diffData = df.use, log2FC.cutoff = 0.5, tile.col = jjAnno::useMyCol("paired",n = length(df.use$cluster %>% unique)), topGeneN = 5, celltypeSize = 2, 
                   size = 3) + ggtitle(gr)
    
    save_gg(p, str_glue("{outdir}/{gr}_manhattan"), 20, 7)
  }
  
  # 火山图
  mk.outdir(str_glue("{outdir}/volcano/"))
  
  for (gr in groups) {  
    cells <- marker.df %>% filter(group == gr) %>% pull(cluster) %>% unique()
    for (cell in cells) {
      df.use <- marker.df %>% filter(group == gr) %>% filter(cluster == cell)
      p <- plot.vol(df.use) + ggtitle( str_glue("{gr}_{cell}"))
      save_gg(p, str_glue("{outdir}/volcano/{gr}_{cell}"), 6, 4)
    }
  }
  
  # 对角线图
  mk.outdir(str_glue("{outdir}/volcano_line/"))
  for (gr in groups) {  
    cells <- marker.df %>% filter(group == gr) %>% pull(cluster) %>% unique()
    for (cell in cells) {
      df.use <- marker.df %>% filter(group == gr) %>% filter(cluster == cell)
      p <- plot.diff.vol.line(df.use, str_glue("{gr}_{cell}")) + ggtitle( str_glue("{gr}_{cell}"))
      save_gg(p, str_glue("{outdir}/volcano_line/{gr}_{cell}"), 4,4)
    }
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


parser <- ArgumentParser(description = "Marker analysis script")
parser$add_argument("--rds", type = "character", required = TRUE, help = "Path to the Seurat object in rds format")
parser$add_argument("--compare_str", type = "character", required = TRUE, help = "List of comparison groups, each group is a vector of two elements")
parser$add_argument("--celltype_col", type = "character", required = TRUE, help = "Column name in meta.data that defines cell types")
parser$add_argument("--group_col", type = "character", required = FALSE, help = "Column name used for sample grouping")
parser$add_argument("--default_data", type = "character", required = FALSE, default = "None")
parser$add_argument("--outdir", type = "character", required = FALSE, default='./diff/') 
args <- parser$parse_args()


rds_path <- args$rds
comparisons_str <- args$compare_str
sample_col <- args$group_col
cluster_col <- args$celltype_col
#fc <- args$fc
#min.pct <- args$min_pct 
#p.val <- args$pval
outdir <- args$outdir
default_data <- args$default_data

comparisons <- parse_comparisons(comparisons_str)
print(str_glue("\n\n{comparisons_str} -->> {comparisons}"))
print(str_glue("rds: {rds_path}"))
print(str_glue("sample_col: {sample_col}"))
print(str_glue("cluster_col: {cluster_col}"))
#print(str_glue("fc: {fc}"))
#print(str_glue("min.pct: {min.pct}"))
#print(str_glue("p.val: {p.val}"))
print(str_glue("outdir: {outdir}\n\n"))


rds <- readRDS(rds_path)

if (default_data != "None") {
  print("change default data")
  DefaultAssay(rds) <- default_data
}

cells <- rds@meta.data[[cluster_col]] %>% unique()
df.all <- list()
inx <- 1

for (cell in cells) {
  for (s in comparisons) {
    df <- get_markers(
      rds, 
      cell, 
      g1 = s[1], 
      g2 = s[2], 
      group.by = sample_col, 
      only.p = F, 
      celltype_col =  cluster_col
      #min.pct= min.pct
    )
    inx <- inx + 1
    df.all[[inx]] <- df
    print(str_glue("finish {cell} {s[1]}  {s[2]}"))
  }
}

summarize_diff_genes(df.all,  rds, cluster_col, outdir = str_glue("{outdir}/degs/"))
