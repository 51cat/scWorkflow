library(clusterProfiler)
library(ReactomePA)
library(ggplot2)
library(tidyverse)
library(argparse)
library(aplot)
library(ggprism)
library(foreach)
library(parallel)
library(doParallel)


# 小鼠：org.Mm.eg.db
# mmu


MIN_GENES = 3

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

##### plot function ######

plot.GO <- function(df, ntop = 30) {
  top10 <- df %>%
    group_by(ONTOLOGY) %>%
    arrange(pvalue) %>%
    slice_head(n = as.integer(ntop/3))
  
  df_top10 <- rbind(
    subset(top10, ONTOLOGY == "BP"),
    subset(top10, ONTOLOGY == "CC"),
    subset(top10, ONTOLOGY == "MF")
  )
  
  df_top10$ONTOLOGY <- factor(df_top10$ONTOLOGY, levels = c('BP', 'CC', 'MF'))
  df_top10$Description <- factor(df_top10$Description, levels = unique(rev(df_top10$Description)))
  
  mycol3 <- c('BP' = '#FF6666', 'CC' = '#6BA5CE', 'MF' = '#F5AA5F')
  p.bar <-  ggplot(data = df_top10, aes(x = -log10(pvalue), y = Description, fill = ONTOLOGY)) +
    geom_bar( stat = 'identity') 
  
  
  p.bar <- p.bar + scale_x_continuous(expand = c(0, 0.5)) +
    scale_fill_manual(values = mycol3) +
    geom_text(
      aes(x =0.2, label = Description),
      hjust =0, size =3
    )+
    labs(
      title = str_glue('Top {ntop} Enriched GO Terms (BP, CC, MF)'),
      x = '-log10(pvalue)',
      y = 'GO Term',
      fill = 'Ontology'
    ) +
    theme_prism() + 
    theme(
      axis.text.y = element_blank(),
      axis.line = element_blank(),
      axis.ticks.y = element_blank(),
      legend.title = element_text(size = 10),  # 图例标题字体大小
      legend.text = element_text(size = 10)
    )+ theme(legend.position = "left") + labs(y='')
  
  p1 <-ggplot(df_top10, aes(x='', Description, size = RichFactor, color = RichFactor)) + geom_point() + 
    scale_color_viridis_c() + labs(y = 'GO Term', x = '') +  theme_prism() +
    theme(
      axis.text.y = element_blank(),
      axis.line = element_blank(),
      axis.ticks.y = element_blank(),
      legend.title = element_text(size = 10),  # 图例标题字体大小
      legend.text = element_text(size = 10))
  
  p.bar %>% insert_left(p1, width=0.05)
  
  
  
  
}

plot.KEGG.Reactome <- function(df2, ntop = 30) {
  top30<- df2 %>%
    arrange(pvalue) %>%
    slice_head(n = ntop)
  
  
  top30$Description <- factor(top30$Description, levels = unique(rev(top30$Description)))
  
  p.bar <-  ggplot(data = top30, aes(x = -log10(pvalue), y = Description)) +
    geom_bar( stat = 'identity',  fill = 'gray') 
  
  p.bar <- p.bar + scale_x_continuous(expand = c(0, 0.5)) +
    geom_text(
      aes(x =0.2, label = Description),
      hjust =0, size =3
    )+
    labs(
      title = str_glue('Top {ntop} Enriched KEGG Terms'),
      x = '-log10(pvalue)',
      y = 'KEGG Term',
    ) +
    theme_prism() + 
    theme(
      axis.text.y = element_blank(),
      axis.line = element_blank(),
      axis.ticks.y = element_blank(),
      legend.title = element_text(size = 10),  # 图例标题字体大小
      legend.text = element_text(size = 10)
    )+ theme(legend.position = "left") + labs(y = '')
  p.bar
  
  
  p1 <-ggplot(top30, aes(x='', Description, size = RichFactor, color = RichFactor)) + geom_point() + 
    scale_color_viridis_c() + labs(y = 'KEGG Term', x = '') +  theme_prism()+
    theme(
      axis.text.y = element_blank(),
      axis.line = element_blank(),
      axis.ticks.y = element_blank(),
      legend.title = element_text(size = 10),  # 图例标题字体大小
      legend.text = element_text(size = 10))
  
  p.bar %>% insert_left(p1, width=0.05)
}

##########################

analysis_GO <- function(markers, orgDb = 'org.Hs.eg.db', output_prefix = "GO_Enrichment") {
  if (nrow(markers) < MIN_GENES) {
    message(str_glue("Gene number very small {nrow(markers)} !"))
    return(NULL)
  }
  
  # Step 1: Convert gene symbols to gene IDs
  geneID <- bitr(markers$gene, fromType = "SYMBOL", 
                 toType = c("ENTREZID", "SYMBOL"), OrgDb = orgDb)
  
  
  # Step 2: Perform GO enrichment analysis
  all.GO <- enrichGO(gene = geneID$ENTREZID, OrgDb = orgDb, 
                     keyType = "ENTREZID", ont = "ALL", 
                     pAdjustMethod = "BH", readable = TRUE)
  
  # Check if there are any significant results
  if (is.null(all.GO) || nrow(as.data.frame(all.GO)) == 0) {
    message("No significant GO enrichment results found.")
    return(NULL)
  }
  
  # Step 3: Save enrichment results
  allgg <- as.data.frame(all.GO)
  output_table <- paste0(output_prefix, "_Results.csv")  
  write.csv(allgg, file = output_table, row.names = FALSE)
  message("GO enrichment results saved as: ", output_table)
  
  p <- plot.GO(allgg)
  output_plot <- paste0(output_prefix, "_GO_bar")
  save_gg(p, output_plot, 6.4, 8.5)
}


analysis_KEGG <- function(markers, orgDb = 'org.Hs.eg.db',organism='hsa', output_prefix = "KEGG_Enrichment"){
  # Step 1: Convert gene symbols to Entrez IDs
  geneID <- bitr(markers$gene, fromType = "SYMBOL", 
                 toType = c("ENTREZID", "SYMBOL"), OrgDb = orgDb)
  
  # Step 2: Perform KEGG enrichment analysis
  kegg_result <- enrichKEGG(gene = geneID$ENTREZID, 
                            organism = organism, 
                            pvalueCutoff = 0.05, use_internal_data =TRUE)
  
  # Check if there are any significant results
  if (is.null(kegg_result) || nrow(as.data.frame(kegg_result)) == 0) {
    message("No significant KEGG enrichment results found.")
    return(NULL)
  }
  
  kegg_result <- setReadable(kegg_result,
                             OrgDb = orgDb,
                             keyType = "ENTREZID")
  
  # Step 3: Save enrichment results
  kegg_df <- as.data.frame(kegg_result)
  output_table <- paste0(output_prefix, "_Results.csv")
  write.csv(kegg_df, file = output_table, row.names = FALSE)
  message("KEGG enrichment results saved as: ", output_table)
  
  p <- plot.KEGG.Reactome(kegg_df)
  output_plot <- paste0(output_prefix, "_KEGG_bar")
  save_gg(p, output_plot, 6.4, 8.5)
}


analysis_Reactome <- function(markers, orgDb = 'org.Hs.eg.db',organism = 'human', output_prefix = "Reactome_Enrichment") {
  if (nrow(markers) < MIN_GENES) {
    message(str_glue("Gene number very small {nrow(markers)} !"))
    return(NULL)
  }
  
  # Step 1: Convert gene symbols to Entrez IDs
  geneID <- bitr(markers$gene, fromType = "SYMBOL", 
                 toType = c("ENTREZID", "SYMBOL"), OrgDb = orgDb)
  
  # Step 2: Perform Reactome enrichment analysis
  reactome_result <- enrichPathway(gene = geneID$ENTREZID, 
                                   organism = organism, 
                                   pvalueCutoff = 0.05, readable = TRUE)
  
  # Check if there are any significant results
  if (is.null(reactome_result) || nrow(as.data.frame(reactome_result)) == 0) {
    message("No significant Reactome enrichment results found.")
    return(NULL)
  }
  # Step 3: Save enrichment results
  reactome_df <- as.data.frame(reactome_result)
  output_table <- paste0(output_prefix, "_Results.csv")
  write.csv(reactome_df, file = output_table, row.names = FALSE)
  message("Reactome enrichment results saved as: ", output_table)
  
  p <- plot.KEGG.Reactome(reactome_df)
  output_plot <- paste0(output_prefix, "_Reactome_bar")
  save_gg(p, output_plot, 6.4, 8.5)
}

split_diff_df <- function(df, fc=0.25, p=0.05) {
  df.filter <- df %>% 
    filter(abs(avg_log2FC) >= fc) %>% 
    filter(p_val_adj <= p) 
  
  # up
  df.up <- df.filter %>% filter(avg_log2FC > 0)
  # down
  df.down <- df.filter %>% filter(avg_log2FC < 0)
  
  
  out.list <- list(
    up = df.up,
    down = df.down,
    all = df.filter
  )
  
  return(out.list)
}

parser <- ArgumentParser(description = "Marker analysis script")
parser$add_argument("--diff_table", type = "character", required = TRUE, help = "Path to the Seurat object in rds format")
parser$add_argument("--group_col", default = 'None', type = "character", required = FALSE, help = "Column name used for sample grouping")
parser$add_argument("--celltype_col", default = 'cluster', type = "character", required = FALSE, help = "Column name used for cluster grouping")
parser$add_argument("--spec", default = 'human', type = "character", required = FALSE, help = "human or mouse")

parser$add_argument("--use_celltype", default = 'all', type = "character", required = FALSE, help = "")
parser$add_argument("--use_group", default = 'all', type = "character", required = FALSE, help = "")

parser$add_argument("--fc", type = "double", required = FALSE, default=0.25)
parser$add_argument("--pval", type = "double", required = FALSE, default=0.05)
parser$add_argument("--outdir", type = "character", required = FALSE, default='./diff_enrich/') 
parser$add_argument("--geneset_use", type = "character", required = FALSE, default='KEGG,GO,Reactome') 
args <- parser$parse_args()

diff_table <- args$diff_table
cluster_col <- args$celltype_col
sample_col <- args$group_col
fc <- args$fc
pval <- args$pval
outdir <- args$outdir
spec <- args$spec

all_geneset <- c('KEGG','GO','Reactome')

diff.df.all <- read_csv(diff_table)

# 支持仅仅输入gene list
if (ncol(diff.df.all) == 1) {
  sample_col <- 'group'
  cluster_col <- 'cluster'
  one_col_input <- TRUE
  
  diff.df.all[['avg_log2FC']] <- 10
  diff.df.all[['p_val_adj']] <- 0.0001
  diff.df.all[[cluster_col]] <- 'all'
  diff.df.all[[sample_col]] <- 'all'
}else {
  one_col_input <- FALSE
}

print(str_glue("only gene list: {one_col_input}"))

#avg_log2FC p_val_adj

if (!'avg_log2FC' %in%  colnames(diff.df.all)) {
  print("avg_log2FC must in column names!!")
  q()
}

if (!'p_val_adj' %in%  colnames(diff.df.all)) {
  print("p_val_adj must in column names!!")
  q()
}



if (args$use_celltype != "all") {
  keep.celltype <- str_split(args$use_celltype, ",")[[1]]
  diff.df.all <- diff.df.all[diff.df.all[[cluster_col]] %in% keep.celltype, ]
  print(diff.df.all[[cluster_col]] %>% unique())
}

if (args$use_group != "all") {
  keep.group <- str_split(args$use_group, ",")[[1]]
  diff.df.all <- diff.df.all[diff.df.all[[use_group]] %in% keep.group, ]
  print(diff.df.all[[sample_col]] %>% unique())
}

if (sample_col == 'None') {
  groups <- 'None'
  diff.df.all[[sample_col]] <- groups
}else {
  groups <- diff.df.all[[sample_col]] %>% unique()
}

geneset.use.analysis <- str_split(args$geneset_use, ",")[[1]]


print(diff.df.all)

inx <- 1

if (spec == "human") {
  orgDb <- 'org.Hs.eg.db'
  organism <- 'hsa'
   ref <- 'human'

}else {
  orgDb <- 'org.Mm.eg.db'
  organism <- 'mmu'
  ref <- 'mouse'
}


for (gr in groups) {
  print(gr)
  diff.df.use <- diff.df.all[diff.df.all[[sample_col]] == gr, ] 
  cells.all <- diff.df.use %>% pull(cluster) %>% unique()
  
  cl <- makeCluster(min(length(cells.all), 16), outfile="")
  registerDoParallel(cl)
  foreach(cell=cells.all, 
          .packages=c('clusterProfiler', 'ReactomePA','tidyverse','aplot','ggprism'), 
          .export=c('diff.df.use', 'cluster_col','fc', 'pval','outdir','geneset.use.analysis', "ref", "organism", "orgDb"), 
          .verbose=TRUE) %dopar%  {
  #for (cell in cells.all) {
            print(cell)
            # out go
            gr_sufx <- str_replace_all(gr, " ", "_")
            
            df.use <-  diff.df.use %>% filter(cluster == cell) 
            df.use <- diff.df.use[diff.df.use[[cluster_col]] == cell,]
            diff.split <- split_diff_df(df.use, fc, pval)
            
            if (one_col_input) {
              diff.split[['up']] <- tibble()
              diff.split[['down']] <- tibble()
            }
            
            if ("GO" %in% geneset.use.analysis)  {
              out.p.go <- mk.outdir(str_glue('{outdir}/go/{cell}/{gr}/'))
              analysis_GO(diff.split[['up']], orgDb = orgDb, output_prefix = str_glue("{out.p.go}/{gr_sufx}_{cell}_up"))
              analysis_GO(diff.split[['down']], orgDb = orgDb, output_prefix = str_glue("{out.p.go}/{gr_sufx}_{cell}_down"))
              analysis_GO(diff.split[['all']], orgDb = orgDb, output_prefix = str_glue("{out.p.go}/{gr_sufx}_{cell}_all"))
              diff.split[['all']] %>% write_csv(str_glue('{outdir}/go/{cell}/{gr}/all_degs.csv'))
              diff.split[['up']] %>% write_csv(str_glue('{outdir}/go/{cell}/{gr}/up_degs.csv'))
              diff.split[['down']] %>% write_csv(str_glue('{outdir}/go/{cell}/{gr}/down_degs.csv'))
              
            }
            
            if ("KEGG" %in% geneset.use.analysis)  {
              out.p.kegg <- mk.outdir(str_glue('{outdir}/kegg/{cell}/{gr}/'))
              analysis_KEGG(diff.split[['up']], orgDb = orgDb, organism = organism, output_prefix = str_glue("{out.p.kegg}/{gr_sufx}_{cell}_up"))
              analysis_KEGG(diff.split[['down']], orgDb = orgDb, organism = organism,output_prefix = str_glue("{out.p.kegg}/{gr_sufx}_{cell}_down"))
              analysis_KEGG(diff.split[['all']], orgDb = orgDb, organism = organism, output_prefix = str_glue("{out.p.kegg}/{gr_sufx}_{cell}_all"))
              diff.split[['all']] %>% write_csv(str_glue('{outdir}/kegg/{cell}/{gr}/all_degs.csv'))
              diff.split[['up']] %>% write_csv(str_glue('{outdir}/kegg/{cell}/{gr}/up_degs.csv'))
              diff.split[['down']] %>% write_csv(str_glue('{outdir}/kegg/{cell}/{gr}/down_degs.csv'))
            }
            
            if ("Reactome" %in% geneset.use.analysis)  {
              out.p.Reactome <- mk.outdir(str_glue('{outdir}/Reactome/{cell}/{gr}/'))
              analysis_Reactome(diff.split[['up']], orgDb =orgDb, organism = ref, output_prefix = str_glue("{out.p.Reactome}/{gr_sufx}_{cell}_up"))
              analysis_Reactome(diff.split[['down']], orgDb = orgDb, organism = ref, output_prefix = str_glue("{out.p.Reactome}/{gr_sufx}_{cell}_down"))
              analysis_Reactome(diff.split[['all']], orgDb = orgDb, organism = ref, output_prefix = str_glue("{out.p.Reactome}/{gr_sufx}_{cell}_all"))
              diff.split[['all']] %>% write_csv(str_glue('{outdir}/Reactome/{cell}/{gr}/all_degs.csv'))
              diff.split[['up']] %>% write_csv(str_glue('{outdir}/Reactome/{cell}/{gr}/up_degs.csv'))
              diff.split[['down']] %>% write_csv(str_glue('{outdir}/Reactome/{cell}/{gr}/down_degs.csv'))
            }
          }
    stopCluster(cl)
}
