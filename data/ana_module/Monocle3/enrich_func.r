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
plot.enrich <- function(enrich_obj, top = 30) {
  enrichplot::dotplot(enrich_obj,showCategory = top, font.size = 8,  label_format = 50, color = "p.adjust" )+
    scale_fill_viridis_c()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 
}

analysis_GO <- function(markers, orgDb = 'org.Hs.eg.db', output_prefix = "./GO_Enrichment", ont = 'BP', gr_name = 'group') {
  if (length(markers) < MIN_GENES) {
    message(str_glue("Gene number very small {nrow(markers)} !"))
    return(NULL)
  }
  
  # Step 1: Convert gene symbols to gene IDs
  geneID <- bitr(markers, fromType = "SYMBOL", 
                 toType = c("ENTREZID", "SYMBOL"), OrgDb = orgDb)
  
  
  # Step 2: Perform GO enrichment analysis
  all.GO <- enrichGO(gene = geneID$ENTREZID, OrgDb = orgDb, 
                     keyType = "ENTREZID", ont = ont, 
                     pAdjustMethod = "BH", readable = TRUE)
  
  # Check if there are any significant results
  if (is.null(all.GO) || nrow(as.data.frame(all.GO)) == 0) {
    message("No significant GO enrichment results found.")
    return(NULL)
  }
  
  # Step 3: Save enrichment results
  allgg <- as.data.frame(all.GO) %>% mutate(group = gr_name)
  output_table <- paste0(output_prefix, str_glue("{ont}_Results.csv")) 
  write.csv(allgg, file = output_table, row.names = FALSE)
  message("GO enrichment results saved as: ", output_table)
  p <- plot.enrich(all.GO)
  
  w <- 6
  h <- 10
  
  save_gg(p, paste0(output_prefix, str_glue("_{ont}_top30")), w, h )
}

analysis_KEGG <- function(markers, orgDb = 'org.Hs.eg.db',organism='hsa', output_prefix = "KEGG_Enrichment", gr_name = 'group'){
  
  if (length(markers) < MIN_GENES) {
    message(str_glue("Gene number very small {nrow(markers)} !"))
    return(NULL)
  }
  
  # Step 1: Convert gene symbols to Entrez IDs
  geneID <- bitr(markers, fromType = "SYMBOL", 
                 toType = c("ENTREZID", "SYMBOL"), OrgDb = orgDb)
  
  # Step 2: Perform KEGG enrichment analysis
  kegg_result <- enrichKEGG(gene = geneID$ENTREZID, 
                            organism = organism, 
                            pvalueCutoff = 0.05)
  
  # Check if there are any significant results
  if (is.null(kegg_result) || nrow(as.data.frame(kegg_result)) == 0) {
    message("No significant KEGG enrichment results found.")
    return(NULL)
  }
  
  kegg_result <- setReadable(kegg_result,
                             OrgDb = orgDb,
                             keyType = "ENTREZID")
  
  # Step 3: Save enrichment results
  kegg_df <- as.data.frame(kegg_result)%>% mutate(group = gr_name)
  output_table <- paste0(output_prefix, "_kegg_Results.csv")
  write.csv(kegg_df, file = output_table, row.names = FALSE)
  message("KEGG enrichment results saved as: ", output_table)
  
  p <- plot.enrich(kegg_result)
  
  w <- 6
  h <- 10
  
  save_gg(p,paste0(output_prefix, str_glue("_kegg_top30")) , w, h )
}


analysis_Reactome <- function(markers, orgDb = 'org.Hs.eg.db',organism = 'human', output_prefix = "Reactome_Enrichment",gr_name = 'group') {
  if (length(markers) < MIN_GENES) {
    message(str_glue("Gene number very small {nrow(markers)} !"))
    return(NULL)
  }
  
  # Step 1: Convert gene symbols to Entrez IDs
  geneID <- bitr(markers, fromType = "SYMBOL", 
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
  reactome_df <- as.data.frame(reactome_result) %>% mutate(group = gr_name)
  output_table <- paste0(output_prefix, "_reactome_Results.csv")
  write.csv(reactome_df, file = output_table, row.names = FALSE)
  message("Reactome enrichment results saved as: ", output_table)
  
  p <- plot.enrich(reactome_result)
  
  w <- 6
  h <- 10
  
  save_gg(p, paste0(output_prefix, str_glue("_reactome_top30")) , w, h )
  
}







