library(tidyverse)
library(argparse)
library(pheatmap)


parser <- ArgumentParser(description = "R script with command line arguments for corr olink cellchat")
parser$add_argument("--corr_file", type = "character", default = "",
                    help = "")
parser$add_argument("--olink_id_file", type = "character", default = "None",
                    help = "")    

parser$add_argument("--target_olink", type = "character", default = "None",
                    help = "")    

parser$add_argument("--outdir", type = "character", default = "",
                    help = "")    

parser$add_argument("--show_gene_name", type = "character", default = "False",
                    help = "") 

parser$add_argument("--target_gene", type = "character", default = "None",
                    help = "")   

parser$add_argument("--group_sort", type = "character", default = "None",
                    help = "")  

parser$add_argument("--gcpal", type = "character", default = "npg",
                    help = "")  

args <- parser$parse_args()


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


replace_matrix_values <- function(mat) {
  mat[mat <= 0.05] <- "*" 
  mat[mat > 0.05] <- "" 
  return(mat)
}


##### 读入配色库/富集模块 #####

source(str_glue("{dirnameScript()}/COLORS/load_color.r"))


dir.create(args$outdir ,recursive = T)

df <- vroom::vroom(args$corr_file)

if (args$group_sort!='None') {
  group_sort <- str_split(args$group_sort, ",")[[1]]
  df$group <- factor(df$group, levels = group_sort)
}

if (args$target_gene != 'None') {
  target.gene <-  read_csv(args$target_gene) %>% pull(gene) %>% unique
  df <- df %>% filter(gene %in% target.gene ) 

}


if (args$olink_id_file == 'None'){
    df$olink_gene <- df$UniProt
}else {
    olink_gene_df <- vroom::vroom(args$olink_id_file)
    df <- df %>% left_join(olink_gene_df)
}

if (args$target_olink!= 'None') {
  olink_use <- read_csv(args$target_olink) %>% pull(gene) %>% unique
  df <- df %>% filter(olink_gene %in% olink_use)
  print(df)
}



df %>% 
  select(-spearman_corr, -spearman_p) %>% 
    #filter(pearson_p <= 0.05) %>% 
    #filter(abs(pearson_corr) >= mean(abs(df$pearson_corr))) %>% 
    write_csv(str_glue('{args$outdir}/corr_TF_filter.csv'))

#colors <- colorRampPalette(c("blue", "white", "red"))(100) 
colors <- get_color(c(1:1000), n = 1000, palette = "RdBu")

gr.use <- df$group %>% unique
group_cpal <- get_color(gr.use, length(gr.use), palette = args$gcpal)
annotation_colors <- list(group=group_cpal)

cells <- df %>% pull(cell) %>% unique()

for (cell.use in cells) {

  df.use <- df %>% 
    filter(cell == cell.use) %>% 
    select(olink_gene,gene,group, pearson_corr, pearson_p) #%>% 
    #filter(abs(pearson_corr) >= 0.1)
  
  mtx <- df.use  %>%
    select(-pearson_p) %>% 
    mutate(gene = str_c(gene, "_", group)) %>% 
    select(-group) %>% 
    spread(key = gene, value = pearson_corr) %>% 
    column_to_rownames(var = 'olink_gene')
  
  p.mtx <- df.use  %>%
    select(-pearson_corr) %>% 
    mutate(gene = str_c(gene, "_", group)) %>% 
    select(-group) %>% 
    spread(key = gene, value = pearson_p) %>% 
    column_to_rownames(var = 'olink_gene')

  print(p.mtx)
  p.mtx <- replace_matrix_values(p.mtx)


  print(p.mtx)

  colann <- df.use  %>%
    mutate(gene = str_c(gene, "_", group)) %>% select(gene, group) %>% distinct() %>% column_to_rownames(var = 'gene')
  
  
  mtx[is.na(mtx)] <- 0
  
  if (args$show_gene_name == 'True'){
    show_colnames = T
  }else{
    show_colnames = F
  }

  pheatmap(mtx[, rownames(colann)], show_colnames = show_colnames, cluster_cols = F,
           annotation_col = colann,annotation_colors = annotation_colors,
           legend_breaks = c(-1, 0, 1),color = colors, fontsize_row  = 7, fontsize_col  = 6, border_color = NA, 
           filename = str_glue('{args$outdir}/{cell.use}_heatmap.pdf'), width = 8 , height = min(10, nrow(mtx)*0.4), legend = T, treeheight_row = 12)

  pheatmap(mtx[, rownames(colann)], show_colnames = show_colnames, cluster_cols = F, 
           annotation_col = colann,annotation_colors = annotation_colors,
           legend_breaks = c(-1, 0, 1),color = colors, fontsize_row  = 7, border_color = NA, fontsize_col  = 6,
           filename = str_glue('{args$outdir}/{cell.use}_heatmap.png'), width =8 , height = min(10, nrow(mtx)*0.4), legend = T, treeheight_row = 12)

  # p矩阵
  if (ncol(mtx) <= 50) {

    pheatmap(mtx[, rownames(colann)], show_colnames = show_colnames, cluster_cols = F, display_numbers =  p.mtx[, rownames(colann)], number_color = "black",fontsize_number = 18,
            annotation_col = colann,annotation_colors = annotation_colors,
            legend_breaks = c(-1, 0, 1),color = colors, fontsize_row  = 7, fontsize_col  = 6, border_color = NA, 
            filename = str_glue('{args$outdir}/{cell.use}_heatmap_p.pdf'), width = 8 , height = min(10, nrow(mtx)*0.4), legend = T, treeheight_row = 12)

    pheatmap(mtx[, rownames(colann)], show_colnames = show_colnames, cluster_cols = F, display_numbers =  p.mtx[, rownames(colann)], number_color = "black",fontsize_number = 18,
            annotation_col = colann,annotation_colors = annotation_colors,
            legend_breaks = c(-1, 0, 1),color = colors, fontsize_row  = 7, border_color = NA, fontsize_col  = 6,
            filename = str_glue('{args$outdir}/{cell.use}_heatmap_p.png'), width =8 , height = min(10, nrow(mtx)*0.4), legend = T, treeheight_row = 12)
  
  }

}