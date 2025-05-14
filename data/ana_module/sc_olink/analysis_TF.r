library(tidyverse)
library(Seurat)

library(Seurat)
library(dplyr)
library(argparse)
calculate_median_expression <- function(seurat_obj, group_column, cell_type_column,  name = "tfsulm") {
  # 获取表达矩阵
  expression_matrix <- as.matrix(seurat_obj[[name]]@data)
  
  # 获取元数据并检查列是否存在
  metadata <- seurat_obj@meta.data
  if (!(group_column %in% colnames(metadata))) {
    stop("The specified group column does not exist in metadata.")
  }
  if (!(cell_type_column %in% colnames(metadata))) {
    stop("The specified cell type column does not exist in metadata.")
  }
  
  # 结果存储列表
  median_expression_list <- list()
  
  # 按 group 和 cell_type 分组计算
  for (group in unique(metadata[[group_column]])) {
    for (cell_type in unique(metadata[[cell_type_column]])) {
      # 获取该分组和细胞类型的细胞
      cell_ids <- rownames(metadata[metadata[[group_column]] == group & metadata[[cell_type_column]] == cell_type, ])
      
      if (length(cell_ids) > 0) {
        # 提取该分组和细胞类型的表达数据
        group_expr <- expression_matrix[, cell_ids, drop = FALSE]
        
        # 计算中位数（去除0值）
        median_values <- apply(group_expr, 1, function(x) {
          x_nonzero <- x[x != 0]
          if (length(x_nonzero) > 0) {
            median(x_nonzero)
          } else {
            NA  # 如果所有值都是 0，则返回 NA
          }
        })
        
        # 存储结果，列名格式：group_celltype
        #col_name <- paste(group, cell_type, sep = "_")
        col_name <- str_c(group, "___",cell_type)
        median_expression_list[[col_name]] <- median_values
      }
    }
  }
  
  # 转换为数据框
  median_expression_df <- as.data.frame(median_expression_list)
  rownames(median_expression_df) <- rownames(expression_matrix)
  median_expression_df$gene <- rownames(median_expression_df)
  
  return(median_expression_df)
}


parser <- ArgumentParser(description = "R script with command line arguments for corr olink cellchat")
parser$add_argument("--rds", type = "character", default = "",
                    help = "Path to the RDS file")
parser$add_argument("--group_col", type = "character", default = "",
                    help = "")
parser$add_argument("--sample_col", type = "character", default = "",
                    help = "")
parser$add_argument("--celltype_col", type = "character", default = "",
                    help = "")
args <- parser$parse_args()


rdsp <- args$rds
group_col <- args$group_col
sample_col <- args$sample_col
celltype_col <- args$celltype_col

rds <- readRDS(rdsp)

gr_sample_df <- rds@meta.data %>% 
 select(all_of(c(group_col, sample_col))) %>%
 as_tibble() %>% distinct()

names(gr_sample_df) <- c('group', 'sample')

df <- calculate_median_expression(rds,  group_column = sample_col, cell_type_column = celltype_col) %>% as_tibble()

df2 <- df %>% gather(key = "sample", value= 'exp',-gene) %>% 
  separate(sample, into = c("sample", "cell"), sep = "___") %>% 
  left_join(gr_sample_df)

df2 %>% write_tsv('./TF_input.tsv')


