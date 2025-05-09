#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparser)
  library(Seurat)
  library(DoubletFinder)
  library(tidyverse)
})

## ========== 定义参数 ==========
p <- arg_parser("Run DoubletFinder on a Seurat object")

p <- add_argument(p, "--rds", help = "输入 RDS 文件路径")
p <- add_argument(p, "--outfile", help = "输出 RDS 文件路径")
p <- add_argument(p, "--pcs", help = "主成分数，例如 '1:10'", default = "1:10")
p <- add_argument(p, "--pN", help = "DoubletFinder 参数 pN", default = 0.25)
p <- add_argument(p, "--pK", help = "DoubletFinder 参数 pK", default = 0.09)
p <- add_argument(p, "--method", help = "双细胞检测方法，默认 DoubletFinder", default = "DoubletFinder")

argv <- parse_args(p)

## ========== 参数解析 ==========
pcs_range <- eval(parse(text = argv$pcs))  

## ========== 函数定义 ==========

save_gg <- function(p, filename, width, height, format = c('pdf', 'png')) {
    for (d in format) {
        ggsave(filename = str_glue("{filename}.{d}"), plot = p, width = width, height = height, device = d)
    }
}

run_DoubletFinder <- function(obj, pcs = 1:10, pN = 0.25, pK = 0.09) {
  seu_kidney <- obj

  sweep.res.list_kidney <- paramSweep(seu_kidney, PCs = pcs, sct = FALSE)
  sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
  bcmvn_kidney <- find.pK(sweep.stats_kidney)

  annotations <- seu_kidney@meta.data[['seurat_clusters']]
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(0.075 * nrow(seu_kidney@meta.data))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

  seu_kidney <- doubletFinder(seu_kidney, PCs = pcs, pN = pN, pK = pK,
                              nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

  colnames(seu_kidney@meta.data) <- gsub("^DF\\.classifications_.*", "doublets_stat", colnames(seu_kidney@meta.data))

  out.df <- seu_kidney@meta.data %>%
    select(doublets_stat) %>%
    mutate(barcode_tmp = rownames(.)) %>%
    as_tibble()

  return(out.df)
}

find_doublets <- function(obj, method = "DoubletFinder", pcs = 1:10, pN = 0.25, pK = 0.09) {
  if (method == "DoubletFinder") {
    total_cell <- nrow(obj@meta.data)
    if (total_cell <= 100000) {
      doub_df <- run_DoubletFinder(obj, pcs, pN, pK)
    } else {
      obj_list <- SplitObject(obj, split.by = "sample")
      doub_df <- do.call('rbind', lapply(obj_list, run_DoubletFinder, pcs, pN, pK))
    }

    obj@meta.data$barcode_tmp <- rownames(obj@meta.data)
    obj@meta.data <- obj@meta.data %>%
      left_join(doub_df, by = "barcode_tmp")
    rownames(obj@meta.data) <- obj@meta.data$barcode_tmp
    obj@meta.data$barcode_tmp <- NULL
  }
  return(obj)
}

## ========== 主流程 ==========
message("读取对象：", argv$rds)
seu_obj <- readRDS(argv$rds)

message("开始运行方法：", argv$method)
message("参数 PCs: ", argv$pcs, " | pN: ", argv$pN, " | pK: ", argv$pK)

seu_obj <- find_doublets(seu_obj,
                         method = argv$method,
                         pcs = pcs_range,
                         pN = as.numeric(argv$pN),
                         pK = as.numeric(argv$pK))

p <- DimPlot(seu_obj, group.by = "doublets_stat", raster = FALSE)

dir.create(dirname(argv$outfile), recursive = T)
ggsave(filename=str_glue('{dirname(argv$outfile)}/doub_stat.pdf'), plot=p , width=4, height=4, device = 'pdf')
ggsave(filename=str_glue('{dirname(argv$outfile)}/doub_stat.png'), plot=p , width=4, height=4, device = 'png')
saveRDS(seu_obj, argv$outfile)


