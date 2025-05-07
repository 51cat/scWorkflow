#!/public/home/zhliu/anaconda3/envs/seurat_v5/bin/Rscript
library(Seurat)
library(tidyverse)
library(argparse)
library(Seurat) # v5.0
# https://satijalab.org/seurat/articles/announcements.html
library(tidyverse)
library(argparser)
library(harmony)

downsample.seurat <- function(obj, n = NULL, p = NULL, g = 'seurat_clusters') {
  obj@meta.data$barcode <- rownames(obj@meta.data)
  
  if (!is.null(n)) {
    cells.use <- obj@meta.data %>%
      group_by(across(all_of(g))) %>%
      slice_sample(n = n) %>% pull(barcode)
    message(stringr::str_glue("downsample n - {n} - {length(cells.use)}"))
  }
  
  if (!is.null(p)) {
    cells.use <- obj@meta.data %>%
      group_by(across(all_of(g))) %>%
      slice_sample(prop = p) %>% pull(barcode)
    message(stringr::str_glue("downsample prop - {p} - {length(cells.use)}"))
  }
  
  obj.use <- subset(obj, cells = cells.use)
  return(obj.use)
}

recluster.seurat <- function(rds, celltype_col, is_recluster, target_cells, resolution,  batch_col = 'None') {

  cell.use <- rownames(rds@meta.data[rds@meta.data[[celltype_col]] %in% target_cells,])
  rds <- subset(rds.all, cells = cell.use)

  if (is_recluster == FALSE) {
    print("Fisish! only subset data")
    return(rds)

  }


  # normalize & scale
  rds <- NormalizeData(
    rds, 
    normalization.method = "LogNormalize",
    scale.factor = 10000
    )

  rds <- FindVariableFeatures(
    rds, 
    selection.method = "vst", 
    nfeatures = 2000,
    verbose = FALSE
    )

  use.genes <- VariableFeatures(rds)

  rds <- ScaleData(
    rds, 
    vars.to.regress = c("nCount_RNA", "percent.mt"), 
    features = use.genes, 
    model.use = "linear"
    )

  rds <- RunPCA(
    object = rds, 
    features = use.genes, 
    do.print = FALSE
    )

  if (batch_col == "None") {
    reduction <- "pca"
  }else{
    print("remove batch")
    print(batch_col)
    rds <- RunHarmony(rds, argv$batch_col)
    reduction <- "harmony"
  }


  rds <- FindNeighbors(
    rds, 
    dims = 1:20, 
    force.recalc = TRUE, 
    reduction = reduction)

  rds <- FindClusters(
    rds, 
    resolution = resolution
    )


  # tsne and umap
  rds <- RunTSNE(
    rds, 
    dims = 1:20, 
    do.fast = TRUE, 
    check_duplicates = FALSE
    )

  rds <- RunUMAP(rds, dims=1:20, reduction = reduction)
  return(rds)

}




argv <- arg_parser('')
argv <- add_argument(argv,"--rds", help="")
argv <- add_argument(argv,"--target_cell", help="")
argv <- add_argument(argv,"--is_recluster", help="", default = 'True')
argv <- add_argument(argv,"--celltype_col", help="")
argv <- add_argument(argv,"--resolution", help="", default = 0.5)
argv <- add_argument(argv,"--batch_col", help="", default = "None")

argv <- add_argument(argv, "--prop", type = "double", default = NULL,
                    help = "Proportion of cells to sample from each Seurat cluster.")
argv <- add_argument(argv, "--group_by", type = "character", default = 'seurat_clusters',
                    help = "")

argv <- add_argument(argv, "--do", type = "character", default = 'None',
                    help = "")

argv <- parse_args(argv)


rds <- readRDS(argv$rds)


if (argv$do == 'recluster') {
  
  if (argv$is_recluster == 'True') {
    is_recluster <- TRUE
  }
  target_cells <- str_split(argv$target_cell, ",")[[1]]
  
  rds <- recluster.seurat(
    rds, argv$celltype_col, is_recluster, target_cells, resolution = as.numeric(argv$resolution), batch_col = argv$batch_col
  )

  saveRDS(rds, "./out_recluster.rds")
}

if (argv$do == 'downsample') {
  
  rds <- downsample.seurat(rds,p = argv$prop, g = argv$group_by)

  saveRDS(rds, "./out_downsample.rds")
}