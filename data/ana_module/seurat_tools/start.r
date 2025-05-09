#!/public/home/zhliu/anaconda3/envs/seurat_v5/bin/Rscript
library(Seurat)
library(tidyverse)
library(argparse)
library(Seurat) # v5.0
# https://satijalab.org/seurat/articles/announcements.html
library(tidyverse)
library(argparser)
library(harmony)

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


#### QC ####

calculate_metadata_stats <- function(seurat_obj) {
  metadata_cols <- c("nFeature_RNA", "percent.mt", "nCount_RNA")
  metadata <- seurat_obj@meta.data[, metadata_cols]
  summary_stats <- function(x) {
    return(data.frame(
      min = min(x, na.rm = TRUE),
      max = max(x, na.rm = TRUE),
      mean = mean(x, na.rm = TRUE),
      median = median(x, na.rm = TRUE),
      q1 = quantile(x, 0.25, na.rm = TRUE),
      q3 = quantile(x, 0.75, na.rm = TRUE)
    ))
  }
  stats_result <- do.call(rbind, lapply(metadata, summary_stats))
  
  return(stats_result)
}

save_violin_plot <- function(rds, feature, group, file_name, gcpal ='npg') {
  
  if (!feature %in% names(rds@meta.data)) (
    return(NULL)
  )

  vec <- rds@meta.data[[group]] %>% unique
  pal3 <- get_color(vec, length(vec), palette = gcpal)
  
  plot <- VlnPlot(rds, pt.size=0, features=feature, group.by=group) + scale_fill_manual(values=pal3) + 
    NoLegend() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  ggsave(filename=file_name, plot=plot, width=length(vec)*0.25, height=6)
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
argv <- add_argument(argv, "--outdir", type = "character", default = '',
                    help = "")

argv <- add_argument(argv, "--feature_show", type = "character", default = 'nCount_RNA,nFeature_RNA,percent.mt,percent.rb,percent.HB',
                    help = "")

argv <- add_argument(argv, "--gcpal", type = "character", default = 'npg',
                    help = "")                   

argv <- add_argument(argv, "--do", type = "character", default = 'None',
                    help = "")

argv <- parse_args(argv)


rds <- readRDS(argv$rds)

#### 重聚类 
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


#### 抽样
if (argv$do == 'downsample') {
  rds <- downsample.seurat(rds,p = argv$prop, g = argv$group_by)
  saveRDS(rds, "./out_downsample.rds")
}


#### qc show
feature_show <- stringr::str_split(argv$feature_show, ',')[[1]]
dir.create(argv$outdir, recursive = TRUE)
if (argv$do == 'qc_display') {
  for (f in feature_show) {
    save_violin_plot(rds,f, argv$group_by, file.path(argv$outdir, str_glue("{f}.png")), gcpal = argv$gcpal)
  }
  umap_plot <- FeaturePlot(rds, features=feature_show,  raster = FALSE, ncol = length(feature_show))
  ggsave(filename=file.path(argv$outdir, "umap.png"), plot=umap_plot, width=4*length(feature_show), height=4)
  

  stat.filter <- calculate_metadata_stats(rds)
  write.csv(stat.filter, file.path(argv$outdir, "qc.csv"))
}