library(irGSEA)
library(Seurat)
library(tidyverse)
library(patchwork)
library(argparse)
library(foreach)
library(doParallel)
library(scCustomize)
library(clusterProfiler)
library(qs)
print(sessionInfo())

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

DB_DIR <- str_glue("{dirnameScript()}/database/")

gmt2rds <- function(gmt_file) {
  df <- clusterProfiler::read.gmt(gmt_file) %>% filter(gene != "")
  df$term <- str_replace_all(df$term, "/", "") %>% 
    str_replace_all(" ", "_")
  gene_sets <- split(df$gene, df$term)
  print(gene_sets)
  return(gene_sets)
}


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

load_geneset <- function(geneset_name) {
  
  if (str_detect(geneset_name ,".gmt")) {
    return(gmt2rds(geneset_name))
  }

  if (str_detect(geneset_name ,".rds")) {
    return(readRDS(geneset_name))
  }
  
  db_path <- str_glue("{DB_DIR}/{geneset_name}.rds")
  return(readRDS(db_path))
}

parse_compare_cols <- function(input_str, rds, group_col) {
  if (input_str != 'all') {
    pairs <- strsplit(input_str, ",")[[1]]
      result <- lapply(pairs, function(pair) {
        elements <- unlist(strsplit(pair, "vs"))
        return(elements)
      })
      return(result)
  }else {
    groups <- rds@meta.data[[group_col]] %>% unique
    return(combn(groups, 2, simplify = FALSE))

  }
}

process_data <- function(
  rds, 
  ncore = 31,
  geneset_name = "msigdb/H", 
  methods = c("AUCell","UCell","singscore","ssgsea", "JASMINE", "viper", "GSVA"), 
  minGSSize = 1, 
  maxGSSize = 500
  ) {
  

  rds <- irGSEA.score(
    object = rds, 
    assay = "RNA",
    slot = "counts", 
    seeds = 123, 
    ncores = ncore,
    min.cells = 3, 
    min.feature = 0,
    custom = T, 
    geneset = load_geneset(geneset_name), 
    msigdb = F, 
    species = "Homo sapiens",
    subcategory = NULL,
    geneid = "symbol",
    method = methods,
    aucell.MaxRank = NULL, 
    ucell.MaxRank = NULL, 
    kcdf = 'Gaussian',
    minGSSize = minGSSize,
    maxGSSize = maxGSSize)
  
  methods_use <- intersect(Seurat::Assays(rds), methods)
  
  
  rds@assays$RNA@data <- as.matrix(0)
  rds@assays$RNA@scale.data <- as.matrix(0)
  
  return(list(
    obj = rds,
    methods_use = methods_use
  ))}



parser <- ArgumentParser(description = "R script with command line arguments for irGSEA")

parser$add_argument("--rds", type = "character", default = "",
                    help = "Path to the RDS file")

parser$add_argument("--geneset", type = "character", default = "msigdb/H",
                    help = "Geneset value")

parser$add_argument("--celltype_col", type = "character", default = "",
                    help = "Cell type column name")

parser$add_argument("--group_col", type = "character", default = "",
                    help = "Group column name")

parser$add_argument("--target_cell", type = "character", default = "all",
                    help = "")

parser$add_argument("--outdir", type = "character", default = "",
                    help = "Output directory")

parser$add_argument("--methods", type = "character", 
                    default = "AUCell,UCell,singscore,ssgsea,GSVA,viper", #"AUCell,UCell,singscore,ssgsea,JASMINE,viper,GSVA", 
                    help = "Methods to use")

# 差异分析相关参数
parser$add_argument("--compare_str", type = "character", 
                    default = "compare_all",
                    help = "avsb,avsc,bvsc...")

parser$add_argument("--minGSSize", type = "character", 
                    default = '1', 
                    help = "Minimum number of genes in one gene set")
parser$add_argument("--maxGSSize", type = "character", 
                    default = "500",
                    help = "Maximal number of genes in one gene set.")


# 运行相关
parser$add_argument("--compare_method", type = "character", 
                    default = "within,between",
                    help = "")


args <- parser$parse_args()

rds_path <- args$rds
geneset <- args$geneset
celltype_col <- args$celltype_col
group_col <- args$group_col

target_cell <- args$target_cell
target_cell <- gsub("::", " ", target_cell)

outdir <- args$outdir
methods <- str_split(args$methods, ",")[[1]]
minGSSize <- as.integer(args$minGSSize)
maxGSSize <- as.integer(args$maxGSSize)


# Main code flow
rds <- readRDS(rds_path)
rds <- Convert_Assay(seurat_object = rds, convert_to = "V3")
rds <- UpdateSeuratObject(rds)
print(rds)

compares <- parse_compare_cols(args$compare_str, rds, group_col)
print(compares)

use_group <- compares %>% unlist %>% unique
print(use_group)

cells <- rownames(rds@meta.data[rds@meta.data[[group_col]] %in% use_group,])
rds <- subset(rds, cells =cells)
print(rds)

compare_method <-  str_split(args$compare_method, ",")[[1]]

#####  输入数据  ######
rds.within <- rds

# 指定某细胞类型用于分析
if (target_cell != "all") {
  use_celltype <- str_split(target_cell, ",")[[1]]
  cells <- rownames(rds@meta.data[rds@meta.data[[celltype_col]] %in% use_celltype,])
  rds.between <- subset(rds, cells =cells)
}else {
  rds.between <- rds
}

print(rds.within)
print(rds.between)

rds.p <- mk.outdir(str_glue("{outdir}/score_obj/"))


# 1. 组内分析
if ('within' %in% compare_method) {
  gr_all <- rds.within@meta.data[[group_col]] %>% unique() %>% as.character()
  print(gr_all)

  cl <- makeCluster(min(length(gr_all), 4), outfile="")
  registerDoParallel(cl)
  
  obj_list <- foreach(group_analysis_use=gr_all, .combine = list,.multicombine = TRUE,
        .packages=c('irGSEA', 'Seurat','tidyverse'), 
        .export=c('DB_DIR', 'rds.p', 'rds.within', 'geneset', 'methods', 'minGSSize', 'maxGSSize', 'group_col', 'compares'), 
        .verbose=TRUE) %dopar% {
      
      print(str_glue("start:{group_analysis_use}"))

      cells <- rownames(rds.within@meta.data[rds.within@meta.data[[group_col]] == group_analysis_use,])
      rds.use <- subset(rds.within, cells =cells)

    irgesa_analysis_result <- process_data(
        rds.use, 
        ncore = 32, 
        geneset_name = geneset, 
        methods = methods, 
        minGSSize=minGSSize, 
        maxGSSize=maxGSSize
      )
      
      irgesa_analysis_result[['use_cell']] <- 'all'
      irgesa_analysis_result[['compare_col']] <- celltype_col
      irgesa_analysis_result[['compare_list']] <- NULL
      irgesa_analysis_result[['compare_method']] <- 'within'
      irgesa_analysis_result[['group_col']] <- group_col
      irgesa_analysis_result[['prfx']] <- group_analysis_use

      print(str_glue("finish:{group_analysis_use}"))

      return(irgesa_analysis_result)
  }
stopCluster(cl)

print('save results... within')

for (obj in obj_list) {
  group_analysis_use <- obj[['prfx']]
  print(group_analysis_use)
  qsave(obj, str_glue("{rds.p}/group_{group_analysis_use}_score.qs"))
}

print('save results... within Finish!')

}


# 2. 组间分析
# 按照细胞类型分析

if ('between' %in% compare_method) {
  cell_all <- rds.between@meta.data[[celltype_col]] %>% unique() %>% as.character()
  print(cell_all)

  cl <- makeCluster(min(length(cell_all), 4), outfile="")
  registerDoParallel(cl)

  obj_list <- foreach(cell_analysis_use=cell_all, .combine = list, .multicombine = TRUE,
        .packages=c('irGSEA', 'Seurat','tidyverse'), 
        .export=c('DB_DIR', 'rds.p', 'rds.between', 'geneset', 'methods', 'minGSSize', 'maxGSSize', 'group_col', 'compares'), 
        .verbose=TRUE) %dopar% {

      print(str_glue("start:{cell_analysis_use}"))
      cells <- rownames(rds.between@meta.data[rds.between@meta.data[[celltype_col]] == cell_analysis_use,])
      prfx <- cell_analysis_use
      rds.use <- subset(rds.between, cells =cells)

      irgesa_analysis_result <- process_data(
          rds.use, 
          ncore = 32, 
          geneset_name = geneset, 
          methods = methods, 
          minGSSize=minGSSize, 
          maxGSSize=maxGSSize
      )
      irgesa_analysis_result[['use_cell']] <- cell_analysis_use
      irgesa_analysis_result[['compare_col']] <- group_col
      irgesa_analysis_result[['compare_list']] <- compares
      irgesa_analysis_result[['compare_method']] <- 'between'
      irgesa_analysis_result[['group_col']] <- group_col
      irgesa_analysis_result[['prfx']] <- cell_analysis_use
      print(str_glue("finish:{cell_analysis_use}"))
      return(irgesa_analysis_result)
  }
  stopCluster(cl)

  print('save results... between')
  for (obj in obj_list) {
    cell_analysis_use <- obj[['prfx']]
    print(cell_analysis_use)
    qsave(obj, str_glue("{rds.p}/{cell_analysis_use}_score.qs"))
  }

  print('save results... between Finish!')

  print('all Finish!')
}
