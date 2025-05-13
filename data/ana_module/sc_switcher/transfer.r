library(scCustomize)
library(stringr)
library(tidyverse)
library(argparse)


#Sys.setenv(PYTHONUSERBASE = "/usr/bin/python")
#Sys.getenv("PYTHONUSERBASE")

library(reticulate)
#print(py_config())

# Seurat to h5ad

get_dir_filename <- function(file_path) {
    file_dir <- dirname(outfile)
    file_name <- str_split(outfile, "/")[[1]] %>% rev() 
    file_name <- file_name[1]
    return(fdir = file_dir, fname = file_name)
} 


convert_seurat_to_h5ad_scCustomize <- function(in_rds, outdir, pyexec = 'None', split_col = 'None', keep = 'None', keep_str='None') {  
  rds <- readRDS(in_rds)
  reticulate::use_python(pyexec)
  print(py_config())
  #reticulate::py_module_available(module = 'anndata')
  dir.create(outdir, recursive = T)

  if (keep_str != 'None') {
    cells <- rds@meta.data %>% filter(!!rlang::parse_expr(keep_str)) %>% rownames(.)
    rds <- subset(rds, cells = cells)
  }

  if (split_col == 'None') {
    as.anndata(x = rds, file_path = outdir, file_name = "output.h5ad")
    return(1)
  }

  
  if (keep != 'None') {
    keep <- str_split(keep, ",")[[1]]
    cells <- rownames(rds@meta.data[rds@meta.data[[split_col]] %in% keep, ])
    rds <- subset(rds, cells = cells)
  }

  if (split_col != 'None') {
    gr.all <- unique(rds@meta.data[[split_col]])

    # split
    for (gr in gr.all) {
        cells <- rownames(rds@meta.data[rds@meta.data[[split_col]] == gr, ])
        rds.use <- subset(rds, cells = cells)
        as.anndata(x = rds.use, file_path = outdir, file_name = str_glue("{gr}.h5ad"))
    }
  }
}



parser <- ArgumentParser(description = 'Convert Seurat object to h5ad with optional filtering and splitting.')

parser$add_argument('--rds', required = TRUE, help = 'Input .rds file (Seurat object)')
parser$add_argument('--outdir', required = TRUE, help = 'Output directory to save .h5ad files')
parser$add_argument('--pyexec', default = '/usr/bin/python', help = 'Python executable path for reticulate')
parser$add_argument('--split_col', default = 'None', help = 'Metadata column for splitting')
parser$add_argument('--keep', default = 'None', help = 'Comma-separated values to keep in split_col')
parser$add_argument('--keep_str', default = 'None', help = 'Filtering expression for metadata, e.g. "nFeature_RNA > 200 & percent.mt < 5"')

args <- parser$parse_args()

convert_seurat_to_h5ad_scCustomize(
  in_rds = args$rds,
  outdir = args$outdir,
  pyexec = args$pyexec,
  split_col = args$split_col,
  keep = args$keep,
  keep_str = args$keep_str
)