library(scCustomize)
library(stringr)
library(tidyverse)
library(argparse)
library(reticulate)
library(Matrix)

########### function ##############
get_dir_filename <- function(file_path) {
    file_dir <- dirname(outfile)
    file_name <- str_split(outfile, "/")[[1]] %>% rev() 
    file_name <- file_name[1]
    return(fdir = file_dir, fname = file_name)
} 

extract_mtx <- function(rds, outdir) {
  dir.create(outdir, recursive = T)
  sparse <- Matrix(rds@assays$RNA$counts,sparse = T)
  features <- data.frame(ID = sparse@Dimnames[[1]], Name = sparse@Dimnames[[1]], EX = "Gene Expression")
  write.table(x = features,file = str_glue("{outdir}/features.tsv"),sep = "\t",quote = F,col.names = F,row.names = F)
  write(x = sparse@Dimnames[[2]],file = str_glue("{outdir}/barcodes.tsv"))
  writeMM(sparse,file = str_glue("{outdir}/matrix.mtx"))
}

parse_keep_str <- function(file_path) {
  line <- readLines(file_path, warn = FALSE, n = 1)
  return(line)
}

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

to_count_h5ad <- function(mtx_file, meta, outfile, is_normalize = FALSE) {
  cmd <- str_glue("python {dirnameScript()}/to_count_h5ad.py --matrix_file {mtx_file} --meta {meta} --outfile {outfile}")
  if (is_normalize %in% c(TRUE, 'True')) {
    cmd <- str_glue("python {dirnameScript()}/to_count_h5ad.py --matrix_file {mtx_file} --meta {meta} --outfile {outfile} --is_normailze True")
  }
  system(cmd)
}

##################################

convert_seurat_to_h5ad_scCustomize <- function(in_rds, outdir, pyexec = 'None', split_col = 'None', keep = 'None', keep_str='None') {  
  rds <- readRDS(in_rds)
  reticulate::use_python(pyexec)
  print(py_config())
  #reticulate::py_module_available(module = 'anndata')
  dir.create(outdir, recursive = T)

  if (keep_str != 'None') {
    keep_str <- parse_keep_str(keep_str)
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

convert_seurat_to_count <-  function(in_rds, outdir, split_col = 'None', keep = 'None', keep_str='None', to_h5ad = FALSE, is_normailze_h5ad = FALSE) {
  rds <- readRDS(in_rds)
  dir.create(outdir, recursive = T)

  if (keep_str != 'None') {
    keep_str <- parse_keep_str(keep_str)
    cells <- rds@meta.data %>% filter(!!rlang::parse_expr(keep_str)) %>% rownames(.)
    rds <- subset(rds, cells = cells)
  }

  if (split_col == 'None') {
    extract_mtx(rds, str_glue("{outdir}/MTX/"))
    system(str_glue("gzip {outdir}/MTX/*"))

    rds@meta.data$barcode <- rownames(rds@meta.data)
    rds@meta.data %>% write_tsv(str_glue("{outdir}/metadata.tsv"))
    
    # TO count h5ad
    if (to_h5ad %in% c(TRUE, 'True')) {
      to_count_h5ad(str_glue("{outdir}/MTX/"), str_glue("{outdir}/metadata.tsv"),str_glue("{outdir}/rna.h5ad"), is_normalize = is_normailze_h5ad)
    }

    return(0)
  }

  if (keep != 'None') {
    keep <- str_split(keep, ",")[[1]]
    cells <- rownames(rds@meta.data[rds@meta.data[[split_col]] %in% keep, ])
    rds <- subset(rds, cells = cells)
  }

  if (split_col != 'None') {
    gr.all <- unique(rds@meta.data[[split_col]])
    rds@meta.data$barcode <- rownames(rds@meta.data)
    rds@meta.data %>% write_tsv(str_glue("{outdir}/metadata.tsv"))

    # split
    for (gr in gr.all) {
        cells <- rownames(rds@meta.data[rds@meta.data[[split_col]] == gr, ])
        rds.use <- subset(rds, cells = cells)
        extract_mtx(rds.use, str_glue("{outdir}/MTX/{gr}/"))
        system(str_glue("gzip {outdir}/MTX/{gr}/*"))
        # TO count h5ad
      
      if (to_h5ad %in% c(TRUE, 'True')) {
        to_count_h5ad(str_glue("{outdir}/MTX/{gr}"), str_glue("{outdir}/metadata.tsv"),str_glue("{outdir}/{gr}.h5ad"), is_normalize = is_normailze_h5ad)
      }

    }
  }
}


convert_h5ad_to_seurat <- function(in_h5ad_file, outdir, split_col = 'None', keep = 'None', keep_str='None') {
  dir.create(outdir, recursive = T)
  sceasy::convertFormat(in_h5ad_file,  from="anndata", to="seurat", outFile=str_glue('{outdir}/_output.rds'))
  rds <- readRDS(str_glue('{outdir}/_output.rds'))
  rds_V5 <- Convert_Assay(seurat_object = rds, convert_to = "V5")
  #saveRDS(rds_V5, str_glue('{outdir}/output_v5.rds'))
  system(str_glue('rm -rf {outdir}/_output.rds'))
  
  # 拆分

  if (keep_str != 'None') {
    keep_str <- parse_keep_str(keep_str)
    cells <- rds@meta.data %>% filter(!!rlang::parse_expr(keep_str)) %>% rownames(.)
    rds <- subset(rds, cells = cells)
  }

  if (split_col == 'None') {
    saveRDS(rds_V5, str_glue('{outdir}/output.rds'))
    return(1)
  }

  
  if (keep != 'None') {
    keep <- str_split(keep, ",")[[1]]
    cells <- rownames(rds_V5@meta.data[rds_V5@meta.data[[split_col]] %in% keep, ])
    rds_V5 <- subset(rds_V5, cells = cells)
  }

  if (split_col != 'None') {
    gr.all <- unique(rds_V5@meta.data[[split_col]])

    # split
    for (gr in gr.all) {
        cells <- rownames(rds_V5@meta.data[rds_V5@meta.data[[split_col]] == gr, ])
        rds.use <- subset(rds_V5, cells = cells)
        saveRDS(rds.use, str_glue('{outdir}/{gr}.rds'))
    }
  }



}


parser <- ArgumentParser(description = 'Convert Seurat object to h5ad with optional filtering and splitting.')

parser$add_argument('--rds', required = TRUE, help = 'Input .rds file (Seurat object)')
parser$add_argument('--outdir', required = TRUE, help = 'Output directory to save .h5ad files')
parser$add_argument('--pyexec', default = '/usr/bin/python', help = 'Python executable path for reticulate')
parser$add_argument('--split_col', default = 'None', help = 'Metadata column for splitting')
parser$add_argument('--keep', default = 'None', help = 'Comma-separated values to keep in split_col')
parser$add_argument('--keep_str', default = 'None', help = 'A file of filtering expression for metadata, e.g. "nFeature_RNA > 200 & percent.mt < 5"')

parser$add_argument('--from', default = 'None', help = '')
parser$add_argument('--to', default = 'None', help = '')

args <- parser$parse_args()


if (args$from == 'seurat' & args$to == 'h5ad') {
  convert_seurat_to_h5ad_scCustomize(
    in_rds = args$rds,
    outdir = args$outdir,
    pyexec = args$pyexec,
    split_col = args$split_col,
    keep = args$keep,
    keep_str = args$keep_str
  )
}

if (args$from == 'seurat' & args$to == 'mtx') {
  convert_seurat_to_count(
    in_rds = args$rds,
    outdir = args$outdir,
    split_col = args$split_col,
    keep = args$keep,
    keep_str = args$keep_str
  )
}

if (args$from == 'seurat' & args$to == 'count_h5ad') {
  convert_seurat_to_count(
    in_rds = args$rds,
    outdir = args$outdir,
    split_col = args$split_col,
    keep = args$keep,
    keep_str = args$keep_str,
    to_h5ad = TRUE
  )
}

if (args$from == 'seurat' & args$to == 'normalize_h5ad') {
  convert_seurat_to_count(
    in_rds = args$rds,
    outdir = args$outdir,
    split_col = args$split_col,
    keep = args$keep,
    keep_str = args$keep_str,
    to_h5ad = TRUE,
    is_normailze_h5ad = TRUE
  )
}

if (args$from == 'h5ad' & args$to == 'seurat') {
  convert_h5ad_to_seurat(
    in_h5ad_file = args$rds, outdir = args$outdir, split_col = args$split_col, keep = args$keep, keep_str=args$keep_str
  )
}