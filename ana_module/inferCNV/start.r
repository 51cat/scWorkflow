library(infercnv)
library(Seurat)
library(argparse)
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

POS_DIR <- str_glue('{dirnameScript()}/data/')


mk.outdir <- function(dir) {
    if (!dir.exists(dir)){
        dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    }
    return(dir)
}

mk_inferCNV_input <- function(
    seurat_obj, 
    group_col, 
    ref_group_names_str,
    posfile_name,
    exclude_group = NULL) {
    
    if (!is.null(exclude_group)) {
        cells <- rownames(seurat_obj@meta.data[!seurat_obj@meta.data[[group_col]] %in% exclude_group, ])
        seurat_obj <- subset(seurat_obj, cells = cells)
        print(unique(seurat_obj@meta.data[[group_col]]))
    }
    counts_matrix <- GetAssayData(seurat_obj, slot="counts")
    Idents(seurat_obj) <- seurat_obj@meta.data[[group_col]]
    group_df <- data.frame(Idents(seurat_obj))
    ref_pos <- stringr::str_glue("{POS_DIR}/{posfile_name}.txt")
    ref_group_names_use <-  stringr::str_split(ref_group_names_str, ",")[[1]]

    # input
    input_obj <- list(
            raw_counts_matrix = counts_matrix,
            annotations_file = group_df,
            gene_order_file = ref_pos,
            ref_group_names=ref_group_names_use
        )

    return(input_obj)       
}

run_inferCNV <- function(input_obj, outdir, cutoff = 0.1) {
    out.p <- mk.outdir(stringr::str_glue("{outdir}/inferCNV_output/"))

    if (file.exists(stringr::str_glue("{out.p}/preliminary.infercnv_obj"))) {
        print("find infercnv obj file!")
        infercnv_obj <- readRDS(stringr::str_glue("{out.p}/preliminary.infercnv_obj"))
        plot_cnv(infercnv_obj,out_dir = stringr::str_glue("{out.p}/infercnv.png"))

    }else{
        infercnv_obj = CreateInfercnvObject(raw_counts_matrix=input_obj[['raw_counts_matrix']],
                                        annotations_file=input_obj[['annotations_file']],
                                        delim="\t",
                                        gene_order_file=input_obj[['gene_order_file']],
                                        ref_group_names=input_obj[['ref_group_names']],
                                        max_cells_per_group=NULL,
                                        min_max_counts_per_cell=c(1, +Inf))

        infercnv_obj = infercnv::run(infercnv_obj,
                                    cutoff=cutoff, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                                    out_dir=out.p, 
                                    cluster_by_groups=TRUE, 
                                    denoise=TRUE,
                                    HMM=TRUE,
                                    write_expr_matrix=TRUE,
                                    output_format="pdf")
        plot_cnv(infercnv_obj,
                out_dir = stringr::str_glue("{out.p}/infercnv.png"))
    }

    }

parser <- ArgumentParser(description = "")

parser$add_argument("--rds", type = "character", default = "",
                    help = "Path to the RDS file")
parser$add_argument("--group_col", type = "character", default = "",
                    help = "group column name")
parser$add_argument("--ref_group_names_str", type = "character", 
                    default = "",
                    help = "")
parser$add_argument("--posfile", type = "character", 
                    default = "refdata_gex_GRCh38_2024_A_human_pos",
                    help = "")
parser$add_argument("--exclude_group", type = "character", 
                    default = "None",
                    help = "")

parser$add_argument("--outdir", type = "character", default = "",
                    help = "Output directory")


args <- parser$parse_args()

print(args$rds)
print(args$group_col)
print(args$ref_group_names_str)
print(args$posfile)
print(args$exclude_group)
print(args$outdir)

rds <- readRDS(args$rds)


if (args$exclude_group == 'None') {
    exclude_group <- NULL
}else{
    exclude_group <- stringr::str_split(args$exclude_group, ",")[[1]]
}

print('exclude_group')
print(exclude_group)

input_obj <- mk_inferCNV_input(
    seurat_obj = rds, 
    group_col = args$group_col, 
    ref_group_names_str = args$ref_group_names_str,
    posfile_name = args$posfile,
    exclude_group = exclude_group)


print( unique(input_obj[['annotations_file']][,1] ))
print('========')
print( input_obj[['ref_group_names']])

run_inferCNV(input_obj, args$outdir)
print('Finish!')