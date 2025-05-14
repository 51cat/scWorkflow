library(scMetabolism)
library(tidyverse)
library(rsvd)
library(Seurat)
library(pheatmap)
library(ComplexHeatmap)
library(ggsci)
library(argparse)
library(patchwork)



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

source(str_glue("{dirnameScript()}/utils/get_table.r"))
# method supports VISION, AUCell, ssgsea, and gsva, which VISION is the default method.
# metabolism.type supports KEGG and REACTOME, where KEGG contains 85 metabolism pathways and REACTOME contains 82 metabolism pathways.
print(sessionInfo())


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

parse_comparisons <- function(input_str) {
  comparison_strs <- strsplit(input_str, ":")[[1]]
  parsed_list <- lapply(comparison_strs, function(x) {
    elements <- strsplit(trimws(x), "vs")[[1]]
    return(elements)
  })
  return(parsed_list)
}

process_data <- function(rds, pathway_db, ident_col, reduction_name = "umap") {
    if (reduction_name != "umap") {
        rds@reductions$umap@cell.embeddings <- rds@reductions[[reduction_name]]@cell.embeddings
        colnames(rds@reductions$umap@cell.embeddings) <- c("UMAP_1", "UMAP_2")
    }
    Idents(rds) <- rds@meta.data[[ident_col]]
    rds.metabolism <- sc.metabolism.Seurat(obj = rds, imputation = F, ncores = 16, metabolism.type = pathway_db, method = "AUCell")
    return(rds.metabolism)
}


find_diff_pathway_all <- function(rds.scmeb, celltype_col) {
    scmeb.mtx <- rds.scmeb@assays$METABOLISM$score
    rds.scmeb@assays$METABOLISM <- NULL
    rds.scmeb[['SC_METABOLISM']] <- Seurat::CreateAssayObject(scmeb.mtx)
    DefaultAssay(rds.scmeb) <- 'SC_METABOLISM'
    Idents(rds.scmeb) <- rds.scmeb@meta.data[[celltype_col]]
    df <- FindAllMarkers(rds.scmeb,only.pos =  F, 
                    verbose = T, 
                    min.pct= 0.1,
                    logfc.threshold = 0.25)%>% 
        mutate(pathway_name = rownames(.)) %>%
        mutate(group = str_glue("all"))
    return(df)
}


find_diff_pathway <- function(rds.scmeb, group_col, compare_str) {
    scmeb.mtx <- rds.scmeb@assays$METABOLISM$score
    rds.scmeb@assays$METABOLISM <- NULL
    rds.scmeb[['SC_METABOLISM']] <- Seurat::CreateAssayObject(scmeb.mtx)
    DefaultAssay(rds.scmeb) <- 'SC_METABOLISM'

    Idents(rds.scmeb) <- rds.scmeb@meta.data[[group_col]]
    comapares <- parse_comparisons(compare_str)
    df.all.list <- list()
    inx <- 1
    for (s in comapares) {
        df <- FindMarkers(rds.scmeb, 
                    ident.1 = s[1], 
                    ident.2 = s[2], 
                    group.by = group_col,
                    only.pos =  F, 
                    test.use = 'bimod',
                    verbose = T, 
                    min.pct= 0.1,
                    logfc.threshold = 0.25)%>% 
        mutate(pathway_name = rownames(.)) %>%
        mutate(group = str_glue("{s[1]} vs {s[2]}"))
        df.all.list[[inx]] <- df
        inx <- inx + 1 
    }
    df.out <- do.call("rbind", df.all.list)
    return(df.out)
}


plot_all_pathway <- function(rds, pathway, outdir, pathway_db, prfx, ident_col) {
    out.p <- mk.outdir(str_glue("{outdir}/{prfx}/{pathway_db}/"))
    p.dot <- DotPlot.metabolism(obj =  rds, pathway = pathway, phenotype = ident_col, norm = "y")
    save_gg(p.dot, str_glue("{out.p}/pathway_all_dot"), 10, 15)
}


parser <- ArgumentParser(description = "")

parser$add_argument("--rds", type = "character", default = "",
                    help = "Path to the RDS file")
parser$add_argument("--celltype_col", type = "character", default = "",
                    help = "Cell type column name")
parser$add_argument("--group_col", type = "character", default = "",
                    help = "Group column name")
parser$add_argument("--outdir", type = "character", default = "",
                    help = "Output directory")
parser$add_argument("--reduction_name", type = "character", 
                    default = "umap",
                    help = "reduction_name")

parser$add_argument("--compare_str", help = "", default = 'None')

args <- parser$parse_args()

rds_path <- args$rds
celltype_col <- args$celltype_col
group_col <- args$group_col
outdir <- args$outdir
reduction_name <- args$reduction_name
compare_str <- args$compare_str


method.use <- c("VISION")
db.use <- c("KEGG", "REACTOME")

rds <- readRDS(rds_path)

group_all <- rds@meta.data[[group_col]] %>% unique()
cell_all <- rds@meta.data[[celltype_col]] %>% unique()
print(group_all)
print(cell_all)

compare_df_list <- list()
inx <- 1

score_raw_dir <- mk.outdir(str_glue("{outdir}/score_mtx/"))
score_df_list <- list()
inx_score_df <- 1

if (compare_str != 'None') {
    for (cell_analysis_use in cell_all)  {
        
        print(cell_analysis_use)
        cells <- rownames(rds@meta.data[rds@meta.data[[celltype_col]] == cell_analysis_use,])
        rds.use <- subset(rds, cells =cells)

        for (db in db.use) {
            rds.scmeb <- process_data(rds.use, db, group_col, reduction_name = reduction_name)
            saveRDS(rds.scmeb, str_glue("{score_raw_dir}/{cell_analysis_use}_{db}.rds"))
            score_df <- get_table_from_scMe_obj(rds.scmeb, phenotype = group_col)
            score_df <- score_df %>% mutate(group2 = cell_analysis_use, db = db)
            score_df_list[[inx_score_df]] <- score_df
            inx_score_df <- inx_score_df + 1
            
            pathway.all <- rownames(rds.scmeb@assays[["METABOLISM"]][["score"]])

            diff_df <- find_diff_pathway(rds.scmeb, group_col, compare_str)
            diff_df["db"] <- db
            diff_df["cell_group"] <- cell_analysis_use
            compare_df_list[[inx]] <- diff_df
            inx <- inx + 1

            print("process success!")
            plot_all_pathway(rds.scmeb, pathway.all, outdir, db, cell_analysis_use, group_col)
            print(str_glue("plot all pathway success! - {length(pathway.all)}"))
        }
        

    }
    diff.df.all <- do.call("rbind", compare_df_list) 
    diff.df.all %>% write_csv(str_glue("{outdir}/pathway_diff.csv"))
}

for (db in db.use) {
    inx <- inx+1

    rds.scmeb <- process_data(rds, db, celltype_col, reduction_name = reduction_name)
    saveRDS(rds.scmeb, str_glue("{score_raw_dir}/all_cell_{db}.rds"))
    score_df <- get_table_from_scMe_obj(rds.scmeb, phenotype = celltype_col)
    score_df <- score_df %>% mutate(group2 = "all_cell", db = db)
    score_df_list[[inx_score_df]] <- score_df
    inx_score_df <- inx_score_df + 1


    pathway.all <- rownames(rds.scmeb@assays[["METABOLISM"]][["score"]])
    print("process success!")
    plot_all_pathway(rds.scmeb, pathway.all, outdir, db, "all_cell", celltype_col)
    print(str_glue("plot all pathway success! - {length(pathway.all)}"))
    
    diff_df <- find_diff_pathway_all(rds.scmeb, celltype_col)
    diff_df["db"] <- db
    diff_df["cell_group"] <- 'all_cell'

    diff_df %>% write_csv(str_glue("{outdir}/pathway_diff_use_all_cell_{db}.csv"))
}

score_all_df <-  do.call("rbind", score_df_list) 
score_all_df %>% write_csv(str_glue("{outdir}/pathway_score_{method.use}_results.csv"))