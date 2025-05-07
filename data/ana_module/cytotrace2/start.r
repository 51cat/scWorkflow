## new code 
library(CytoTRACE2)
library(tidyverse)
library(argparse)

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


process_data <- function(rds, species, celltype_use = "all") {
    
    if (celltype_use != 'all') {
        celltype_input <- str_split(celltype_use, ",") %>% unlist
        print(celltype_input)
        cells <-  rownames(rds@meta.data[rds@meta.data[[celltype_col]] %in% celltype_input,])
        rds <- subset(rds, cells = cells)
    }

    if (ncol(rds) > 50000) {
        # https://github.com/digitalcytometry/cytotrace2/issues/44
         print("close parallelize..")
         cytotrace2_result <- cytotrace2(rds, species = species, is_seurat = T, parallelize_models = FALSE, parallelize_smoothing = FALSE)
    }else{
        cytotrace2_result <- cytotrace2(rds, species = species, is_seurat = T)
    }
    return(cytotrace2_result)

}

change_reduction <- function(plot_obj, seurat_obj, reduction_name) {

    if (F) {
        print(str_glue("change reduction -> {reduction_name}"))

        emb_1 <- seurat_obj@reductions[[reduction_name]]@cell.embeddings[,1]  # First dimension
        emb_2 <- seurat_obj@reductions[[reduction_name]]@cell.embeddings[,2]  # Second dimension

        plot.all <- c("CytoTRACE2_UMAP", "Phenotype_UMAP", "CytoTRACE2_Potency_UMAP","CytoTRACE2_Relative_UMAP")
        for (plot in plot.all) {
            plot_obj[[plot]][[1]][["data"]]["UMAP_1"] <- emb_1
            plot_obj[[plot]][[1]][["data"]]["UMAP_2"] <- emb_2
        }
    }
        
        return(plot_obj)
}


out.result <- function(cytotrace2_result, seurat_obj, outdir, celltype_col, prfx, group_col = "", reduction_name = 'umap') {
    out.p <- mk.outdir(str_glue("{outdir}/{prfx}/"))

    write.csv(cytotrace2_result@meta.data %>% select(contains(celltype_col) ,contains('TRACE2')), 
            str_glue("{out.p}/{prfx}_cytotrace2_result.txt"),row.name=T)
    
    print(c(celltype_col, 'orig.ident'))

    annotation <- seurat_obj@meta.data[,c(celltype_col, 'orig.ident')]

    #if (celltype_col == 'seurat_clusters') {
    #    # 防止潜在的bug
    #    names(annotation) <- c('seurat_clusters_input', 'orig.ident')
    #}

    print(annotation %>% head(10))

    plots <- plotData(cytotrace2_result = cytotrace2_result, 
                    annotation = annotation,
                    is_seurat = TRUE,
                    expression_data = NULL # Seurat object and is_seurat = TRUE, can be left NULL (default is NULL).
                    )

    
    # 替换指定的降维结果: 如果降维度结果不是umap/UMAP, 则使用指定的降维结果展示(v1.1.0后版本支持)
    # 暂不启用, 未测试新版本
    # https://github.com/digitalcytometry/cytotrace2/issues/49
    if (!reduction_name %in% c('UMAP', 'umap')) {
        plots <- change_reduction(plots, seurat_obj, reduction_name = reduction_name)
    }

    p1 <- plots$CytoTRACE2_Boxplot_byPheno + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

    save_gg(p1, str_glue("{out.p}/{prfx}_CytoTRACE2_Boxplot_byPheno"), 6,6)
    save_gg(plots$CytoTRACE2_Potency_UMAP, str_glue("{out.p}/{prfx}_CytoTRACE2_Potency_UMAP"), 6,6)
    save_gg(plots$CytoTRACE2_Relative_UMAP, str_glue("{out.p}/{prfx}_CytoTRACE2_Relative_UMAP"), 6,6)
    save_gg(plots$CytoTRACE2_UMAP, str_glue("{out.p}/{prfx}_CytoTRACE2_UMAP"), 6,6)
    save_gg(plots$Phenotype_UMAP, str_glue("{out.p}/{prfx}_Phenotype_UMAP"), 6,6)

    if (group_col != "") {
        # group 1
        out.p <- mk.outdir(str_glue("{outdir}/{prfx}/group/")) 
        # 分组展示, 每组作图然后合并
        gr.all <- cytotrace2_result@meta.data[[group_col]] %>% unique()
        
        ann.group <- cytotrace2_result@meta.data[,c(group_col, celltype_col)] %>% 
                    mutate(group_cell = str_c(
                        cytotrace2_result@meta.data[[group_col]], "_", 
                        cytotrace2_result@meta.data[[celltype_col]])) %>% select(group_cell)
        
        plots <- plotData(cytotrace2_result = cytotrace2_result, 
                    annotation = ann.group,
                    is_seurat = TRUE,
                    expression_data = NULL # Seurat object and is_seurat = TRUE, can be left NULL (default is NULL).
                    )
        
        if (!reduction_name %in% c('UMAP', 'umap')) {
            plots <- change_reduction(plots, seurat_obj, reduction_name = reduction_name)
        }

        p1 <- plots$CytoTRACE2_Boxplot_byPheno + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

        save_gg(p1, str_glue("{out.p}/{prfx}_group_{group_col}_CytoTRACE2_Boxplot_byPheno"), 6,6)
        save_gg(plots$CytoTRACE2_Potency_UMAP, str_glue("{out.p}/{prfx}_group_{group_col}_CytoTRACE2_Potency_UMAP"), 6,6)
        save_gg(plots$CytoTRACE2_Relative_UMAP, str_glue("{out.p}/{prfx}_group_{group_col}_CytoTRACE2_Relative_UMAP"),6,6)
        save_gg(plots$CytoTRACE2_UMAP, str_glue("{out.p}/{prfx}_group_{group_col}_CytoTRACE2_UMAP"), 6,6)
        save_gg(plots$Phenotype_UMAP, str_glue("{out.p}/{prfx}_group_{group_col}_Phenotype_UMAP"), 6,6)

        # group2
        out.p <- mk.outdir(str_glue("{outdir}/{prfx}/group2/")) 

        gr.all <- cytotrace2_result@meta.data[[group_col]] %>% unique()
        ann.group <- cytotrace2_result@meta.data[,c(group_col, 'orig.ident')]

        plots <- plotData(cytotrace2_result = cytotrace2_result, 
                    annotation = ann.group,
                    is_seurat = TRUE,
                    expression_data = NULL # Seurat object and is_seurat = TRUE, can be left NULL (default is NULL).
                    )
        
        if (!reduction_name %in% c('UMAP', 'umap')) {
            plots <- change_reduction(plots, seurat_obj, reduction_name = reduction_name)
        }

        p1 <- plots$CytoTRACE2_Boxplot_byPheno + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

        save_gg(p1, str_glue("{out.p}/{prfx}_group2_{group_col}_CytoTRACE2_Boxplot_byPheno"), 6,6)
        save_gg(plots$CytoTRACE2_Potency_UMAP, str_glue("{out.p}/{prfx}_group2_{group_col}_CytoTRACE2_Potency_UMAP"), 6,6)
        save_gg(plots$CytoTRACE2_Relative_UMAP, str_glue("{out.p}/{prfx}_group2_{group_col}_CytoTRACE2_Relative_UMAP"),6,6)
        save_gg(plots$CytoTRACE2_UMAP, str_glue("{out.p}/{prfx}_group2_{group_col}_CytoTRACE2_UMAP"), 6,6)
        save_gg(plots$Phenotype_UMAP, str_glue("{out.p}/{prfx}_group2_{group_col}_Phenotype_UMAP"), 6,6)

    
    }
}


#rds_path = "/home/zhliu/sc_project/pipline/pipline_test/obj.rds"
#celltype_col = 'celltype'
#celltype_use = 'all'
#group_col = "split"
#species = "human"
#outdir = "./output_res/"


parser <- ArgumentParser(description = "")

# 添加命令行参数
parser$add_argument("--rds",
                    help = "")
parser$add_argument("--celltype_col",
                    help = "")
# 弃用
parser$add_argument("--group_col",
                    help = "", default = "")
parser$add_argument("--celltype_use",
                    default = "all",
                    help = "")
parser$add_argument("--spec",
                    default = "human")
                    
parser$add_argument("--reduction_name",
                    default = "umap")

parser$add_argument("--outdir",
                    help = "")

# 解析命令行参数
args <- parser$parse_args()

rds_path <- args$rds
celltype_col <- args$celltype_col
group_col <- args$group_col
celltype_use <- args$celltype_use
species <- args$spec
outdir <- args$outdir
reduction_name <- args$reduction_name

rds <- readRDS(rds_path)

cytotrace2_result <- process_data(rds, celltype_use, species =species)
out.result(cytotrace2_result,rds, outdir, celltype_col, celltype_use, group_col = group_col, reduction_name = reduction_name)
saveRDS(cytotrace2_result, str_glue("{outdir}/cyto_res.rds"))