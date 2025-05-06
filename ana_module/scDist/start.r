library(scDist)
library(Seurat)
library(ggpubr)
library(tidyverse)
library(argparse)


plotFDR <- function(scd.object, title = 'avsb') {
  results <- scd.object$results
  p.value <- p.adjust(scd.object$results$p.val,method="fdr")
  dist <- scd.object$results$Dist.
  
  df.err <- data.frame(label=rownames(results),
                   se_up=results$`95% CI (upper)`,
                   se_down=results$`95% CI (low)`)
  
  df <- data.frame(Dist=dist,log10p=-log10(p.value))
  df$label <- rownames(scd.object$results)
  df <- df %>% left_join(df.err)
  df <- df %>% mutate(dist_rank = dense_rank(Dist)) %>% 
    arrange(dist_rank) %>% mutate(label = as_factor(label))
  
  outer <- df %>% filter(log10p < 1) %>% mutate(label2 = 'X')
  
  p<-df %>%
    arrange(desc(Dist)) %>%
    ggplot(aes(x=label, y=Dist, size=log10p, fill=Dist)) +
    geom_point(size = 1, color = 'red') +
    geom_point(alpha=0.9, shape=21, color="black") +
    geom_errorbar(aes(label, dist, ymin=se_down,ymax=se_up), color = 'gray30', size = 0.5, width = 0.3)+
    geom_text(outer, mapping = aes(x=label, y=Dist, label = label2), 
              color = "red", size = 4, alpha = 0.8)+
    scale_size(range = c(0.01, max(df$log10p)*1.2), name="-log10 FDR") +
    scale_fill_viridis_c(direction = 1, option = 'A') +
    theme_bw()+
    theme(legend.position="bottom") +
    ylab("Estimated distance") +
    xlab("") +
    ylim(c(0, max(df$Dist)*1.1))+
    ggtitle(title)+
    theme(legend.position = "top") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1))
  return(list(
    p = p,
    df = df
  ))
  
}

fetch.seurat <- function(rds, group_col, target) {
    cell.use <- rownames(rds@meta.data[rds@meta.data[[group_col]] %in% target,])
    return(subset(rds, cells = cell.use))
} 

get.scDist.obj <- function(rds, group_col, celltype_col, random.effects = "orig.ident") {
        sim <- list(
     Y = GetAssayData(rds, slot = "scale.data", assay = "RNA") %>% as.data.frame(),
     meta.data = rds@meta.data %>% as.data.frame()
    )

    out <- scDist(normalized_counts = sim$Y %>% as.matrix, # 标准化的数据矩阵
                meta.data = sim$meta.data, # metadata表格
                d = 20,                   # 指定用于PCA分析的维度数量为前20
                fixed.effects = group_col, # 你感兴趣的分组条件，对应metadata表格中的列名称
                random.effects = random.effects,  # 需要去除的潜在的影响因素，比如不同的样本、年龄、批次等等，对应metadata表格中的列名称
                clusters=celltype_col # 待分析的细胞类型，对应metadata表格中的列名称
    )
    return(out)
}



parse_comparisons <- function(comparison_strings) {
  result <- list()
  for (comparison in comparison_strings) {
    parts <- unlist(strsplit(comparison, "vs"))
    if (length(parts) == 2) {
      result[[length(result) + 1]] <- parts
    }
  }
  return(result)
}


save_gg <- function(p, filename, width, height, format = c('pdf', 'png')) {
    for (d in format) {
        ggplot2::ggsave(filename = str_glue("{filename}.{d}"), plot = p, width = width, height = height, device = d)
    }
}


write_res <- function(scDist_obj, outdir, celltypes, prfx, important_pct = 0.20) {

    dir.create(recursive = T, str_glue("{outdir}/{prfx}/"))
    outdir <- str_glue("{outdir}/{prfx}/")

    #p1 <- DistPlot(scDist_obj)
    res_list <- plotFDR(scDist_obj, prfx)

    #save_gg(p1, str_glue("{outdir}/{prfx}_DistPlot"), 9,9)
    ncell <- nrow(res_list$df)
    w <- min(0.25*ncell, 8)
    if (ncell <= 15) {
        w <- 6
    }
    save_gg(res_list$p, str_glue("{outdir}/{prfx}_FDRDistPlot"), w, 4)
    res_list$df %>% write_csv(str_glue("{outdir}/{prfx}_dist.csv"))
    dfs <- list()
    inx <- 1

    for (cell in celltypes) {
        
        inx <- inx + 1
        print(scDist_obj$vals[[cell]]$beta.hat %>% head)
        print(scDist_obj$gene.names %>% head)
        print(cell)

        if (is.null(scDist_obj$vals[[cell]]$beta.hat)) {
            print(str_glue("{gr}...{cell}...no results"))
        } else {
            df.all <- data.frame(value = scDist_obj$vals[[cell]]$beta.hat, 
            label = scDist_obj$gene.names, celltype = cell)

            dfs[[inx]] <- df.all
            
            df <- df.all %>% top_n(30, abs(value))

            df$color <- ifelse(df$value>0, "Positive", "Negative")

            p <- ggbarplot(df,
                    x="label", 
                    y="value", 
                    fill = "color",
                    color = "white",
                    palette = "npg",
                    sort.val = "asc",
                    sort.by.groups = FALSE,
                    xlab = "",
                    legend.title = str_glue("Gene importance score ({cell})"))+ theme_bw() + ylab("Condition difference") + coord_flip()
            p.gene.cell <- distGenes(scDist_obj, cluster = cell)

            save_gg(p, str_glue("{outdir}/{prfx}_{cell}_gene_importance"), 8,10)
            save_gg(p.gene.cell, str_glue("{outdir}/{prfx}_{cell}_distGenes"), 4,4)
            print(str_glue("{gr}...{cell}...FINISH"))
        }
    }

    df.out1 <- do.call('rbind' , dfs) %>% write_csv(str_glue("{outdir}/{prfx}_gene_importance.csv"))
    
    # 筛选用于富集分析的基因
    df.out1 <- df.out1 %>%
        group_by(celltype) %>%
        group_split() %>%
        map_dfr(~ {
            # 计算每组的行数并取 25%
            n <- nrow(.x)
            top_n <- max(1, ceiling(n * important_pct))  # 确保至少取 1 行，向上取整

            # 筛选 value > 0 的数据，取前 25%
            positive <- .x %>%
            filter(value > 0) %>%
            arrange(desc(value)) %>%
            slice_head(n = top_n)
            
            # 筛选 value < 0 的数据，按绝对值排序，取前 25%
            negative <- .x %>%
            filter(value < 0) %>%
            arrange(desc(abs(value))) %>%
            slice_head(n = top_n)
            
            # 合并两部分数据
            bind_rows(positive, negative)
        }) %>%
        ungroup()  # 去掉分组

    df.out2 <- df.out1 %>% rename(gene = label, gene=label, cluster = celltype) %>% mutate(p_val_adj = 0.00001, group = prfx)
    #df.out2$avg_log2FC <- ifelse(df.out2$value>0, 1, -1)
    #df.out2 <- df.out2 %>% select(gene, avg_log2FC, p_val_adj, cluster, group)
    df.out2 <- df.out2 %>% write_csv(str_glue("{outdir}/{prfx}_enrich_input.csv"))

    res <- scDist_obj$results %>% mutate(col = rownames(.)) %>% write_csv(str_glue("{outdir}/{prfx}_dist_results.csv"))

}

 
parser <- ArgumentParser()
parser$add_argument("--rds",default=TRUE, help="")
parser$add_argument("--group_col",default=TRUE, help="")
parser$add_argument("--celltype_col",default=TRUE, help="")
parser$add_argument("--outdir",default=TRUE, help="")
parser$add_argument("--random_effects", help="",default = 'orig_ident')
parser$add_argument("--compare_str",default=TRUE,help = 'compare str: AvsB,AvsC,BvsC')
args <- parser$parse_args()

rds <- args$rds
outdir <- args$outdir
celltype_col <- args$celltype_col
group_col <- args$group_col
random_effects <- str_replace(args$random_effects, "_", ".")
compare <- args$compare_str

print(random_effects)

dir.create(recursive = T, outdir)
rds <- readRDS(rds)
celltypes <- rds@meta.data[[celltype_col]] %>% unique

compare_gr <- parse_comparisons(stringr::str_split(compare, ",")[[1]])

print(compare_gr)

enrich_dfs <- list()
inx <- 1

for (gr in compare_gr) {
    print(gr)
    rds.use <- fetch.seurat(rds, group_col, gr)
    print(str_glue("start: {rds.use@meta.data[[group_col]] %>% unique}"))
    
    scDist_obj <- get.scDist.obj(
        rds.use, group_col, celltype_col, random.effects = random_effects)
    prfx <-  stringr::str_c(gr, collapse = "vs")
    write_res(scDist_obj, outdir, celltypes, prfx =prfx)
    saveRDS(scDist_obj, str_glue("{outdir}/{prfx}.rds"))

    enrich_df <- read_csv(str_glue("{outdir}/{prfx}/{prfx}_enrich_input.csv"))
    enrich_dfs[[inx]] <- enrich_df
    inx <- inx + 1
}

do.call('rbind' , enrich_dfs) %>% write_csv(str_glue("{outdir}/enrich_input.csv"))