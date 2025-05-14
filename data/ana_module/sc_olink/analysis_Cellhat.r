library(CellChat)
library(Seurat)
library(tidyverse)
library(argparse)
library(foreach)
library(parallel)
library(doParallel)


process_cellchat <- function(obj, celltype_col) {
    data <- GetAssayData(object = obj, slot = 'data')
    meta <- obj@meta.data
    cellchat <- createCellChat(object = obj,
                             meta = meta,
                             group.by = celltype_col,
                             assay = 'RNA')
  cellchat@DB <- CellChatDB.human #human
  cellchat <- subsetData(cellchat) #save time and memory
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  options(future.globals.maxSize = 2000 * 1024^2)
  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  #pathway level
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  # get data
  df.net<- subsetCommunication(cellchat, slot.name="netP") %>% 
    select(source, target, pathway_name, prob)

  return(
    list(cellchat = cellchat, df.net = df.net)
    )
}



parser <- ArgumentParser(description = "R script with command line arguments for corr olink cellchat")

parser$add_argument("--rds", type = "character", default = "",
                    help = "Path to the RDS file")

parser$add_argument("--group_col", type = "character", default = "",
                    help = "")

parser$add_argument("--sample_col", type = "character", default = "",
                    help = "")

parser$add_argument("--celltype_col", type = "character", default = "",
                    help = "")

args <- parser$parse_args()


rdsp <- args$rds
group_col <- args$group_col
sample_col <- args$sample_col
celltype_col <- args$celltype_col


mk_input <- function(seurat.rds.path, group_col, celltype_col) {
  
  obj <- readRDS(seurat.rds.path)
  
  objlist <- SplitObject(obj,split.by = group_col)
  print(objlist)

  return(objlist)
  # 下面的不需要
  
  # 筛选共有细胞
  all_celltypes <- lapply(objlist, function(seurat_obj) {
    unique(seurat_obj@meta.data[[celltype_col]]) 
  })
  
  common_celltypes <- Reduce(intersect, all_celltypes)
  
  print(all_celltypes)
  print(common_celltypes)
  
  seurat_list_filtered <- lapply(objlist, function(seurat_obj) {
    cell_ids <- rownames(seurat_obj@meta.data)
    
    # 找到符合共有细胞类型的细胞
    cells_to_keep <- cell_ids[seurat_obj@meta.data[[celltype_col]] %in% common_celltypes]
    
    subset(seurat_obj, cells = cells_to_keep)
  })
  
  # 输出细胞数统计信息
  print(lapply(seurat_list_filtered, function(x) {
    table(x@meta.data[[group_col]], x@meta.data[[celltype_col]])
  }))
  
  return(seurat_list_filtered)
}

rds <- readRDS(rdsp)

df.all <- list()
inx <- 1
all_samples <- rds@meta.data[[sample_col]] %>% unique

gr_sample_df <- rds@meta.data %>% 
  select(all_of(c(group_col, sample_col))) %>%
  as_tibble() %>% distinct()

names(gr_sample_df) <- c('group', 'sample')
print(gr_sample_df)


rds.all.list <- mk_input(
  rdsp,sample_col,celltype_col
)



cl <- makeCluster(min(length(rds.all.list), 4), outfile="")
registerDoParallel(cl)

res_list <- foreach(rds.sub=rds.all.list, .combine = list,.multicombine = TRUE,
        .packages=c('CellChat', 'Seurat','tidyverse'), 
        .export=c('sample_col','celltype_col', 'gr_sample_df'), 
        .verbose=TRUE) %dopar% {
    #print(sample)
    #cells <- rownames(rds@meta.data[rds@meta.data[[sample_col]] == sample,])
    #rds.sub <- subset(rds, cells = cells)

    res <- process_cellchat(rds.sub, celltype_col)[['cellchat']]
    sample_name <- res@meta[[sample_col]] %>% unique

    # 计算每个pathway的贡献度(对每个pathway求和)
    con_pathway_df <- rankNet(res, mode = "single", return.data=T, measure = 'weight', )$signaling.contribution %>% as_tibble()
    # scale: -1/log(psum)
    con_pathway_df <- con_pathway_df %>% select(name, contribution.scaled ) 
    names(con_pathway_df) <- c('gene', 'exp')
    con_pathway_df[['sample']] <- sample_name
    con_pathway_df <- con_pathway_df %>% left_join(gr_sample_df)
    con_pathway_df[['cell']] <- 'cell'
     


    #ccc_rds <- res[['cellchat']]
    #net <- res[['df.net']]
  
    #net <- net %>% mutate(sample = sample) %>% left_join(gr_sample_df)

    #net2 <- net %>% mutate(gene = str_c(source,'---'	,target,	'---',pathway_name)) %>% 
    #  mutate(cell = str_c(source, '---',	target)) %>%
    #  rename(exp = prob) %>%
    #  select(-source,	-target,	-pathway_name)

    #df.all[[inx]] <- net2
    #inx <- inx + 1

    #saveRDS(ccc_rds, str_glue("./cellchat_objs/{gr}.rds"))
    #net2 %>% write_tsv(str_glue('./{sample}_net.tsv'))
    #print(net2)
   
    #return(ccc_rds)
    return(con_pathway_df)
}
stopCluster(cl)

#cellchat.merge<-  mergeCellChat(res_list, add.names = names(res_list))

#saveRDS(cellchat.merge, "./cellchat_merge.rds")
df.out <- do.call('rbind', res_list)
df.out %>% write_tsv('./cellchat_input.tsv')






