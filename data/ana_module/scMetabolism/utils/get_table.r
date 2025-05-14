library(tidyverse)

# 从scmetabolism对象提取分析结果, 包括
# 1. 通路在每个group的打分中位数
# 2. 通路在每个group的打分中位数的标准化结果(绘图数据)

get_table_from_scMe_obj <- function(obj, phenotype, norm = "y"){
  input.norm <- norm
  #input.pathway <- pathway
  input.parameter<-phenotype
  
  metadata<-obj@meta.data
  metabolism.matrix <- obj@assays$METABOLISM$score
  input.pathway <- rownames(metabolism.matrix)
  
  metadata[,input.parameter]<-as.character(metadata[,input.parameter])
  metabolism.matrix_sub<-t(metabolism.matrix[input.pathway,])
  
  #arrange large table
  gg_table<-c()
  for (i in 1:length(input.pathway)){
    gg_table<-rbind(gg_table, cbind(metadata[,input.parameter], input.pathway[i], metabolism.matrix_sub[,i]))
  }
  gg_table<-data.frame(gg_table)
  
  #get median value
  gg_table_median<-c()
  input.group.x<-unique(as.character(gg_table[,1]))
  input.group.y<-unique(as.character(gg_table[,2]))
  # get mean vale 
  # TO DO: 增加均值的输出
  
  
  for (x in 1:length(input.group.x)){
    for (y in 1:length(input.group.y)){
      gg_table_sub<-subset(gg_table, gg_table[,1] == input.group.x[x] & gg_table[,2] == input.group.y[y])
      gg_table_median<-rbind(gg_table_median, cbind(input.group.x[x], input.group.y[y], median(as.numeric(as.character(gg_table_sub[,3])))))
      
    }
  }
  gg_table_median<-data.frame(gg_table_median)
  gg_table_median[,3]<-as.numeric(as.character(gg_table_median[,3]))
  
  
  #normalize
  gg_table_median_norm<-c()
  input.group.x<-unique(as.character(gg_table[,1]))
  input.group.y<-unique(as.character(gg_table[,2]))
  
  
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  
  if (input.norm == "y")
    for (y in 1:length(input.group.y)){
      gg_table_median_sub<-subset(gg_table_median, gg_table_median[,2] == input.group.y[y])
      norm_value<- range01(as.numeric(as.character(gg_table_median_sub[,3])))
      gg_table_median_sub[,3]<-norm_value
      gg_table_median_norm<-rbind(gg_table_median_norm, gg_table_median_sub)
    }
  
  if (input.norm == "x")
    for (x in 1:length(input.group.x)){
      gg_table_median_sub<-subset(gg_table_median, gg_table_median[,1] == input.group.x[x])
      norm_value<- range01(as.numeric(as.character(gg_table_median_sub[,3])))
      gg_table_median_sub[,3]<-norm_value
      gg_table_median_norm<-rbind(gg_table_median_norm, gg_table_median_sub)
    }
  
  if (input.norm == "na") gg_table_median_norm<-gg_table_median
  
  
  gg_table_median_norm<-data.frame(gg_table_median_norm)
  gg_table_median_norm[,3]<-as.numeric(as.character(gg_table_median_norm[,3]))
  
  gg_table_median_norm_tibble <- gg_table_median_norm %>% as_tibble()
  names(gg_table_median_norm_tibble) <- c('group', 'pathway', 'score_median_norm_max_min')
  
  gg_table_median_tibble <- gg_table_median %>% as_tibble()
  names(gg_table_median_tibble) <- c('group', 'pathway', 'score_median')
  
  out.df <- gg_table_median_tibble %>% left_join(gg_table_median_norm_tibble)
  return(out.df)
  
}

#rds <- readRDS('Bcells_KEGG.rds')
#rds
#a <- get_table_from_scMe_obj(rds, "diabetes")
#a