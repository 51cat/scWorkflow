library(clusterProfiler)
library(tidyverse)
df1 <- clusterProfiler::read.gmt("REACTOME_metabolism.gmt")


df1$term <- str_replace_all(df1$term, "/", "") %>% str_replace_all( " ", "_")

reach <- split(df1$gene, df1$term)

df2<- clusterProfiler::read.gmt("KEGG_metabolism_nc.gmt")


df2$term <- str_replace_all(df2$term, "/", "") %>% str_replace_all( " ", "_")

kegg <- split(df2$gene, df2$term)

saveRDS(reach, "REACTOME_metabolism.rds")
saveRDS(kegg, "KEGG_metabolism_nc.rds")