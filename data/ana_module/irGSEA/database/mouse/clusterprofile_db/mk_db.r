library(clusterProfiler)
library(tidyverse)

# KEGG
kegg.spec <- "mmu"
# GO
OrgDb <- "org.Mm.eg.db"

GO.ont <- c("BP", "MF", "CC", "ALL")



kk <- clusterProfiler::gson_KEGG(species = kegg.spec)
gson::write.gson(kk, file = str_glue("./KEGG_{kegg.spec}.gson"))

# read gson file
kk2 <- gson::read.gson(str_glue("./KEGG_{kegg.spec}.gson"))
# Convert to a data frame
kegg.list <- dplyr::left_join(kk2@gsid2name,
                              kk2@gsid2gene,
                              by = "gsid")

gene_name <- clusterProfiler::bitr(kegg.list$gene, 
                                   fromType = "ENTREZID", 
                                   toType = "SYMBOL", 
                                   OrgDb = OrgDb)
kegg.list <- dplyr::full_join(kegg.list,
                              gene_name,
                              by = c("gene"="ENTREZID"))

kegg.list <- kegg.list[complete.cases(kegg.list[, c("gene", "SYMBOL")]), ]

kegg.list$name <- factor(kegg.list$name)
kegg.list <- kegg.list %>% 
  dplyr::group_split(name, .keep = F) %>%
  purrr::map( ~.x %>% dplyr::pull(SYMBOL) %>% unique(.)) %>%
  purrr::set_names(levels(kegg.list$name))

saveRDS(kegg.list, str_glue("./kegg_{kegg.spec}.rds"))
print(str_glue("kegg_{kegg.spec} finish!"))

# GO
for (ont.name in GO.ont) {
    ont.use <- ont.name
    go <- clusterProfiler::gson_GO(OrgDb = OrgDb, ont = ont.use)
    gson::write.gson(go, file = str_glue("./go_{ont.use}.gson"))

    # read gson file
    go2 <- gson::read.gson(str_glue("./go_{ont.use}.gson"))

    # Convert to a data frame
    go.list <- dplyr::left_join(go2@gsid2name,
                                go2@gsid2gene,
                                by = "gsid")

    go.list <- dplyr::full_join(go.list,
                                go2@gene2name,
                                by = c("gene"="ENTREZID"))
    # remove NA value if exist
    go.list <- go.list[complete.cases(go.list[, c("gene", "SYMBOL")]), ]

    go.list$name <- factor(go.list$name)
    go.list <- go.list %>% 
    dplyr::group_split(name, .keep = F) %>%
    purrr::map( ~.x %>% dplyr::pull(SYMBOL) %>% unique(.)) %>%
    purrr::set_names(levels(go.list$name))

    saveRDS(go.list, str_glue("./GO_{ont.name}.rds"))
    print(str_glue("GO_{ont.name} finish!"))
}
