library(clusterProfiler)
library(tidyverse)
library(DBI)
library(RSQLite)

date <- "2024.1"

URL <- str_glue("https://data.broadinstitute.org/gsea-msigdb/msigdb/release/{date}.Mm/msigdb_v{date}.Mm.db.zip")

# mk db (human)

local_zip_path <- str_glue("./msigdb_v{date}.Mm.db.zip")
#download.file(URL, local_zip_path)
print("download success!")
unzip(local_zip_path, exdir = "./")

con <- DBI::dbConnect(RSQLite::SQLite(), dbname = str_glue("./msigdb_v{date}.Mm.db"))
DBI::dbListTables(con)

geneset_db <- dplyr::tbl(con, 'gene_set')                                              # standard_name, collection_name
details_db <- dplyr::tbl(con, 'gene_set_details')                                      # description_brief, description_full
geneset_genesymbol_db <- dplyr::tbl(con, 'gene_set_gene_symbol')                       # meat and potatoes
genesymbol_db <- dplyr::tbl(con, 'gene_symbol')                                        # mapping from ids to gene symbols
collection_db <- dplyr::tbl(con, 'collection') %>% dplyr::select(collection_name, full_name)

# join tables
msigdb <- geneset_db %>%
  dplyr::left_join(details_db, by = c('id' = 'gene_set_id')) %>%
  dplyr::left_join(collection_db, by = 'collection_name') %>%
  dplyr::left_join(geneset_genesymbol_db, by = c('id' = 'gene_set_id')) %>%
  dplyr::left_join(genesymbol_db, by = c('gene_symbol_id' = 'id')) %>%
  dplyr::select(collection = collection_name, subcollection = full_name, geneset = standard_name, description = description_brief, symbol) %>%
  dplyr::as_tibble() 


DBI::dbDisconnect(con)

collection.all <- msigdb %>% pull(collection) %>% unique

for (collection.use in collection.all) {
    msigdb.h <- msigdb %>% 
    dplyr::filter(collection==collection.use) %>% 
    dplyr::select(c("geneset", "symbol"))
    msigdb.h$geneset <- factor(msigdb.h$geneset)
    msigdb.h <- msigdb.h %>% 
    dplyr::group_split(geneset, .keep = F) %>%
    purrr::map( ~.x %>% dplyr::pull(symbol) %>% unique(.)) %>%
    purrr::set_names(levels(msigdb.h$geneset))

    saveRDS(msigdb.h, str_glue("./{collection.use}.rds"))
    print(str_glue("Finish {collection.use}"))
}

