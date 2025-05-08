library(argparse)
library(stringr)

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

print(dirnameScript())

COLOR_DB <- readRDS(str_glue("{dirnameScript()}/colorDB.rds"))

updateDB <- function(name, colors, color_type) {
    if (name %in% names(COLOR_DB)) {
        print('颜色名称重复!!')
        return(-1)
    }

    if (!color_type %in% c('c', 'd')) {
        print('颜色类型错误，请参考文档！')
        return(-1)
    }

    if (color_type == 'c') {
        color_type <- 'continuous'
    }

    if (color_type == 'd') {
        color_type <- 'discrete'
    }


    COLOR_DB[[name]] <- colors
    attr(COLOR_DB[[name]], 'type') <- color_type
    return(COLOR_DB)
}



parser <- ArgumentParser(description = "")

parser$add_argument("--colors", type = "character", required = TRUE, help = "")
parser$add_argument("--name", type = "character", required = TRUE, help = "")
parser$add_argument("--color_type", type = "character", required = FALSE, help = "", default = 'd')


args <- parser$parse_args()

colors.use <- stringr::str_split(args$colors, ",")[[1]]

DB.new <- updateDB(
    args$name,
    colors.use,
    args$color_type
)

search_dir <- dirname(dirnameScript())

print(search_dir)

files <- list.files(
  path = search_dir,
  pattern = "^colorDB\\.rds$", 
  recursive = TRUE,
  full.names = TRUE
)

print(files)

for (f in files) {
    saveRDS(DB.new, f)
}