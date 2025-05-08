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

parser <- ArgumentParser(description = "Wrapper script to run score.r and compare.r")

# 运行 score.r 需要的参数
parser$add_argument("--rds", type = "character", help = "Path to RDS file")
parser$add_argument("--group_col", type = "character", help = "Group column name")
parser$add_argument("--celltype_col", type = "character", help = "Cell type column name")
parser$add_argument("--outdir", type = "character", help = "Output directory")
parser$add_argument("--target_cell", type = "character", default = "all", help = "Target cell type")
parser$add_argument("--geneset", type = "character", help = "Path to geneset .rds or .gmt")
parser$add_argument("--compare_str", type = "character", help = "Comparison string (e.g., casevscontrol)")
parser$add_argument("--methods", type = "character", default = "AUCell,UCell,singscore,ssgsea", help = "Methods to use")
parser$add_argument("--compare_method", type = "character", default = "within,between", help = "Comparison method")
parser$add_argument("--minGSSize", type = "integer", default = 1, help = "Minimum gene set size")
parser$add_argument("--maxGSSize", type = "integer", default = 9999, help = "Maximum gene set size")

# compare.r 专用参数
parser$add_argument("--group_sort", type = "character", default = NULL, help = "Group sort order (e.g., control,case)")
parser$add_argument("--gcpal", type = "character", default = NULL, help = "Group color palette")
parser$add_argument("--ccpal", type = "character", default = NULL, help = "Cell type color palette")

# 控制步骤
parser$add_argument("--step", help = "Which step(s) to run: score, compare, or both (comma-separated)", default = 'score,compare')

args <- parser$parse_args()

score_script <- str_glue("{dirnameScript()}/score.r")
compare_script <- str_glue("{dirnameScript()}/compare.r")

steps <- str_split(args$step, ",")[[1]]

# Step 1: Run score.r
if ('score' %in% steps) {
  cat("Running score.r ...\n")
  score_cmd <- str_glue("Rscript {score_script} ",
                        "--rds \"{args$rds}\" ",
                        "--group_col {args$group_col} ",
                        "--celltype_col {args$celltype_col} ",
                        "--outdir \"{args$outdir}\" ",
                        "--target_cell \"{args$target_cell}\" ",
                        "--geneset \"{args$geneset}\" ",
                        "--compare_str \"{args$compare_str}\" ",
                        "--methods \"{args$methods}\" ",
                        "--compare_method \"{args$compare_method}\" ",
                        "--minGSSize {args$minGSSize} ",
                        "--maxGSSize {args$maxGSSize}")
  system(score_cmd)
}

# Step 2: Run compare.r
if ('compare' %in% steps) {
  cat("Running compare.r ...\n")
  compare_cmd <- str_glue("Rscript {compare_script} ",
                          "--score_step_outdir \"{args$outdir}\" ",
                          "--outdir \"{args$outdir}\"")

  if (!is.null(args$group_sort)) {
    compare_cmd <- str_glue("{compare_cmd} --group_sort \"{args$group_sort}\"")
  }
  if (!is.null(args$gcpal)) {
    compare_cmd <- str_glue("{compare_cmd} --gcpal \"{args$gcpal}\"")
  }
  if (!is.null(args$ccpal)) {
    compare_cmd <- str_glue("{compare_cmd} --ccpal \"{args$ccpal}\"")
  }

  system(compare_cmd)
}
