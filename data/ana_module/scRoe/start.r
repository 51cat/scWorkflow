library(Startrac)
library(ggplot2)
library(tictoc)
library(ggpubr)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(tidyverse)
library(sscVis)
library(argparse)
print(sessionInfo())

mk.outdir <- function(dir) {
    if (!dir.exists(dir)){
        dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    }
    return(dir)
}


parser <- ArgumentParser(description = "")

parser$add_argument("--rds", type = "character", default = "",
                    help = "Path to the RDS file")
parser$add_argument("--group_col", type = "character", default = "",
                    help = "group column name")
parser$add_argument("--batch_col", type = "character", 
                    default = "orig.ident",
                    help = "")
parser$add_argument("--celltype_col", type = "character", 
                    default = "",
                    help = "")

parser$add_argument("--outdir", type = "character", default = "",
                    help = "Output directory")


args <- parser$parse_args()


rds <- readRDS(args$rds)
in_data <- rds@meta.data

group_col <- args$group_col
patient_col <- args$batch_col
celltype_col <- args$celltype_col


Roe <- calTissueDist(in_data,
       byPatient = F,
         colname.cluster = celltype_col,
         colname.patient = patient_col,
         colname.tissue = group_col,
         method = "chisq", 
         min.rowSum = 0) 

col_fun <- colorRamp2(c(min(Roe, na.rm = TRUE), 1, max(Roe, na.rm = TRUE)), 
                      c("blue", "white", "red"))

ncell <- rds@meta.data[[celltype_col]] %>% unique %>% length()
ngroup <- rds@meta.data[[group_col]] %>% unique %>% length()
w <- min(ngroup*1.5, 20)
h <- min(ncell*0.3, 40)

out.p <- mk.outdir(args$outdir)

Roe.df <- Roe %>% 
  as.data.frame() %>% 
  mutate(`Ro/e` = case_when(
    Freq > 1 ~ "+++",
    Freq > 0.8 ~ "++",
    Freq > 0 ~ "+/-",
    .default='-'
  )) 
names(Roe.df) <- c(celltype_col,group_col, 'Ro/e_Index',  "Ro/e")

Roe.df %>% write_csv(str_glue("{out.p}/{group_col}_{celltype_col}_Roe.csv")) 


# 横向
df.mtx <- Roe.df %>% 
  select(-`Ro/e`) %>% 
  spread(key = group_col, value = "Ro/e_Index") %>% 
  column_to_rownames(var = celltype_col) %>% t

df.display <- Roe.df %>% 
  select(-`Ro/e_Index`) %>% 
  spread(key = group_col, value = "Ro/e") %>% 
  column_to_rownames(var = celltype_col) %>% t

print(df.mtx )
w.h <- min(ncell*0.5, 20)
h.h <- min(ngroup*1, 40)

# 颜色映射固定 white 对应 1
col_fun <- colorRamp2(
  c(min(Roe, na.rm = TRUE), 1, max(Roe, na.rm = TRUE)),
  c("blue", "white", "red")
)

# ========== 纵向热图：数值显示 ==========
pdf(str_glue("{out.p}/{group_col}_{celltype_col}_Roe_V.pdf"), width = w, height = h) 
Heatmap(as.matrix(Roe),
        show_heatmap_legend = TRUE, 
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        row_names_side = 'right', 
        show_column_names = TRUE,
        show_row_names = TRUE,
        col = col_fun,
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        heatmap_legend_param = list(
          title = "Ro/e Index",
          at = c(min(Roe, na.rm = TRUE), max(Roe, na.rm = TRUE)),
          labels = c("Low", "High")
        ),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", Roe[i, j]), x, y, gp = gpar(fontsize = 8, col = "black"))
        }
)
dev.off()

# ========== 纵向热图：符号显示（+++、++、+/- 等）==========
pdf(str_glue("{out.p}/{group_col}_{celltype_col}_Roe_label_V.pdf"), width = w, height = h) 
Heatmap(as.matrix(Roe),
        show_heatmap_legend = TRUE,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_names_side = 'right',
        show_column_names = TRUE,
        show_row_names = TRUE,
        col = col_fun,
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        heatmap_legend_param = list(
          title = "Ro/e",
          at = c(min(Roe, na.rm = TRUE), max(Roe, na.rm = TRUE)),
          labels = c("Low", "High")
        ),
        cell_fun = function(j, i, x, y, width, height, fill) {
          value <- Roe[i, j]
          label <- ifelse(value > 1, "+++",
                   ifelse(value > 0.8, "++",
                   ifelse(value >= 0.2, "+",
                   ifelse(value > 0, "+/-", "-"))))
          grid.text(label, x, y, gp = gpar(fontsize = 8, col = "black"))
        }
)
dev.off()


# ========== 横向热图：数值显示 ==========
pdf(str_glue("{out.p}/{group_col}_{celltype_col}_Roe_H_ComplexHeatmap.pdf"), width = w.h, height = h.h)
Heatmap(t(as.matrix(Roe)),  # 转置，横轴为 celltype
        show_heatmap_legend = TRUE, 
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        row_names_side = 'right', 
        show_column_names = TRUE,
        show_row_names = TRUE,
        col = col_fun,
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        heatmap_legend_param = list(
          title = "Ro/e Index",
          at = c(min(Roe, na.rm = TRUE), max(Roe, na.rm = TRUE)),
          labels = c("Low", "High")
        ),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", t(Roe)[i, j]), x, y, gp = gpar(fontsize = 8, col = "black"))
        }
)
dev.off()

# ========== 横向热图：符号显示（+++、++、+/- 等）==========
pdf(str_glue("{out.p}/{group_col}_{celltype_col}_Roe_H_label_ComplexHeatmap.pdf"), width = w.h, height = h.h)
Heatmap(t(as.matrix(Roe)), 
        show_heatmap_legend = TRUE,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_names_side = 'right',
        show_column_names = TRUE,
        show_row_names = TRUE,
        col = col_fun,
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        heatmap_legend_param = list(
          title = "Ro/e",
          at = c(min(Roe, na.rm = TRUE), max(Roe, na.rm = TRUE)),
          labels = c("Low", "High")
        ),
        cell_fun = function(j, i, x, y, width, height, fill) {
          value <- t(Roe)[i, j]
          label <- ifelse(value > 1, "+++",
                   ifelse(value > 0.8, "++",
                   ifelse(value >= 0.2, "+",
                   ifelse(value > 0, "+/-", "-"))))
          grid.text(label, x, y, gp = gpar(fontsize = 8, col = "black"))
        }
)
dev.off()
