分析参数：

- 表格
  - foldchange：0
  - p：不设阈值
  - min.pct = 0.25

- 柱状图/散点图：
  - foldchange：0.25
  - p：0.05
  - min.pct = 0.25

| 参数             | 是否为必须参数 | 默认值 | 描述                                            |
| ---------------- | -------------- | ------ | ----------------------------------------------- |
| `--rds`          | 是             | 无     | seurat对象的rds路径                             |
| `--outdir`       | 是             | 无     | 结果输出路径                                    |
| `--celltype_col` | 是             | 无     | 细胞类型所在列的列名，必须在rds@meta.data中存在 |
| `--group_col`    | 是             | 无     | 分组所在列的列名，必须在rds@meta.data中存在     |
| `--compare_str`  | 是             | 无     | 比较方式, 逗号分隔: AvsB,AvsC,BvsC              |
| `--default_data` | 可选参数       | None   | 使用哪个数据分析，默认用RNA数据分析             |

