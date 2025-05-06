| 参数               | 是否为必须参数 | 默认值      | 描述                                            |
| ------------------ | -------------- | ----------- | ----------------------------------------------- |
| `--rds`            | 是             | 无          | seurat对象的rds路径                             |
| `--outdir`         | 是             | 无          | 结果输出路径                                    |
| `--celltype_col`   | 是             | 无          | 细胞类型所在列的列名，必须在rds@meta.data中存在 |
| `--compare_str`    | 是             | 无          | 比较方式, 逗号分隔: AvsB,AvsC,BvsC              |
| `--tre_col`        | 是             | 无          | 处理所在列的列名，必须在rds@meta.data中存在     |
| `--rep_col`        | 可选参数       | orig.ident  | 重复所在列的列名，必须在rds@meta.data中存在     |
| `--compare_method` | 可选参数       | wilcox.test | 组间比较方法，wilcox.test 或 kruskal.test       |
| `--gcpal`          | 可选参数       | npg         | 不同分组的配色方案                              |
