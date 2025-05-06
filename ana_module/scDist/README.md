| 参数               | 是否为必须参数 | 默认值     | 描述                                            |
| ------------------ | -------------- | ---------- | ----------------------------------------------- |
| `--rds`            | 是             | 无         | seurat对象的rds路径                             |
| `--outdir`         | 是             | 无         | 结果输出路径                                    |
| `--celltype_col`   | 是             | 无         | 细胞类型所在列的列名，必须在rds@meta.data中存在 |
| `--group_col`      | 是             | 无         | 分组所在列的列名，必须在rds@meta.data中存在     |
| `--compare_str`    | 是             | 无         | 比较方式, 逗号分隔: AvsB,AvsC,BvsC              |
| `--random_effects` | 否             | orig_ident | 批次信息所在列，必须在rds@meta.data中存在       |

