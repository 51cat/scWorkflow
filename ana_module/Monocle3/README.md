| 参数                | 是否为必须参数 | 默认值 | 描述                                                         |
| ------------------- | -------------- | ------ | ------------------------------------------------------------ |
| `--rds`             | 是             | 无     | seurat对象的rds路径                                          |
| `--outdir`          | 是             | 无     | 结果输出路径                                                 |
| `--celltype_col`    | 是             | 无     | 细胞类型所在列的列名，必须在rds@meta.data中存在              |
| `--group_col`       | 是             | 无     | 分组所在列的列名，必须在rds@meta.data中存在                  |
| `--root_celltype`   | 是             | 无     | 起点细胞类型，必须在`--celltype_col`列中存在                 |
| `--cds`             | 可选参数       | None   | monocle3分析结果的cds路径，加入此参数可跳过轨迹分析部分      |
| `--use_seurat_umap` | 可选参数       | True   | 是否使用`Seurat`的聚类结果代替`Monocle3`的聚类结果           |
| `--reduction_name`  | 可选参数       | umap   | 使用`Seurat`的聚类结果代替`Monocle3`的聚类结果时，umap结果在seurat中的名称 |
| `--target_genes`    | 可选参数       | None   | 使用指定的gene list进行分析，输入是一个表格（csv），表格中必须包含gene列，自动识别gene列中的基因进行分析 |
| `--spec`            | 可选参数       | human  | 富集分析使用的物种，支持human和mouse                         |
| `--gcpal`           | 可选参数       | npg    | 分组展示使用的配色                                           |
| `--ccpal`           | 可选参数       | Paired | 细胞类型展示使用的配色                                       |
