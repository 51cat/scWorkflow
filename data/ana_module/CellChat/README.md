| 参数                 | 是否为必须参数 | 默认值        | 描述                                                         |
| -------------------- | -------------- | ------------- | ------------------------------------------------------------ |
| `--rds`              | 是             | 无            | seurat对象的rds路径                                          |
| `--outdir`           | 是             | 无            | 结果输出路径                                                 |
| `--celltype_col`     | 是             | 无            | 细胞类型所在列的列名，必须在rds@meta.data中存在              |
| `--group_col`        | 是             | 无            | 分组所在列的列名，必须在rds@meta.data中存在                  |
| `--compare_str`      | 可选参数       | 无            | 比较方式, 逗号分隔: AvsB,AvsC,BvsC                           |
| `--compare_method`   | 可选参数       | pairs         | **pairs**: 两两比较，此时必须提供参数`--compare_str`; **multi**: 多组比较分析，默认所有组放在一起分析 |
| `--cellchat_rds_dir` | 可选参数       | None          | 流程输出的rds所在目录，给定此参数可以跳过前处理过程          |
| `--spec`             | 可选参数       | human         | 富集分析使用的物种，支持human和mouse                         |
| `--gcpal`            | 可选参数       | npg           | 分组展示使用的配色                                           |
| `--ccpal`            | 可选参数       | Paired        | 细胞类型展示使用的配色                                       |
| `--group_sort`       | 可选参数       | None          | 当compare_method为multi时，每个组别在输出图片的排列顺序，逗号分隔 |
| `--py_exec`          | 可选参数       | /usr/bin/py39 | python路径                                                   |

实例：

```shell
singularity exec CellChat.sif Rscript start.r \
    --rds "/public/home/zhliu/singlecell/test_data/obj_v5.rds" \
    --outdir ./outs \
    --celltype_col subcelltype \
    --group_col split \
    --compare_method pairs \
    --compare_str casevscontrol \
    --cellchat_rds_dir "/public/home/zhliu/PIPLINE/singlecell/Cellchat/demo/demo_results_new666/innput_obj/"
# multi
singularity exec CellChat.sif Rscript start.r \
    --rds "/public/home/zhliu/singlecell/test_data/obj_v5.rds" \
    --outdir ./outs \
    --celltype_col subcelltype \
    --group_col split \
    --compare_method multi \
    --cellchat_rds_dir "/public/home/zhliu/PIPLINE/singlecell/Cellchat/demo/demo_results_new666/innput_obj/"
```

