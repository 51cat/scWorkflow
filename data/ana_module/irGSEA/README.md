| 参数               | 是否为必须参数 | 默认值                        | 描述                                                         |
| ------------------ | -------------- | ----------------------------- | ------------------------------------------------------------ |
| `--rds`            | 是             | 无                            | Seurat 对象的 RDS 路径                                       |
| `--outdir`         | 是             | 无                            | 结果输出路径                                                 |
| `--celltype_col`   | 是             | 无                            | 细胞类型所在列名，必须存在于 `rds@meta.data` 中              |
| `--group_col`      | 是             | 无                            | 分组所在列名，必须存在于 `rds@meta.data` 中                  |
| `--geneset`        | 是             | 无                            | 富集分析使用的基因集路径，可以使用自带数据库，详情见下，支持 `.rds` 或 `.gmt` 格式 |
| `--compare_str`    | 是             | 无                            | 分组比较字符串，例如 `casevscontrol`，用于指定比较对，多个比较对用逗号分开 |
| `--target_cell`    | 可选参数       | all                           | 指定用于分析的细胞类型，多个用逗号分隔，all 表示全部细胞类型 |
| `--methods`        | 可选参数       | AUCell,UCell,singscore,ssgsea | 使用的打分方法，多个方法用逗号分隔                           |
| `--compare_method` | 可选参数       | within,between                | 差异分析的方法，支持组内（within）与组间（between）          |
| `--minGSSize`      | 可选参数       | 1                             | 基因集最小基因数                                             |
| `--maxGSSize`      | 可选参数       | 500                           | 基因集最大基因数                                             |
| `--group_sort`     | 可选参数       | 无                            | 指定分组顺序（例如 `control,case`），影响绘图顺序            |
| `--gcpal`          | 可选参数       | npg                           | 分组颜色配色方案，用于 `compare.r` 中的绘图                  |
| `--ccpal`          | 可选参数       | Paired                        | 细胞类型颜色配色方案，用于 `compare.r` 中的绘图              |
| `--step`           | 可选参数       | score,compare                 | 指定执行的步骤，支持 `score`、`compare` 或两者组合（用逗号分隔） |

- 使用指定的基因集以及scMetabolism的基因集时，建议设定：
  - `--minGSSize 1`
  - `--maxGSSize 99999`
- 数据库信息

1. 使用`--geneset`可以指定`gmt`文件或者rds文件，rds文件的对象需要是一个`list`每个元素是一个通路对应的基因，需要是`gene symbol`

2. 软件自带数据库human和mouse，直接传入`--geneset`即可，例如 `--geneset clusterprofile_db/GO_BP` 

human

- clusterprofile_db/GO_ALL
- clusterprofile_db/GO_BP
- clusterprofile_db/GO_CC
- clusterprofile_db/GO_MF
- clusterprofile_db/kegg_hsa
- scMetabolism/KEG_metabolism_nc
- scMetabolism/REACTOME_metabolism
- msigdb/C1
- msigdb/C2_ALL
- msigdb/C2:CGP
- msigdb/C2:CP:BIOCARTA
- msigdb/C2:CP:KEG_LEGACY
- msigdb/C2:CP:KEG_MEDICUS
- msigdb/C2:CP:PID
- msigdb/C2:CP:REACTOME
- msigdb/C2:CP:WIKIPATHWAYS
- msigdb/C3:MIR:MIRDB
- msigdb/C3:MIR:MIR_LEGACY
- msigdb/C3:TFT:GTRD
- msigdb/C3:TFT:TFT_LEGACY
- msigdb/C4:3CA
- msigdb/C4:CGN
- msigdb/C4:CM
- msigdb/C5_ALL
- msigdb/C5:GO:BP
- msigdb/C5:GO:CC
- msigdb/C5:GO:MF
- msigdb/C5:HPO
- msigdb/C6
- msigdb/C7:IMMUNESIGDB
- msigdb/C7:VAX
- msigdb/C8
- msigdb/H

mouse

- clusterprofile_db/GO_ALL
- clusterprofile_db/GO_BP
- clusterprofile_db/GO_CC
- clusterprofile_db/GO_MF
- clusterprofile_db/kegg_mmu
- msigdb/M1
- msigdb/M2:CGP
- msigdb/M2:CP:BIOCARTA
- msigdb/M2:CP:REACTOME
- msigdb/M2:CP:WIKIPATHWAYS
- msigdb/M3:GTRD
- msigdb/M3:MIRDB
- msigdb/M5:GO:BP
- msigdb/M5:GO:CC
- msigdb/M5:GO:MF
- msigdb/M5:MPT
- msigdb/M8
- msigdb/MH