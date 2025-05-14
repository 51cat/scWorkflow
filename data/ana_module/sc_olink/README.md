
| 参数名              | 是否为必须参数 | 默认值 | 说明                                                         |
| ------------------- | -------------- | ------ | ------------------------------------------------------------ |
| `--npx`             | 是             | 无     | 输入的NPX数据文件路径                                        |
| `--rds`             | 是             | 无     | 输入的Seurat对象rds文件路径                                  |
| `--tf_rds`          | 可选           | None   | 输入的TF数据rds路径，如果方法是TF时使用                      |
| `--method`          | 是             | 无     | 选择分析方法，可选值为`cellchat`或`TF`                       |
| `--target_pathway`  | 可选           | None   | CellChat分析中使用的目标通路路径，列名为gene                 |
| `--target_TF`       | 可选           | None   | TF分析中使用的目标转录因子，列名为gene                       |
| `--target_olink`    | 可选           | None   | 分析中使用的目标Olink ID文件路径，列名为gene                 |
| `--group_sort`      | 可选           | None   | 分组展示的顺序，多个分组之间用逗号分隔                       |
| `--gcpal`           | 可选           | npg    | 可选配色方案，默认为npg，指定不同的配色方案展示分组          |
| `--show_gene_name`  | 可选           | False  | 是否显示基因名，默认为`False`                                |
| `--outdir`          | 是             | 无     | 结果输出目录路径                                             |
| `--group_col`       | 是             | 无     | 分组所在列的列名，必须存在于rds@meta.data中                  |
| `--sample_col`      | 是             | 无     | 样本/重复所在列的列名，必须存在于rds@meta.data中             |
| `--celltype_col`    | 是             | 无     | 细胞类型所在列的列名，必须存在于rds@meta.data中              |
| `--olink_id_file`   | 可选           | None   | Olink ID文件路径，用于数据可视化                             |
| `--drop_group_info` | 可选           | False  | 是否丢弃分组信息，默认值为`False`，如果为`True`则不会保留分组信息 |

- 需要增加 `--no_env True`参数

- NPX数据格式要求

  - 示例

  ```
  SampleID	UniProt	NPX	sample_name	inx	MAP_name	group	sample
  KBMC23072905BP-01	P13236	7.85	LHC240905128	1	MAP24091440H07	Ctrl	Ctrl_16
  KBMC23072905BP-01	P01137	8.87	LHC240905128	1	MAP24091440H07	Ctrl	Ctrl_16
  KBMC23072905BP-01	O43508	7.41	LHC240905128	1	MAP24091440H07	Ctrl	Ctrl_16
  KBMC23072905BP-01	P07585	-0.16	LHC240905128	1	MAP24091440H07	Ctrl	Ctrl_16
  ```

  1. 必须包含group, sample, NPX, UniPro四列，其他列可以是rds@meta.data中的其他信息
  2. group和sample列的信息分别为`--group_col`和`--sample_col`的信息

- olink_id_file: 列名需要一致

  - 示例

  ```
  UniProt olink_gene
  P13236  CCL4
  P01137  LAP TGF-beta-1
  O43508  TWEAK
  ```

- `--tf_rds`:  转录因子分析流程会输出此rds