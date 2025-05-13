| 参数                    | 是否为必须参数 | 默认值          | 描述                                                         |
| ----------------------- | -------------- | --------------- | ------------------------------------------------------------ |
| `--ref_h5ad`            | 可选参数       | `control.h5ad`  | 参考数据的 `.h5ad` 文件路径（通常为对照组）                  |
| `--query_h5ad`          | 可选参数       | `case.h5ad`     | 查询数据的 `.h5ad` 文件路径（通常为实验组）                  |
| `--ref_name`            | 可选参数       | `control`       | 参考数据的名称，用于输出文件命名                             |
| `--query_name`          | 可选参数       | `case`          | 查询数据的名称，用于输出文件命名                             |
| `--ann_col`             | 可选参数       | 无              | 注释列名，表示细胞类型等信息，需在 `obs` 中存在              |
| `--pse_df`              | 是             | 空字符串        | 包含伪时间信息的 CSV 文件路径（index 为细胞名，需含 time 列） |
| `--target_gene_file`    | 可选参数       | `None`          | 包含要展示的目标基因名的文本文件，每行一个基因               |
| `--compare_file`        | 可选参数       | `None`          | 若需批量对比多个样本对，提供包含四列的 TSV：ref_name、ref_h5ad、query_name、query_h5ad（无需列名） |
| `--find_cluster_method` | 可选参数       | `distance,auto` | 聚类方法，格式为 `方式,参数`，例如 `cluster,5` 表示指定聚类数为 5 |
| `--outdir`              | 是             | `./test_out2/`  | 输出目录                                                     |
| `--n_bins`              | 可选参数       | `14`            | 伪时间轴的分箱数，用于对齐分析                               |
| `--is_scale_pse`        | 可选参数       | `True`          | 是否对伪时间进行归一化                                       |
| `--is_downsample`       | 可选参数       | `False`         | 是否下采样基因数量（1500个基因）                             |

- `--compare_file`： 示例

  - 提供包含四列的 TSV：ref_name、ref_h5ad、query_name、query_h5ad（无需列名）

  ```
  control1	/path/to/control1.h5ad	case1	/path/to/case1.h5ad
  control2	/path/to/control2.h5ad	case2	/path/to/case2.h5ad
  control3	/path/to/control3.h5ad	case3	/path/to/case3.h5ad
  ```

  