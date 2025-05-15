| 参数          | 是否为必须参数 | 默认值          | 描述                                                         |
| ------------- | -------------- | --------------- | ------------------------------------------------------------ |
| `--input`     | 是             | 无              | seurat对象的rds路径                                          |
| `--outdir`    | 是             | 无              | 结果输出路径                                                 |
| `--split_col` | 可选参数       | None            | 根据哪列拆分数据，传入此参数会自动拆分产生h5ad               |
| `--keep`      | 可选参数       | None            | `--split_col`指定的列中保留哪些数据，例如，cell1,cell2,cell3 |
| `--keep_str`  | 可选参数       | None            | 以一个文件的形式输入，文件内容只有一行，内容是保留哪些数据的表达式, 可以搭配`--split_col`和`--keep`实现任意部分的数据提取和格式转换和拆分 |
| `--fr`        | 可选参数       | seurat          | 输入数据格式(seurat/h5ad)                                    |
| `--to`        | 可选参数       | h5ad            | 输出数据格式(h5ad/mtx/count_h5ad/normalized_h5ad)            |
| `--pyexec`    | 可选参数       | /usr/bin/python | python解释器路径                                             |

- 需要增加 `--no_env True`参数

- 支持的格式转换：

  1. seurat -> h5ad:  

     ```python
     --fr seurat
     --to h5ad
     ```

  2. seurat  -> mtx:

     ```python
     --fr seurat
     --to mtx
     ```

  3. seurat  -> count_h5ad:

     ```python
     # 一个只包含count矩阵的h5ad文件
     --fr seurat
     --to count_h5ad
     ```

  4. seurat -> normailze_h5ad

     ```python
     # 一个包含count和normalize矩阵的h5ad文件
     --fr seurat
     --to normalize_h5ad
     ```

  5. h5ad -> seurat 

     ```python
     --fr h5ad
     --to seurat
     ```

     

- `--keep_str`例子

1. 根据`nFeature_RNA`/`percent.mt`/`ann`筛选数据

   ```
   nFeature_RNA > 200 & percent.mt < 5 & ann == 'use'
   ```

2. 根据`celltype`等选择

   ```
   subcelltype %in% c('M0', 'M1') & orig.ident %in% c('library1','library2','library3' )
   ```

3. 筛选特定`barcode `(`rds`的`meta.data`中需要有一列信息储存了barcode)

   ```
   barcode_colnames %in% read_tsv('/path/to/barcode_list.tsv')$barcode
   ```

