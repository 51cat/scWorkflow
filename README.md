# 单细胞转录组分析流程

单细胞转录组分析流程：

1. 目前包含了22个分析流程，均支持单样本/多样本/对比分析
2. 基于`singularity`，方便迁移
3. 统一的调用接口
4. 支持`slurm`

## 安装

1. clone 本仓库

   ```shell
   git clone https://github.com/51cat/scWorkflow.git 
   ```

2. 创建分析环境

   ```shell
   cd scWorkflow
   conda create -n myenv python=3.12
   conda activate myenv
   conda install -y --file conda_pkgs.txt
   ```

3. 安装软件

   ```shell
   pip install .
   ```

4. 测试

   ```shell
   scwf --help
   ```

   

## 安装流程

1. 分析流程安装

   ```shell
   scwf_tk install-module
   ```

2. 安装分析环境（非常耗时）

   ```shell
   scwf_tk mk-config
   scwf_tk install-env
   ```

## 测试

查看分析流程

```shell
scwf ls
```

查看分析环境

```shell
scwf lsenv
```

## 查看帮助信息

```shell
scwf help [流程名]
# eg
scwf help CellChat
```

## 执行流程

使用`run`子命令执行分析流程

```shell
scwf run --module CellChat [cellchat流程的选项参数, 使用scwf help CellChat查看]
```

- 一个例子

  ```shell
   scwf run --module scMetabolism \
       --rds /path/to/my.rds \
       --outdir ./debug_tmp/03/outs/ \
       --celltype_col subcelltype  \
       --group_col split \
       --compare_str casevscontrol
  ```

- 提交`slurm`

  ```shell
   scwf run --module scMetabolism \
       --rds /path/to/my.rds \
       --outdir ./debug_tmp/03/outs/ \
       --celltype_col subcelltype  \
       --group_col split \
       --compare_str casevscontrol \
       --exec_method slurm \
       --cpu 1 \
       --mem 30G \
       --wn my_work_name
  ```

## 查看slurm脚本执行状态

```shell
scwf task-stat
```

