import pandas as pd
from scipy.stats import pearsonr, spearmanr
import tqdm
import argparse
import os
import subprocess
from scwf import get_env
from collections import namedtuple

mk_TF_input_sc = f"{os.path.dirname(__file__)}/analysis_TF.r"
mk_cellchat_input_sc = f"{os.path.dirname(__file__)}/analysis_Cellhat.r"
plot_corr_sc = f"{os.path.dirname(__file__)}/plot_corr.r"

def calculate_correlations(df, groupby_cols, col_x, col_y):
    results = []
    compare_detail = []
    # ('O95727', 'ZFP64', 'Classical_monocytes', 'Ctrl')
    grouped = df.groupby(groupby_cols)
    for group, data in grouped:
        if (len(data)) <3 :
            print(group)
            print(data)

        if len(data) > 1: 
            pearson_corr, pearson_p = pearsonr(data[col_x], data[col_y])
            spearman_corr, spearman_p = spearmanr(data[col_x], data[col_y])
            
            results.append({
                **{col: val for col, val in zip(groupby_cols, group)},
                'pearson_corr': pearson_corr,
                'pearson_p': pearson_p,
                'spearman_corr': spearman_corr,
                'spearman_p': spearman_p
            })
        nPro = len(data[col_x])
        nGene = len(data[col_y])
        group = list(group)
        group.append(nPro)
        group.append(nGene)

        compare_detail.append(list(group))

    return pd.DataFrame(results)


def get_target_npx(
        df, 
        group, 
        group_col_name = 'group',
        sample_col_name = 'sample'):
    df = df.loc[(df[group_col_name] == group), ['UniProt',sample_col_name, 'NPX', group_col_name]]
    return df


def get_target_tf(
        df, 
        cell, 
        group, 
        group_col_name = 'group',
        sample_col_name = 'sample'):
    
    f1 = df['cell'] == cell
    f2 = df[group_col_name] == group
    df_out = df[f1 & f2][['gene', sample_col_name, 'exp','cell', group_col_name]]
    return df_out
    

def get_analysis_df(df_npx, df_tf, cell, group, sample_col_name = 'sample', group_col_name = 'group'):
    df1 = get_target_npx(df_npx, group=group)
    df2 = get_target_tf(df_tf,cell =cell, group=group)

    df3 = pd.merge(df1, df2, how='inner', on = [sample_col_name, group_col_name])
    print(df3)
    return(df3)

def mk_input(
    rds,
    outdir,
    method,
    group_col,
    sample_col,
    celltype_col,
    script_path,
    env='CellChat',
    verbose=False,
    **kwargs
):
    if method not in ['TF', 'cellchat']:
        raise KeyError(f"method must in ['TF', 'cellchat']!")
    
    run_env = get_env(name=env, load=True)
    args = [
        f"--rds {rds}",
        f"--group_col {group_col}",
        f"--sample_col {sample_col}",
        f"--celltype_col {celltype_col}"
    ]
    extra_args = [f"--{k} {v}" for k, v in kwargs.items()]
    args.extend(extra_args)
    cmd = f"{run_env} Rscript {script_path} " + " ".join(args)
    if verbose:
        print(cmd)
    subprocess.check_call(cmd, shell=True)
    if not os.path.exists(outdir):
        os.system(f'mkdir -p {outdir}')
    os.system(f'mv ./{method}_input.tsv {outdir} ')
    
    return f"{outdir}/{method}_input.tsv"


def mk_analysis_data(df_sc_path, df_olink_path, group_sort, method, is_drop_group_info=False, is_drop_cell_info=False):
    analysisData = namedtuple('analysisData', ['sc', 'olink', 'all_cell', 'all_group', 'group_sort', 'method'])
    df_sc = pd.read_table(df_sc_path).dropna()
    df_olink = pd.read_table(df_olink_path).dropna()

    if is_drop_group_info not in [False, 'False']:
        df_olink['group'] = 'all'
        df_sc['group'] = 'all'
        group_sort = 'all'
    
    if is_drop_cell_info not in [False, 'False']:
        df_olink['cell'] = 'all'
        df_sc['cell'] = 'all'

    all_cells = df_sc['cell'].unique()
    all_group = df_olink['group'].unique()
    print(all_cells)
    print(all_group)

    return analysisData(df_sc, df_olink, all_cells, all_group, group_sort, method)


def get_corr_data(analysis_data, outdir):
    res = []
    for cell in tqdm.tqdm(analysis_data.all_cell):
        for gr in analysis_data.all_group:
            df3 = get_analysis_df(analysis_data.olink, analysis_data.sc, cell, gr,  sample_col_name = 'sample')
            correlation_results = calculate_correlations(df3, ['UniProt', 'gene', 'cell', 'group'], 'NPX', 'exp')
            res.append(correlation_results)
    df_out = pd.concat(res).dropna()
        
    df_out.to_csv(f'{outdir}/corr_olink_{analysis_data.method}.csv', index= None)
    print(f'finish {analysis_data.method}!')
    return f'{outdir}/corr_olink_{analysis_data.method}.csv'

def plot_corr(
    corr_file,
    method,
    olink_id_file,
    outdir,
    env='_common',
    verbose=False,
    **kwargs
):
    run_env = get_env(name=env, load=True)

    corr_file = os.path.join(outdir, f"corr_olink_{method}.csv")
    if not os.path.exists(corr_file):
        raise FileNotFoundError(f"Correlation file not found: {corr_file}")
    if not os.path.exists(olink_id_file):
        raise FileNotFoundError(f"Olink ID file not found: {olink_id_file}")

    analysis_dir = os.path.join(outdir, f"analysis_{method}")
    os.makedirs(analysis_dir, exist_ok=True)

    args = [
        f"--corr_file {corr_file}",
        f"--olink_id_file {olink_id_file}",
        f"--outdir {analysis_dir}"
    ]

    extra_args = [f"--{k} {v}" for k, v in kwargs.items()]
    args.extend(extra_args)

    cmd = f"{run_env} Rscript {plot_corr_sc} " + " ".join(args)

    if verbose:
        print("Running command:")
        print(cmd)

    subprocess.check_call(cmd, shell=True)

def main():
    parser = argparse.ArgumentParser(description="Process some arguments.")
        
    # 输入文件和输出目录
    parser.add_argument('--npx', type=str, default='', help='')
    parser.add_argument('--rds', type=str, default='', help='')
    parser.add_argument('--tf_rds', type=str, default='None', help='')

    parser.add_argument('--method', type=str, default='', help='')
    parser.add_argument('--target_pathway', type=str, default='None', help='')
    parser.add_argument('--target_TF', type=str, default='None', help='')
    parser.add_argument('--target_olink', type=str, default='None', help='')
    

    parser.add_argument('--group_sort', type=str, default='None', help='')
    parser.add_argument('--gcpal', type=str, default='npg', help='')

    parser.add_argument('--show_gene_name', type=str, default='False', help='')
    parser.add_argument('--outdir', type=str, default='', help='')

    parser.add_argument('--group_col', type=str, default='', help='')
    parser.add_argument('--sample_col', type=str, default='', help='')
    parser.add_argument('--celltype_col', type=str, default='', help='')
    parser.add_argument('--olink_id_file', type=str, default='None', help='')
    parser.add_argument('--drop_group_info', type=str, default='False', help='')

    args = parser.parse_args()

    if not os.path.exists(args.outdir):
        os.system(f'mkdir -p {args.outdir}')

    if args.method == 'cellchat':
        input_path = mk_input(
            rds=args.rds,
            outdir=args.outdir,
            method='cellchat',
            group_col=args.group_col,
            sample_col=args.sample_col,
            celltype_col=args.celltype_col,
            script_path=mk_cellchat_input_sc,
            env='CellChat',
        )
    elif args.method == 'TF':
        input_path = mk_input(
            rds=args.tf_rds,
            outdir=args.outdir,
            method='TF',
            group_col=args.group_col,
            sample_col=args.sample_col,
            celltype_col=args.celltype_col,
            script_path=mk_TF_input_sc,
            env='DoRothEA',
        )

    # Create the analysis data
    analysis_data = mk_analysis_data(
        df_sc_path=input_path,
        df_olink_path=args.npx,
        group_sort=args.group_sort,
        method=args.method,
        is_drop_group_info=args.drop_group_info
    )

    # Get correlation data
    corr_file = get_corr_data(analysis_data, args.outdir)

    # Set up arguments for plot_corr
    plot_args = {
        'target_gene': args.target_pathway if args.method == 'cellchat' else args.target_TF,
        'show_gene_name': args.show_gene_name,
        'group_sort': args.group_sort,
        'gcpal': args.gcpal,
        'target_olink': args.target_olink
    }

    # Plot the correlations
    plot_corr(
        corr_file=corr_file,
        method=args.method,
        olink_id_file=args.olink_id_file,
        outdir=args.outdir,
        env='_common',
        **plot_args
    )

if __name__ == '__main__':
    main()