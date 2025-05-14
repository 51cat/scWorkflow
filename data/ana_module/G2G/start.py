import anndata
import numpy as np
import seaborn as sb
import numpy as np
import pandas as pd
import warnings
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('agg')
import blitzgsea

from tqdm import tqdm
warnings.filterwarnings("ignore")
import os
import sys
from genes2genes import Main
from genes2genes import ClusterUtils
from genes2genes import TimeSeriesPreprocessor
from genes2genes import PathwayAnalyser
from genes2genes import VisualUtils

from .find import Find_Cluster_thr

from optbinning import ContinuousOptimalBinning
import argparse

import pickle
from .add_log import add_log

DEFAULT_GENESET = [
                'WikiPathway_2021_Human', 
                'KEGG_2021_Human',
                'Reactome_Pathways_2024',
                'GO_Biological_Process_2021',
                'GO_Molecular_Function_2023',
                'GO_Cellular_Component_2017b'
                ]


def find_steady_points(array, output):

    # 计算相邻元素的差异（增量）
    data = np.cumsum(array)
    diff = np.abs(np.diff(data))
    
    threshold = np.median(diff)
    # 找到增量小于阈值的点
    steady_indices = np.where(diff < threshold)[0] + 1  # 索引从1开始
    
    # 如果找到了平稳点，取第一个平稳点的索引
    if len(steady_indices) > 0:
        first_steady_index = steady_indices[0]
    else:
        first_steady_index = None

    # 绘图
    #plt.plot(data, label='Cumulative Data')
    #plt.scatter(steady_indices, data[steady_indices], color='red', label='Steady Points')
    
    #if first_steady_index is not None:
    #    plt.axvline(x=first_steady_index, color='green', linestyle='--', label='First Steady Point')
    
    #plt.xlabel('Index')
    #plt.ylabel('Value')
    #plt.title('Steady Points in Cumulative Data')
    #plt.legend()
    #plt.savefig(f"{output}", dpi=300, bbox_inches='tight')
    #plt.close()
    
    return steady_indices[0]


def get_pathway_alignment_stat(aligner, GENE_LIST, pathway_name, cluster=False, FIGSIZE = (14,7)):
    
    print('Gene set: ======= ', pathway_name)
    perct_A = []
    perct_S = []
    perct_T = []
    for gene in GENE_LIST:
        series_match_percent = aligner.results_map[gene].get_series_match_percentage()
        perct_A.append(series_match_percent[0])
        perct_S.append(series_match_percent[1])
        perct_T.append(series_match_percent[2])

    print('mean matched percentage: ', round(np.mean(perct_A),2),'%' )
    #print('mean matched percentage wrt ref: ',round(np.mean(perct_S),2),'%'  )
    #print('mean matched percentage wrt query: ', round(np.mean(perct_T),2),'%' )
    average_alignment, alignment_path =  ClusterUtils.get_cluster_average_alignments(aligner, GENE_LIST)
    mat = ClusterUtils.get_pairwise_match_count_mat(aligner,GENE_LIST )
    print('Average Alignment: ', VisualUtils.color_al_str(average_alignment), '(cell-level)')
    print('- Plotting average alignment path')

    return alignment_path, aligner,  mat



def redirect_print_output(func, output_file, *args, **kwargs):
    original_stdout = sys.stdout
    try:
        with open(output_file, 'w') as f:
            sys.stdout = f
            func(*args, **kwargs)
    finally:
        sys.stdout = original_stdout

def extract_genes(adata1, adata2, seed=None):
    if seed is not None:
        np.random.seed(seed)
    gene_indices = np.random.choice(adata1.shape[1], size=1500, replace=False)
    sub_adata1 = adata1[:, gene_indices].copy()
    sub_adata2 = adata2[:, gene_indices].copy()
    return sub_adata1, sub_adata2

def merge_on_index_intersection(df1, df2):
    common_index = df1.index.intersection(df2.index)
    df1_filtered = df1.loc[common_index]
    df2_filtered = df2.loc[common_index]
    merged_df = pd.concat([df1_filtered, df2_filtered], axis=1)
    return merged_df

def mk_manual_color(fearture_lst):
    fearture_list_unique = list(set(fearture_lst))
    col = np.array(sb.color_palette('colorblind', len(fearture_list_unique)))[range(len(fearture_list_unique))]
    joint_cmap={}
    for inx, f in enumerate(fearture_list_unique):
        joint_cmap.update({f: col[inx]})
    return joint_cmap


# 构建G2G分析数据
class prepG2Ginput:
    def __init__(self, ref_h5ad, query_h5ad, ann_col, is_downsample =False):
        self.ref_h5ad = ref_h5ad
        self.query_h5ad = query_h5ad
        self.ann_col = ann_col
        self.is_downsample = is_downsample
    @add_log
    def load_h5ad(self):
        self.ref = anndata.read_h5ad(self.ref_h5ad)
        self.query = anndata.read_h5ad(self.query_h5ad)
        if self.is_downsample == "True":
            print("downsample! - 1500 * genes")
            self.ref,  self.query = extract_genes(self.ref, self.query)

    def set_psetime(self, pse_df_path):
        self.pse_df = pd.read_csv(pse_df_path, index_col=0)
    
    def add_pse_time(self):
        self.ref.obs = merge_on_index_intersection(self.ref.obs, self.pse_df)
        self.query.obs = merge_on_index_intersection(self.query.obs, self.pse_df)

        # 将inf值替换为NaN
        self.ref.obs.replace([np.inf, -np.inf], np.nan, inplace=True)
        self.query.obs.replace([np.inf, -np.inf], np.nan, inplace=True)

        # 删除包含NaN的行
        self.ref.obs.dropna(inplace=True)
        self.query.obs.dropna(inplace=True)

        self.ref = self.ref[list(self.ref.obs.index), :].copy()
        self.query = self.query[list(self.query.obs.index), :].copy()

        print(self.ref)
        print(self.query)

    def scale_psetime(self):
        if "time" not in self.ref.obs.columns:
            raise KeyError(f"time column not in ref {self.ref.obs.columns}")

        if "time" not in self.ref.obs.columns:
            raise KeyError(f"time column not in query {self.query.obs.columns}")

        self.ref.obs['time'] = TimeSeriesPreprocessor.Utils.minmax_normalise(np.asarray(self.ref.obs['time']))
        self.query.obs['time'] = TimeSeriesPreprocessor.Utils.minmax_normalise(np.asarray(self.query.obs['time']))
        print(min(self.ref.obs['time']), max(self.ref.obs['time']))
        print(min(self.query.obs['time']), max(self.query.obs['time']))

    @add_log
    def get_input_data(self):
        return self.ref, self.query


# G2G
class G2GAnalysiser:
    def __init__(self, 
            ref, 
            query, 
            ref_name, 
            query_name,
            ann_col, 
            outdir,
            target_gene = None,
            find_cluster_method = 'distance,auto'):
        self.ref = ref
        self.query = query
        self.ref_name = ref_name
        self.query_name = query_name
        self.ann_col = ann_col
        self.outdir = outdir 
        self.target_gene = target_gene
        self.find_ncluster_method_use, self.find_ncluster_method_th = find_cluster_method.split(",")
        
        if not os.path.exists(outdir):
            os.system(f'mkdir -p {outdir}')
    
    def mk_step_outdir(self, name):
        os.system(f"mkdir -p {self.outdir}/{name}")
        return f"{self.outdir}/{name}"
    
    @add_log
    def analysis_pseudotime_distributions(self, nbins):
        # 分析ref和query之间的伪时间的分布差异
        # 创建输出目录
        outdir_step = self.mk_step_outdir('01.pseudotime_distributions')

        # 1. # Visualize the pseudotime distributions
        sb.kdeplot(self.ref.obs['time'], fill=True, label=self.ref_name, color='forestgreen') 
        sb.kdeplot(self.query.obs['time'], fill=True, label=self.query_name, color='midnightblue'); 
        plt.xlabel('pseudotime')
        plt.legend(ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.3), frameon=False)
        plt.tight_layout()
        plt.savefig(f'{outdir_step}/pseudotime_distributions.png', dpi=300, bbox_inches='tight')
        plt.savefig(f'{outdir_step}/pseudotime_distributions.pdf', dpi=300, bbox_inches='tight')

        x = np.asarray(self.ref.obs.time)
        optb = ContinuousOptimalBinning(name='pseudotime', dtype="numerical")
        optb.fit(x, x)
        x = np.asarray(self.query.obs.time)
        optb = ContinuousOptimalBinning(name='pseudotime', dtype="numerical")
        optb.fit(x, x)
        # define annotation column name in the adata obs
        annotation_colname = self.ann_col
        self.annotation_colname = self.ann_col

        self.ref.obs[annotation_colname] = [x for x in self.ref.obs[self.ann_col]] 
        self.query.obs[annotation_colname] = [x for x in self.query.obs[self.ann_col]] 
        self.joint_cmap = mk_manual_color(list(self.ref.obs[annotation_colname]) + list(self.query.obs[annotation_colname]))

        print(self.joint_cmap)

        VisualUtils.plot_pseudotime_dists_with_interpolation_points(self.ref, self.query, nbins)
        plt.savefig(f"{outdir_step}/pseudotime_dists_with_interpolation_points.png", dpi=300, bbox_inches='tight')
        plt.savefig(f"{outdir_step}/pseudotime_dists_with_interpolation_points.pdf", dpi=300, bbox_inches='tight')
        plt.close()

        VisualUtils.plot_celltype_barplot(self.ref, nbins, annotation_colname, self.joint_cmap, legend=True)
        plt.title(self.ref_name)
        plt.tight_layout()
        plt.savefig(f"{outdir_step}/{self.ref_name}_celltype_barplot.png", dpi=300, bbox_inches='tight')
        plt.savefig(f"{outdir_step}/{self.ref_name}_celltype_barplot.pdf", dpi=300, bbox_inches='tight')
        plt.close()

        VisualUtils.plot_celltype_barplot(self.query, nbins, annotation_colname, self.joint_cmap, legend=True)
        plt.title(self.query_name)
        plt.tight_layout()
        plt.savefig(f"{outdir_step}/{self.query_name}_celltype_barplot.png", dpi=300, bbox_inches='tight')
        plt.savefig(f"{outdir_step}/{self.query_name}_celltype_barplot.pdf", dpi=300, bbox_inches='tight')
        plt.close()
    
    @add_log
    def run_G2G_alignment(self, nbins, threshold_similarity = 0.3):
        # 进行伪时间水平的基因比对分析
        outdir_step = self.mk_step_outdir('02.G2G_alignment_stat')
        gene_list = self.ref.var_names 
        print(len(gene_list),'genes')

        aligner = Main.RefQueryAligner(self.ref, self.query, gene_list, nbins)
        # 执行比对, 最慢的一步
        aligner.align_all_pairs() 

        ###Aggregate (average) cell-level alignment across all aligned genes
        # 整体比对结果展示
        aligner.get_aggregate_alignment()
        plt.savefig(f"{outdir_step}/alignment_heatmap.png", dpi=300, bbox_inches='tight')
        plt.savefig(f"{outdir_step}/alignment_heatmap.pdf", dpi=300, bbox_inches='tight')
        plt.close()
        # 输出每个基因的比对结果
        #3. Analysing gene-level alignments
        #Ranking genes based on their alignment similarities
        df = aligner.get_stat_df() # ordered genes according to alignment similarity statistics 
        df.to_csv(f"{outdir_step}/gene_stat.csv")
        plt.savefig(f"{outdir_step}/gene_stat.png",dpi=300, bbox_inches='tight')
        plt.savefig(f"{outdir_step}/gene_stat.pdf",dpi=300, bbox_inches='tight')
        plt.close()
        
        # 输出logfc展示
        VisualUtils.plot_alignmentSim_vs_l2fc(df)
        plt.savefig(f"{outdir_step}/plot_alignmentSim_vs_l2fc.png",dpi=500, bbox_inches='tight')
        plt.savefig(f"{outdir_step}/plot_alignmentSim_vs_l2fc.pdf",dpi=500, bbox_inches='tight')
        plt.close()

        self.aligner = aligner
        
        # 展示在伪时间上差异最大的gene, 并做GSEA分析
        # Gene-set overrepresentation analysis on the top dissimilar genes
        topDEgenes = df[df['alignment_similarity_percentage'] <=threshold_similarity]
        print(topDEgenes)

        if self.target_gene not in ['None', None]:
            
            # 绘制指定基因
            new_target_gene = []
            all_gene = set(list(aligner.results_map.keys()))

            for g in self.target_gene:
                if g not in all_gene:
                    self.run_G2G_alignment.logger.warning(f"{g} not found in data!")
                else:
                    new_target_gene.append(g)

            print('plot target gene!')
            os.system(f"mkdir -p {outdir_step}/target_gene_show/")
            outdir_step_gene_target = f"{outdir_step}/target_gene_show/"
            
            for show in tqdm(list(set(new_target_gene))):
                gene_obj = aligner.results_map[show]
                # Alignment landscape of costs (Note: dashed black path is the optimal alignment)
                gene_obj.landscape_obj.plot_alignment_landscape()
                plt.savefig(f"{outdir_step_gene_target}/{show}_heatmap.png", dpi=300, bbox_inches='tight')
                plt.savefig(f"{outdir_step_gene_target}/{show}_heatmap.pdf", dpi=300, bbox_inches='tight')
                plt.close()

                VisualUtils.show_gene_alignment(show, aligner, self.ref, self.query, self.annotation_colname, self.joint_cmap)
                plt.savefig(f"{outdir_step_gene_target}/{show}_trace_compare.png", dpi=300, bbox_inches='tight')
                plt.savefig(f"{outdir_step_gene_target}/{show}_trace_compare.pdf", dpi=300, bbox_inches='tight')
                plt.close()

                VisualUtils.visualize_gene_alignment(aligner.results_map[show], self.ref, self.query, self.annotation_colname, cmap=self.joint_cmap)
                plt.savefig(f"{outdir_step_gene_target}/{show}_compare_bar.png", dpi=300, bbox_inches='tight')
                plt.savefig(f"{outdir_step_gene_target}/{show}_compare_bar.pdf", dpi=300, bbox_inches='tight')
                plt.close()
            plt.close()



        
        if topDEgenes.empty or isinstance(topDEgenes, str):
            print("empty topDEgenes!!!!!!!")
        else:
            # 展示top基因
            os.system(f"mkdir -p {outdir_step}/diff_gene_show/")
            outdir_step_gene = f"{outdir_step}/diff_gene_show/"
            
            show_all = topDEgenes['Gene'].to_list()
            
            for show in tqdm(list(set(show_all))):
                gene_obj = aligner.results_map[show]
                # Alignment landscape of costs (Note: dashed black path is the optimal alignment)
                gene_obj.landscape_obj.plot_alignment_landscape()
                plt.savefig(f"{outdir_step_gene}/{show}_heatmap.png", dpi=300, bbox_inches='tight')
                plt.savefig(f"{outdir_step_gene}/{show}_heatmap.pdf", dpi=300, bbox_inches='tight')
                plt.close()

                VisualUtils.show_gene_alignment(show, aligner, self.ref, self.query, self.annotation_colname, self.joint_cmap)
                plt.savefig(f"{outdir_step_gene}/{show}_trace_compare.png", dpi=300, bbox_inches='tight')
                plt.savefig(f"{outdir_step_gene}/{show}_trace_compare.pdf", dpi=300, bbox_inches='tight')
                plt.close()

                VisualUtils.visualize_gene_alignment(aligner.results_map[show], self.ref, self.query, self.annotation_colname, cmap=self.joint_cmap)
                plt.savefig(f"{outdir_step_gene}/{show}_compare_bar.png", dpi=300, bbox_inches='tight')
                plt.savefig(f"{outdir_step_gene}/{show}_compare_bar.pdf", dpi=300, bbox_inches='tight')
                plt.close()
            plt.close()

            # gsea
            # Calling wrapper function for GSEAPy enrichr inferface
            pathway_df = PathwayAnalyser.run_overrepresentation_analysis(list(topDEgenes), TARGET_GENESETS=DEFAULT_GENESET) 
            print(pathway_df)
            pathway_df.to_csv(f"{outdir_step}/top_gene_GSEA.csv")
            

    @add_log
    def clust_alignment(self):
        outdir_step = self.mk_step_outdir('03.cluster_alignment')
        # Clustering alignments
        # important: Running experiment to determine the distance threshold for alignment clusters from hierarchical clustering. 
        df = ClusterUtils.run_clustering(self.aligner, metric='levenshtein', experiment_mode=True) 
        df.to_csv(f"{outdir_step}/cluster_results_experiment.csv")
        plt.savefig(f"{outdir_step}/cluster_stat_experiment.png",dpi=500, bbox_inches='tight')
        plt.savefig(f"{outdir_step}/cluster_stat_experiment.pdf",dpi=500, bbox_inches='tight')
        plt.close()
        # find best distance 
        # We aim to select a locally optimal threshold that gives a good trade-off between high mean Silhouette score and low number of clusters which can be biologically meaningful.
        #inx = find_steady_points(df["Number of clusters"], f"{outdir_step}/steady_points.png")
        #best_th = df.iloc[inx]['Distance threshold']

        #if self.cluster_th == -1:
        #    pass
        #else:
        #    best_th = self.cluster_th
        
        finder = Find_Cluster_thr(df)
        # self.find_ncluster_method_use, self.find_ncluster_method_th
        if self.find_ncluster_method_use == 'cluster':
            best_th = finder.find_by_cluster(self.find_ncluster_method_th)
        if self.find_ncluster_method_use == 'distance':
            best_th = finder.find_by_distance(self.find_ncluster_method_th)

        print(f"cluster th:  {best_th}")

        # cluster
        ClusterUtils.run_clustering(self.aligner, metric='levenshtein', DIST_THRESHOLD=best_th)
        ClusterUtils.visualise_clusters(self.aligner,n_cols = 4, figsize= (10,6))
        plt.savefig(f"{outdir_step}/cluster_result.png",dpi=500, bbox_inches='tight')
        plt.savefig(f"{outdir_step}/cluster_result.pdf",dpi=500, bbox_inches='tight')
        plt.close()
        
        #VisualUtils.plot_distmap_with_clusters(self.aligner)
        #plt.savefig(f"{outdir_step}/distmap_with_clusters.png",dpi=300, bbox_inches='tight')
        #plt.savefig(f"{outdir_step}/distmap_with_clusters.pdf",dpi=300, bbox_inches='tight')
        #plt.close()

        # 输出每个cluster的基因
        rows = []
        for key, values in self.aligner.gene_clusters.items():
            for value in values:
                rows.append({'Category': key, 'Gene': value})
        
        # 创建DataFrame
        df = pd.DataFrame(rows)
        df.to_csv(f"{outdir_step}/cluster_genes.csv")

        for key, values in self.aligner.gene_clusters.items():
            self.aligner.get_aggregate_alignment_for_subset(list(values))
            plt.savefig(f"{outdir_step}/cluster_{key}_align.png",dpi=300, bbox_inches='tight')
            plt.savefig(f"{outdir_step}/cluster_{key}_align.pdf",dpi=300, bbox_inches='tight')
            plt.close()
        plt.close()
    



def main():
    parser = argparse.ArgumentParser(description="Process some arguments.")
        
    # 输入文件和输出目录
    parser.add_argument('--compare_file', type=str, default='None', help='compare file')
    parser.add_argument('--ref_h5ad', type=str, default='control.h5ad', help='Reference h5ad file')
    parser.add_argument('--query_h5ad', type=str, default='case.h5ad', help='Query h5ad file')
    parser.add_argument('--target_gene_file', type=str, default='None', help='')
        
    parser.add_argument('--ref_name', type=str, default="control", help='Name for the reference')
    parser.add_argument('--query_name', type=str, default="case", help='Name for the query')
        
    parser.add_argument('--ann_col', type=str, default="", help='Annotation column name')
    parser.add_argument('--find_cluster_method', type=str, default='distance,auto', help='')
        
    parser.add_argument('--pse_df', type=str, default="", help='Path to pse_df')
        
    parser.add_argument('--outdir', type=str, default="./test_out2/", help='Output directory')
        
    parser.add_argument('--n_bins', type=int, default=14, help='Number of bins')
    parser.add_argument('--is_scale_pse', default='True', help='Whether to scale pse')
    parser.add_argument('--is_downsample', default='False', help='Whether to downsample the data')
    #parser.add_argument('--threshold_similarity', type=float, default=0.3, help='Threshold for similarity')
        
    args = parser.parse_args()


    #### G2G analysis #### 

    ## 
    # 前处理
    compare_file = args.compare_file
    ref_h5ad = args.ref_h5ad
    query_h5ad = args.query_h5ad
    find_cluster_method =  args.find_cluster_method

    ref_name = args.ref_name
    query_name = args.query_name

    ann_col = args.ann_col
    pse_df = args.pse_df
    outdir = args.outdir
    n_bins = args.n_bins
    is_scale_pse = args.is_scale_pse
    target_gene_file = args.target_gene_file
    
    compare_list= []

    if compare_file in [None, 'None']:
        compare_list = [(ref_name, ref_h5ad, query_name, query_h5ad)]
    else:
        with open(compare_file) as fd:
            for line in fd.readlines():
                line = line.strip('\n')
                ref_name, ref_h5ad, query_name, query_h5ad = line.split('\t')
                compare_list.append((ref_name, ref_h5ad, query_name, query_h5ad))

    
    for ref_name, ref_h5ad, query_name, query_h5ad in compare_list:
        print(ref_name, ref_h5ad, query_name, query_h5ad)
        preper = prepG2Ginput(
            ref_h5ad, query_h5ad, ann_col, is_downsample = args.is_downsample
        )
        preper.load_h5ad()

        if pse_df != '':
            preper.set_psetime(pse_df)
            preper.add_pse_time()
            preper.scale_psetime()

        if is_scale_pse == 'True' :
            preper.scale_psetime()

        ref, query = preper.get_input_data()

        if target_gene_file not in ['None', None]:
            target_gene = []
            with open(target_gene_file) as fd:
                for line in fd.readlines():
                    target_gene.append(line.strip('\n'))
        else:
            target_gene = None

        print(f"target gene show: {target_gene}")

        runner = G2GAnalysiser(
                ref, 
                query, 
                ref_name, 
                query_name,
                ann_col, 
                f"{outdir}/ref_{ref_name}_vs_query_{query_name}/", 
                find_cluster_method = find_cluster_method,
                target_gene = target_gene)

        runner.analysis_pseudotime_distributions(nbins=n_bins)
        runner.run_G2G_alignment(n_bins)
        runner.clust_alignment()

if __name__ == '__main__':
    main()