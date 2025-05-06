import pandas as pd


'''
脑膜瘤拷贝数变异：
1：常见的染色体臂丢失：22q、14q、1p、18q、6q、10q、11p、7p 和 4p，肿瘤级别越高，染色体丢失出现频率越大；
2：NF2突变，22q12.2，一半脑膜瘤中有，1-3级中都有；
3：TRAF7突变，16p13，五分之一脑膜瘤中有，一般是1级；
4：AKT1突变，c.49G>A；一般是1级；
5：KLF4突变，c.1225A>C；一般是1级；
6：SMO突变，一般是1级；
7：PIK3CA突变，3q26.3，一般是1级；
8：SUFU突变，TERT突变，通常是2级或3级
'''

def parse_name_input(s):
    if 'p' in s:
        s_attr = s.split('p')
        chrom = f"chr{s_attr[0]}"
        try:
            chrname = f"p{s_attr[2]}"
        except IndexError:
            chrname = 'p'
    
    if 'q' in s:
        s_attr = s.split('q')
        chrom = f"chr{s_attr[0]}"
        try:
            chrname = f"q{s_attr[1]}"
        except IndexError:
            chrname = 'q'
    
    return chrom, chrname


def extract(df, chrom, chrname):
    df_ext = df[(df['#chrom'] == chrom) & (df['name'].str.contains(chrname))]
    return df_ext

def load_cytoband(file):
    df = pd.read_table(file).dropna()
    return df

df = load_cytoband("/public/home/zhliu/singlecell/sc_PIPLINE2/inferCNV/data/cytoband/hg38/cytoband.txt")

chrom, chrname = parse_name_input('22q')
target_df = extract(df, chrom, chrname)
target_df.to_csv("target.bed", index = None, sep = '\t', header=False)