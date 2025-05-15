import scanpy as sc
import numpy as np
import pandas as pd
import os 
import argparse

# scanpy mitochonrial variable name
MITO_VAR = "mito"
NORMALIZED_LAYER = "normalized"

class Scanpy_wrapper:
    def __init__(self, matrix_file, meta, outfile=None, run_normailze = False):
        # data
        self.outfile = outfile
        self.meta = meta
        self.matrix_file = matrix_file
        self.adata = sc.read_10x_mtx(
            self.matrix_file,
            var_names="gene_symbols",
        )
        self.adata.layers["raw"] = self.adata.X.copy()
        self.h5ad_file = self.outfile
        self.run_normalize = run_normailze

        print('load data success!')

    def calculate_qc_metrics(self):

        self.adata.var[MITO_VAR] = self.adata.var_names.str.upper().str.startswith(
                "MT-"
            )

        sc.pp.calculate_qc_metrics(
            self.adata,
            qc_vars=[MITO_VAR],
            percent_top=None,
            use_raw=False,
            log1p=False,
            inplace=True,
        )
    
    
    def normalize(self):
        """
        sc.pp.normalize_per_cell() and sc.pp.log1p()
        """

        sc.pp.normalize_total(
            self.adata,
            target_sum=1e4,
            inplace=True,
        )
        sc.pp.log1p(
            self.adata,
        )
        self.adata.layers[NORMALIZED_LAYER] = self.adata.X

    def write_h5ad(self):
        self.meta_df = pd.read_table(self.meta)
        self.meta_df = self.meta_df.set_index('barcode', drop=False)

        # 只保留在adata中的barcode
        self.meta_df = self.meta_df.loc[self.meta_df.index.isin(self.adata.obs_names)]

        print(self.adata.X.data)
        print(self.adata.var_names)
        self.adata.obs = self.adata.obs.join(self.meta_df)
        print(self.adata.obs)
        self.adata.write(self.h5ad_file)


    def run(self):
        self.calculate_qc_metrics()
        if self.run_normalize in [True, 'True', 'TRUE']:
            self.normalize()
        self.write_h5ad()
    

def main():
    parser = argparse.ArgumentParser(description='Scanpy wrapper command line tool')
    parser.add_argument('--matrix_file', type=str, help='Path to the 10x matrix file')
    parser.add_argument('--meta', type=str, help='Path to the metadata file')
    parser.add_argument('--outfile', type=str, help='Output directory')
    parser.add_argument('--is_normailze', type=str, help='', default='False')

    args = parser.parse_args()

    wrapper = Scanpy_wrapper(args.matrix_file, args.meta, args.outfile, args.is_normailze)
    wrapper.run()

if __name__ == '__main__':
    main()