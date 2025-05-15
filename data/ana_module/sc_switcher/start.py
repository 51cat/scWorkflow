import subprocess 
import os
import argparse
from scwf import get_env



class Transfer:
    def __init__(
            self,
            input,
            outdir,
            split_col,
            keep,
            keep_str,
            pyexec = '/usr/bin/python',
            fr = 'seurat',
            to = 'h5ad',
            
            ):
        self.input = input
        self.outdir = outdir

        self.split_col = split_col
        self.keep = keep
        self.keep_str = keep_str
        self.pyexec = pyexec
        self.pyexec = pyexec

        self.fr = fr
        self.to = to
        self.env = get_env("sc_switcher", load = True)
    
    def tarsfer(self):
        sc_path = f'{os.path.dirname(__file__)}/transfer.r'

        cmd = (
            f"{self.env} Rscript {sc_path} "
            f"--rds {self.input} "
            f"--outdir {self.outdir} "
            f"--pyexec {self.pyexec} "
            f"--split_col {self.split_col} "
            f"--keep {self.keep} "
            f"--keep_str {self.keep_str} "
            f"--from {self.fr} "
            f"--to {self.to} "
        )
        subprocess.check_call(cmd, shell=True)


def main():
    parser = argparse.ArgumentParser(
        description=""
    )
    parser.add_argument("--input", required=True, help="Input Seurat RDS file")
    parser.add_argument("--outdir", required=True, help="Output directory for .h5ad files")
    parser.add_argument("--pyexec", default="/usr/bin/python", help="Path to Python executable for reticulate")
    parser.add_argument("--split_col", default='None', help="Column in metadata to split by")
    parser.add_argument("--keep", default='None', help="Comma-separated values to keep from split_col")
    parser.add_argument("--keep_str", default='None', help='Filter expression, e.g., "nFeature_RNA > 200 & percent.mt < 5"')
    parser.add_argument("--fr", default='seurat', help='')
    parser.add_argument("--to", default='h5ad', help='')

    args = parser.parse_args()


    converter = Transfer(
        input=args.input,
        outdir=args.outdir,
        pyexec=args.pyexec,
        split_col=args.split_col,
        keep=args.keep,
        keep_str=args.keep_str,
        fr = args.fr,
        to = args.to,
    )

    converter.tarsfer()

if __name__ == '__main__':
    main()