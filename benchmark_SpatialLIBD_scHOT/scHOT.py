import os
from os.path import join
import subprocess
import sys
import pandas as pd

def run(df, genes, out_pref, tmp_dir='./scHOT_tmp', script_dir='.', verbose=False):
    
    df_meta = df[['imagerow', 'imagecol']]
    df_expr = df[genes]

    meta_f = join(tmp_dir, '{}.meta.tsv'.format(out_pref))
    expr_f = join(tmp_dir, '{}.expr.tsv'.format(out_pref))
    out_f = join(tmp_dir, '{}.out.tsv'.format(out_pref))
    df_meta = df_meta.reset_index().rename(columns={'index': 'barcode'})
    df_meta.to_csv(meta_f, sep='\t', index=None)
    df_expr = df_expr.reset_index().rename(columns={'index': 'barcode'})
    df_expr.to_csv(expr_f, sep='\t', index=None)

    cmd = 'Rscript {}/scHOT.R {} {} {}'.format(
        script_dir, expr_f, meta_f, out_f        
    )
    if verbose:
        print("Running: ", cmd)
    out = subprocess.run(
        cmd, 
        shell=True, 
        universal_newlines=True, 
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE
    )
    if out.returncode == 1:
        print(out.stdout)
        print(out.stderr)
        return
    else:
        df_res = pd.read_csv(out_f, sep='\t', index_col=0)
    return df_res
    

