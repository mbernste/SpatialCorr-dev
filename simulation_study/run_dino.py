import subprocess
import sys
import os
from os.path import join
from optparse import OptionParser

# This is pretty dirty. Need to just write a single
# script for any number of genes...
N_GENES_TO_CMD = {
    'standard_pair': './dino_scripts/DinoTestScript.R',
    2: './dino_scripts/Dino_two_genes.R',
    3: './dino_scripts/Dino_three_genes.R',
    4: './dino_scripts/Dino_four_genes.R',
    5: './dino_scripts/Dino_five_genes.R',
    8: './dino_scripts/Dino_eight_genes.R',
    9: './dino_scripts/Dino_nine_genes.R',
}

def main():
    parser = OptionParser()
    parser.add_option("-g", "--num_genes", help='Number of genes')
    (options, args) = parser.parse_args()

    in_dir = args[0]
    out_dir = args[1]
    if options.num_genes is not None:
        n_genes = int(options.num_genes)
    else:
        n_genes = 'standard_pair'

    for elem in sorted(os.listdir(in_dir)):
        in_fname = join(in_dir, elem)
        out_fname = join(out_dir, elem)
        cmd = f'Rscript {N_GENES_TO_CMD[n_genes]} {in_fname} {out_fname}'
        print(f"Running: {cmd}")
        subprocess.run(cmd, shell=True)

if __name__ == '__main__':
    main()
