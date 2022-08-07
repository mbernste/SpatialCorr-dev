import pandas as pd
import numpy as np

df = pd.read_csv('top_genes.tsv', sep='\t', index_col=0)

pairs = []
while len(pairs) < 20:
    g1 = np.random.choice(df.index, size=1)[0]
    g2 = np.random.choice(df.index, size=1)[0]
    if g1 != g2 and frozenset([g1, g2]) not in [frozenset([x[0], x[1]]) for x in pairs]:
        pairs.append((g1, g2))

df_out = pd.DataFrame(
    data=pairs,
    columns=['gene_1', 'gene_2']
)
df_out.to_csv('random_pairs_1.tsv', sep='\t')

