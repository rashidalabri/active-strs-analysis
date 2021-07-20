import pandas as pd
from common import apply_filters

novel = pd.read_csv(snakemake.input[0], sep='\t')

former_cols = ['chr', 'start', 'stop', 'motif']
latter_cols = [col for col in novel.columns if col not in former_cols]

novel = novel[former_cols + latter_cols]

novel = apply_filters(novel, snakemake.config)

novel.to_csv(snakemake.output[0], sep='\t', index=False)