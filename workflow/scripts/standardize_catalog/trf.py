import pandas as pd
from common import apply_filters

trf = pd.read_csv(snakemake.input[0], sep='\t')
trf.columns = ['bin', 'chr', 'start', 'stop', 'name', 'period', 'copyNum', 'consensusSize', 'perMatch', 'perIndel', 'score', 'A', 'C', 'G', 'T', 'entropy', 'motif']

former_cols = ['chr', 'start', 'stop', 'motif']
latter_cols = [col for col in trf.columns if col not in former_cols]

trf = trf[former_cols + latter_cols]

trf = apply_filters(trf, snakemake.config)

trf.to_csv(snakemake.output[0], sep='\t', index=False)