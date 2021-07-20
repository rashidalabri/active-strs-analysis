import pandas as pd
from common import apply_filters

gangstr = pd.read_csv(snakemake.input[0], sep='\t')
gangstr = gangstr.rename({
    'CHROM': 'chr',
    'POS': 'start',
    'END': 'stop',
    'RU': 'motif'
}, axis='columns')

former_cols = ['chr', 'start', 'stop', 'motif']
latter_cols = [col for col in gangstr.columns if col not in former_cols]

gangstr = gangstr[former_cols + latter_cols]

# Remove loci that do not pass
gangstr = gangstr[gangstr['NA12892_FILTER'] == 'PASS']

# Zero-index start positions
gangstr['start'] = gangstr['start'] - 1

gangstr = apply_filters(gangstr, snakemake.config)

gangstr.to_csv(snakemake.output[0], sep='\t', index=False)