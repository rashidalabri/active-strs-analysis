import pandas as pd
from pyliftover import LiftOver
from common import apply_filters

estrs = pd.read_csv(snakemake.input[0], sep='\t')
estrs = estrs.rename({
    'chrom': 'chr',
    'str.start': 'start',
    'str.end': 'stop',
    'str.motif.forward': 'motif'
}, axis='columns')

former_cols = ['chr', 'start', 'stop', 'motif']
latter_cols = [col for col in estrs.columns if col not in former_cols]

estrs = estrs[former_cols + latter_cols]

# Only include finely-mapped eSTRs (CAVIAR score >0.3)
estrs = estrs[estrs['score'] > 0.3]

estrs_hg19 = estrs.copy()
estrs_hg38 = estrs.copy()

# Liftover
lo = LiftOver('hg19', 'hg38')
to_drop = []
for i, row in estrs_hg38.iterrows():
    start = lo.convert_coordinate(row['chr'], row['start'])
    stop = lo.convert_coordinate(row['chr'], row['stop'])
    if len(start) > 0 and len(stop) > 0:
        estrs_hg38.loc[i, 'start'] = start[0][1]
        estrs_hg38.loc[i, 'stop'] = stop[0][1]
    else:
        to_drop.append(i)

estrs_hg38 = estrs_hg38.drop(to_drop, axis='rows')

estrs_hg19 = apply_filters(estrs_hg19, snakemake.config)
estrs_hg38 = apply_filters(estrs_hg38, snakemake.config)

estrs_hg19.to_csv(snakemake.output['hg19'], sep='\t', index=False)
estrs_hg38.to_csv(snakemake.output['hg38'], sep='\t', index=False)
