import json
import pandas as pd
import re
from utils import LooseDataFrame

catalog = json.load(open(snakemake.input[0], 'r'))
df = pd.DataFrame(catalog)

old_cols = list(df.columns)
new_cols = ['chr', 'start', 'stop']

# Expand locus IDs into separate chr, start, stop columns
df[new_cols] = df['LocusId'].str.split('_', expand=True)

df['motif'] = ''
new_cols.append('motif')

# Place chr, start, stop, motif columns at the beginning
df = df[new_cols + old_cols]

# Expand variants
begin = re.escape('(')
end = re.escape(')*')
p = re.compile('{}[GCAT]+{}'.format(begin, end))

new_df = LooseDataFrame(df.columns)
for i, row in df.iterrows():
    structure_motifs = p.findall(row['LocusStructure'])
    
    if len(structure_motifs) > 1:
        ref_regions = row['ReferenceRegion']
        variant_ids = row['VariantId']
        variant_types = row['VariantType']
    else:
        ref_regions = [row['ReferenceRegion']]
        variant_ids = [row['VariantId']]
        variant_types = [row['VariantType']]
    
    assert len(structure_motifs) == len(ref_regions) == len(variant_ids) == len(variant_types)

    for i in range(len(structure_motifs)):
        chr, pos = ref_regions[i].split(':')
        start, stop = pos.split('-')

        new_row = row.to_dict()
        new_row['start'] = start
        new_row['stop'] = stop
        # Slice motif to remove regex and keep repeating unit motif
        new_row['motif'] = structure_motifs[i][1:-2]
        new_row['VariantId'] = variant_ids[i]
        new_row['VariantType'] = variant_types[i]
        new_df.append(new_row)

new_df = new_df.to_df()

new_df.to_csv(snakemake.output[0], sep='\t', index=False)
