import re
import pandas as pd
import pysam
from collections import defaultdict
from time import time

extract_nodes_p = re.compile('(\d+)\[((?:\d+[A-Z])+)\]', re.IGNORECASE)
extract_operations_p = re.compile('(\d+)([A-Z])', re.IGNORECASE)

class LocusFlankData(object):

    def __init__(self) -> None:
        self.flank_data = {
            'left': {},
            'right': {}
        } 
    
    def add(self, flank, position, operation) -> None:
        if position not in self.flank_data[flank]:
            self.flank_data[flank][position] = {op: 0 for op in 'MXID'}
        self.flank_data[flank][position][operation] += 1
    
    def add_range(self, flank, position, operation, count) -> None:
        for i in range(count):
            self.add(flank, position + i, operation)

    def parse_graph_cigar(self, cigar_string):
        nodes = extract_nodes_p.findall(cigar_string)
        nodes = [(int(node_id), cigar_string) for node_id, cigar_string in nodes]
        
        for node in nodes:
            node_id, cigar_string = node
            if node_id == 0:
                self.parse_cigar_string('left', cigar_string)
            elif node_id == 2:
                self.parse_cigar_string('right', cigar_string)
    
    def parse_cigar_string(self, flank, cigar_string) -> None:
        operations = extract_operations_p.findall(cigar_string)

        if flank == 'left':
            operations.reverse()

        current_position = 1
        for count, op in operations:
            count = int(count)
            if op not in 'S':
                self.add_range(flank, current_position, op, count)
            if op in 'MX':
                current_position += count

    def to_df(self) -> pd.DataFrame:
        dfs = []
        for flank in self.flank_data:
            data = self.flank_data[flank]
            df = pd.DataFrame({
                'position': data.keys(),
                'matches': [data[pos]['M'] for pos in data],
                'mismatches': [data[pos]['X'] for pos in data],
                'insertions': [data[pos]['I'] for pos in data],
                'deletions': [data[pos]['D'] for pos in data],
            })
            # df['flank'] = flank
            if flank == 'left':
                df['position'] *= -1
            dfs.append(df)
        return pd.concat(dfs, ignore_index=True)


bam_file_path = snakemake.input['bam']
out_file_path = snakemake.output[0]
sample = snakemake.wildcards['sample']

# bam_file_path = 'resources/realigned_bam/active/HG00118/HG00118_realigned.bam'
# out_file_path = 'mutations.h5'
# sample = 'HG00118'

# Get all cigar strings
loci_cigars = defaultdict(list)
bam = pysam.AlignmentFile(bam_file_path)
for read in bam.fetch(until_eof=True):
    locus_id, _, cigar_string = read.get_tag('XG').split(',')
    loci_cigars[locus_id].append(cigar_string)

print('Finished fetching {} reads. Parsing...'.format(sample))

# Parse and store, while continously freeing up memory
all_locus_ids = list(loci_cigars.keys())
for locus_id in all_locus_ids:
    t0 = time()
    locus_flank_data = LocusFlankData()
    for cigar_string in loci_cigars[locus_id]:
        locus_flank_data.parse_graph_cigar(cigar_string)
    df = locus_flank_data.to_df()
    df['locus_id'] = locus_id
    df = df.set_index(['locus_id', 'position'])
    df.to_hdf(out_file_path, sample, format='table', append=True, min_itemsize={ 'locus_id' : 32 })
    del loci_cigars[locus_id]
    print('Took {}s to proccess one locus.'.format(time()-t0))
