import re
import pandas as pd
import pysam

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

loci_flank_data = {}
bam = pysam.AlignmentFile(bam_file_path)
for read in bam.fetch(until_eof=True):
    locus_id, _, cigar_string = read.get_tag('XG').split(',')
    if locus_id not in loci_flank_data:
        loci_flank_data[locus_id] = LocusFlankData()
    locus_flank_data = loci_flank_data[locus_id]
    locus_flank_data.parse_graph_cigar(cigar_string)

print('Finished processing {} reads. Storing...'.format(sample))

for locus_id in loci_flank_data:
    df = loci_flank_data[locus_id].to_df()
    df['locus_id'] = locus_id
    df = df.set_index(['locus_id', 'position'])
    df.to_hdf(out_file_path, sample, format='table', append=True)
