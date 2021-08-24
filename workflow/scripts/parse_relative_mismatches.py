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
            df['flank'] = flank
            if flank == 'left':
                df['position'] *= -1
            dfs.append(df)
        return pd.concat(dfs, ignore_index=True)


# bam_file_path = snakemake.input['bam']
# out_file_path = snakemake.output[0]

bam_file_path = 'resources/realigned_bam/active/HG00118/HG00118_realigned.bam'
out_file_path = 'mutations.tsv'


loci_flank_data = {}
bam = pysam.AlignmentFile(bam_file_path)
for read in bam.fetch(until_eof=True):
    locus_id, _, cigar_string = read.get_tag('XG').split(',')
    if locus_id not in loci_flank_data:
        loci_flank_data[locus_id] = LocusFlankData()
    locus_flank_data = loci_flank_data[locus_id]
    locus_flank_data.parse_graph_cigar(cigar_string)

print('Done reading')


loci_ids = list(loci_flank_data.keys())

if len(loci_ids) > 0:
    locus_id = loci_ids[0]
    df = loci_flank_data[locus_id].to_df()
    df['sample'] = 'HG00118_realigned.bam'
    # df['sample'] = snakemake.wildcards['sample']
    df['locus_id'] = locus_id
    df.to_csv(out_file_path, sep='\t', index=False)
    del loci_flank_data[locus_id]

for locus_id in loci_ids[1:]:
    df = loci_flank_data[locus_id].to_df()
    df['sample'] = 'HG00118_realigned.bam'
    # df['sample'] = snakemake.wildcards['sample']
    df['locus_id'] = locus_id
    df.to_csv(out_file_path, mode='a', header=False)
    del loci_flank_data[locus_id]