import re
import pandas as pd
import pysam

class LocusAlignment(object):
    extract_nodes_p = re.compile('(\d+)\[((?:\d+[A-Z])+)\]', re.IGNORECASE)
    extract_operations_p = re.compile('(\d+)([A-Z])', re.IGNORECASE)

    def __init__(self) -> None:
        self.alignments = {}
        self.alignments_node_id = {}
    
    def add(self, node_id, position, operation) -> None:
        if position not in self.alignments:
            self.alignments[position] = {op: 0 for op in 'MXID'}
        self.alignments[position][operation] += 1

        # if position in self.alignments_node_id:
            # if self.alignments_node_id[position] != node_id:
            #     print(self.alignments_node_id[position], node_id)
            # assert self.alignments_node_id[position] == node_id
        self.alignments_node_id[position] = node_id
    
    def add_range(self, node_id, position, operation, count) -> None:
        for i in range(count):
            self.add(node_id, position + i, operation)

    def parse_graph_cigar(self, cigar_string, start_position, only_flanks=True):
        nodes = self.extract_nodes_p.findall(cigar_string)
        nodes = [(int(node_id), cigar_string) for node_id, cigar_string in nodes]
        
        if only_flanks:
            nodes = [nodes[i] for i in [0, -1]]

        current_position = start_position
        for node in nodes:
            current_position = self.parse_cigar_node(node, current_position)

    def parse_cigar_node(self, node, start_position) -> int:
        node_id, cigar_string = node
        operations = self.extract_operations_p.findall(cigar_string)
        current_position = start_position
        for count, op in operations:
            count = int(count)
            if op not in 'S':
                self.add_range(node_id, current_position, op, count)
            if op in 'MX':
                current_position += count
        return current_position

    def to_df(self) -> pd.DataFrame:
        df = pd.DataFrame({
            'node_id': self.alignments_node_id.values(),
            'position': self.alignments.keys(),
            'matches': [self.alignments[pos]['M'] for pos in self.alignments],
            'mismatches': [self.alignments[pos]['X'] for pos in self.alignments],
            'insertions': [self.alignments[pos]['I'] for pos in self.alignments],
            'deletions': [self.alignments[pos]['D'] for pos in self.alignments],
        })
        return df


# bam_file_path = snakemake.input[0]
# out_file_path = snakemake.output[0]

locus_alignments = {}
bam = pysam.AlignmentFile('HG002_realigned.sorted.bam')
for read in bam.fetch(until_eof=True):
    locus_id, _, cigar_string = read.get_tag('XG').split(',')
    if locus_id not in locus_alignments:
        locus_alignments[locus_id] = LocusAlignment()
    
    alignment = locus_alignments[locus_id]
    alignment.parse_graph_cigar(cigar_string, read.reference_start)


dfs = []
for locus_id in locus_alignments:
    df = locus_alignments[locus_id].to_df()
    # df['name'] = snakemake.wilcards['sample']
    df['locus_id'] = locus_id
    dfs.append(df)

concated = pd.concat(dfs, ignore_index=True)
# concated.to_csv(out_file_path, sep='\t', index=False)
concated.to_csv('out.tsv', sep='\t', index=False)


