import pandas as pd

class LooseDataFrame():
    
    def __init__(self, columns):
        self.columns = columns
        self.table = {column: [] for column in columns}
    
    def append(self, row):
        for column in row:
            self.table[column].append(row[column])
    
    def to_df(self):
        return pd.DataFrame(self.table)