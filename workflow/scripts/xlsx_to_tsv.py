import pandas as pd

xls = pd.ExcelFile(snakemake.input[0])
df = pd.read_excel(xls, snakemake.params['sheet'])
df.to_csv(snakemake.output[0], sep='\t', index=False)
