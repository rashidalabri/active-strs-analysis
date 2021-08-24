import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import pandas as pd
import seaborn as sns
# import scipy.stats as stats

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

sns.set_context('talk', font_scale=1)

dfs = []
for bed in snakemake.input:
    df = pd.read_table(bed, header=None)
    last_col = list(df.columns)[-1]
    df = df[[last_col]]
    df.columns = ['distance']
    df['catalog'] = bed.split('/')[0].split('.')[0]
    dfs.append(df)
    
data = pd.concat(dfs, ignore_index=True)

pdf = PdfPages(snakemake.output[0])

p = snakemake.params

fig, ax = plt.subplots(figsize=(8,6))
sns.kdeplot(x='distance', data=data[data['distance'].abs() < 20000], ax=ax, hue='catalog', bw_adjust=0.75)
sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
ax.set_title('Distance to {} {} Peaks in GM12878\n'.format(p.target, p.assay))
ax.set_xlabel('Distance (bp)')
ax.set_xlim(-10000, 10000)
ax.text(0, -0.1, 'Audit Information\nWarning: {}\nCompliance: {}\nError: {}'.format(p.audit_warning, p.audit_not_compliant, p.audit_error), transform=fig.transFigure);
plt.tight_layout()

pdf.savefig(fig)
pdf.close()