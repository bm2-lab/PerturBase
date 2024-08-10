import scanpy as sc, pandas as pd, anndata as ad
import warnings, os
warnings.filterwarnings('ignore')
from matplotlib import pyplot as plt

pd.set_option('display.float_format', lambda x: '%.4f' % x)

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
plt.rcParams.update({'font.size': 20})
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
plt.rc('font',family='Times New Roman')
plt.rcParams.update({'font.size': 10})

os.chdir(dirName)

dat = sc.read_h5ad('mixscape_hvg.h5ad')


geneList = ['SLC7A11', 'CYTIP', 'CXCL3', 'TNFAIP6', 'ADAMDEC1', 'IL8', 'CDK5RAP2', 'PGD', 'TXN', 'TXNRD1', 'KYNU', 'EREG', 'CXCL2', 'DCSTAMP', 'TMEM255A', 'CCL4', 'SLC12A8',
           'SOD2', 'SH3BP5', 'SAT1', 'PDE4DIP', 'GCLM', 'FTH1', 'TDP2', 'CCL8', 'GK', 'LRP12', 'IDO1', 'PLA2G7', 'CYP27A1', 'CES1', 'ME1', 'DUSP6', 'C5AR1', 'CCL3', 'CDKN1A',
           'RIT1', 'TCHH', 'SLAMF7', 'ABCA1', 'MARCKSL1', 'CEBPB', 'LUCAT1', 'EBI3', 'INSIG1', 'PLEK', 'TMEM38B', 'TTLL4', 'IL32', 'COL6A2', 'CXCL9', 'PPP1R15A', 'MSC', 'RAI14', 
           'PIM1', 'SQSTM1', 'NEAT1', 'C1QA', 'NINJ1', 'CCL13', 'KLF9', 'CD274', 'CAV1', 'SERPINE2']
geneList = [i for i in geneList if i in dat.var_names]

dat2 = dat[dat.obs['mixscape_class'].isin(['CTRL'])][:100]
dat3 = dat[dat.obs['mixscape_class'].isin(['CUL3 NP'])][:100]
dat4 = dat[dat.obs['mixscape_class'].isin(['CUL3 KO'])]
dat1 =  ad.concat([dat2, dat3, dat4])
dat1.obs['mixscape_class'] = dat1.obs['mixscape_class'].cat.reorder_categories(['CTRL', 'CUL3 NP', 'CUL3 KO'])
sc.pp.scale(dat1, max_value=2)
sc.pl.heatmap(dat1, geneList, groupby='mixscape_class', swap_axes=True, show_gene_labels=True)

    
