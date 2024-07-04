from myUtil import *
import scanpy as sc
import warnings
import anndata as ad
from scipy import sparse
from scipy import stats
import signal
import multiprocessing
import pickle
from tqdm import tqdm
warnings.filterwarnings('ignore')

### #### 


def cal_result(dirName):
    os.chdir(dirName)
    if not os.path.isfile('mixscape_filter.h5ad'): return
    adata1 = sc.read_h5ad('mixscape_filter.h5ad')
    adata1.uns['log1p']["base"] = None
    sc.pp.highly_variable_genes(adata1, n_top_genes=10000)
    all_hvgs = adata1[:, adata1.var['highly_variable']].var_names  ### 直接保留的基因

    perts = adata1.obs['gene'].unique()
    perts = [i for i in perts if i !='CTRL']
    with open('missingDEG_ratio1.tsv', 'w') as fout:
        fout.write('Pertb\tMissingNums\ttotalNums\tMissingRatio\n')
        for pert in tqdm(perts):
            mylist = []
            adata3 = adata1[adata1.obs['gene'].isin(['CTRL', pert])]
            adata3.uns['log1p']["base"] = None
            sc.pp.highly_variable_genes(adata3, subset=True, n_top_genes=4000)
            sc.tl.rank_genes_groups(adata3, 'gene', groups=[pert], reference='CTRL')
            for pvalue, geneName, logfoldchanges in zip(adata3.uns['rank_genes_groups']['pvals_adj'], adata3.uns['rank_genes_groups']['names'],  adata3.uns['rank_genes_groups']['logfoldchanges']):
                if pvalue[0] <= 0.05 and abs(logfoldchanges[0]) >=1:
                    mylist.append(geneName[0])
            tmp = [i for i in mylist if i not in all_hvgs]
            if len(mylist) !=0:
                fout.write('{}\t{}\t{}\t{:.3f}\n'.format(pert, len(tmp), len(mylist) ,len(tmp)/ len(mylist)))
            else:
                fout.write('{}\t{}\t{}\t{:.3f}\n'.format(pert, len(tmp), 0, 0))


### 分析数据
def analysis1():
    mylist = []
    dat1 = pd.read_excel('/home/wzt/project/HC-CrisprScreen/poolSC_data/Sheet10_modif.xlsx')
    dat1 = dat1[((dat1['QC'] == 'Pass') & (dat1['Modality'] !='ATAC')  & (dat1['Perturbation Type'] == 'Genetic'))]
    dat1.sort_values('Filter_PerturbationNums', ascending=True, inplace=True)
    for dirName in tqdm(dat1['Datapath']):
        os.chdir(dirName)
        if os.path.isfile('missingDEG_ratio.tsv'):
            dat = pd.read_csv('missingDEG_ratio.tsv', sep='\t')
            dat['PertNums']= dat.shape[0]
            if dat.shape[0] >=10000:continue
            mylist.append(dat)
    datAll = pd.concat(mylist)




#### 
if __name__== '__main__':
    dat = pd.read_excel('/home/wzt/project/HC-CrisprScreen/poolSC_data/Sheet10_modif.xlsx')
    dat = dat[((dat['QC'] == 'Pass') & (dat['Modality'] !='ATAC'))]
    dat.sort_values('Filter_PerturbationNums', ascending=True, inplace=True)
    for dirName in tqdm(dat['Datapath']):
        print (dirName)
        cal_result(dirName)


