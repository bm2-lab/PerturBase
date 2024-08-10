import scanpy as sc
import warnings, os, subprocess
warnings.filterwarnings('ignore')
import seaborn as sns
import pandas as pd
pd.set_option('display.max_colwidth', 100)
from scipy import stats
from statsmodels.stats import multitest
import numpy as np
import os
import sys
import pertpy as pt #type: ignore
import subprocess
import matplotlib.pyplot as plt
from collections import Counter


from sklearn.metrics.pairwise import cosine_similarity
def calCosine(Xtr, Xte):
    dat_cor = pd.DataFrame(cosine_similarity(Xte,Xtr)) 
    dat_cor.columns = Xtr.index
    dat_cor.index = Xte.index
    return dat_cor

def doPCA():
    adata = sc.read_h5ad('filtered.h5ad')
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata)
    adata.write_h5ad('rawPCA.h5ad')

    mydict = {}
    distance = pt.tl.Distance('edistance', 'X_pca')
    pertGenes = adata.obs['gene'].unique()
    for pertGene in pertGenes:
        if pertGene == 'CTRL': continue
        X = adata.obsm['X_pca'][adata.obs['gene'] == 'CTRL']
        Y = adata.obsm['X_pca'][adata.obs['gene'] == pertGene]
        result = distance(X, Y)
        mydict[pertGene] = result
    dat = pd.DataFrame({'perturbation': mydict.keys(), 'distance': mydict.values()})
    dat.sort_values(by='distance', ascending=False, inplace=True)
    dat.to_csv('perturbationClass.tsv', sep='\t', index=False)

def doEdist():
    if not os.path.isdir('figures'): os.makedirs('figures')
    adata2 = sc.read_h5ad('rawPCA.h5ad')

    perturbationClass = pd.read_csv('perturbationClass.tsv', sep='\t')
    strongP = list(perturbationClass['perturbation'][:25])
    strongP1 = list(perturbationClass['perturbation'])
    if 'CTRL' not in strongP: strongP.append('CTRL')
    if 'CTRL' not in strongP1: strongP1.append('CTRL')

    adata3 = adata2[adata2.obs['gene'].isin(strongP)]
    adata4 = adata2[adata2.obs['gene'].isin(strongP1)] 
    distance = pt.tl.Distance('edistance', 'X_pca')
    df = distance.pairwise(adata3, groupby="gene", verbose=True)
    df1 = distance.pairwise(adata4, groupby="gene", verbose=True)
    df = df / df.max().max()
    df1 = df1 / df1.max().max()
    df.columns.name = '';  df.index.name = ''
    df1.columns.name = '';  df1.index.name = ''
    df = df.round(5)
    df1 = df1.round(5)
    df.to_csv('corBetweenPerturb/edistance_order.tsv', sep='\t', index=True, header=True)
    df.to_csv('corBetweenPerturb/corEdistance_order.tsv', sep='\t', index=True, header=True)
    df1.to_csv('corBetweenPerturb/edistance.tsv', sep='\t', index=True, header=True)
    df1.to_csv('corBetweenPerturb/corEdistance.tsv', sep='\t', index=True, header=True)


        

def CorE():
    perturbationClass = pd.read_csv('perturbationClass.tsv', sep='\t')
    strongP = list(perturbationClass['perturbation'][:25])
    strongP1 = list(perturbationClass['perturbation'])

    if 'CTRL' not in strongP: strongP.append('CTRL')
    if 'CTRL' not in strongP1: strongP1.append('CTRL')
    adata = sc.read_h5ad('rawPCA.h5ad')
    dat = pd.DataFrame(adata.obsm['X_pca'])
    df = dat.groupby(list(adata.obs['gene']), axis=0).mean()
    
    df1 = df.loc[strongP, :]
    df = df.loc[strongP1, :]
    cor_dat1 = calCosine(df1, df1)
    cor_dat = calCosine(df, df)
    cor_dat = cor_dat.round(5)
    cor_dat1 = cor_dat1.round(5)
    if not os.path.isdir('corBetweenPerturb'): os.makedirs('corBetweenPerturb')
    cor_dat1.to_csv('corBetweenPerturb/corExp_order.tsv', sep='\t', header=True, index=True) 
    cor_dat.to_csv('corBetweenPerturb/corExp.tsv', sep='\t', header=True, index=True) 


def ClusterDistribution(filein, fileout):
    adata1 = sc.read_h5ad('rawPCA.h5ad')
    count1 = pd.crosstab(index=adata1.obs["leiden"], columns=adata1.obs["gene"])
    count1.to_csv(filein, sep='\t', index=True)
    
    dat = pd.read_csv(filein, sep='\t', index_col=0)
    pertGenes = [i for i in dat.columns]
    columns = [str(i) for i in dat.index] + ['allPvalue', 'allScore']
    results = pd.DataFrame(1, columns=columns, index=pertGenes)
    for gene in pertGenes:
        for cluster in dat.index:
            sumGene = dat[gene].sum()
            sumCTRL = dat['CTRL'].sum()
            sumCluster = dat.loc[cluster, :].sum()
            clusterGene = dat.loc[cluster, gene]
            clusterCTRL = dat.loc[cluster, 'CTRL']
            if  clusterGene / sumCluster <= 0.1: continue
            if clusterGene == 0 and clusterCTRL == 0: continue
            score, pvalue, *_  = stats.chi2_contingency([[sumGene, clusterGene], [sumCTRL, clusterCTRL]])
            results.loc[gene, str(cluster)] = pvalue
    for gene in pertGenes:
        clusterGene = []
        clusterCTRL = []
        for i, j in zip(dat[gene], dat['CTRL']):
            if i +j != 0:
                clusterGene.append(i)
                clusterCTRL.append(j)
        score, pvalue, *_  = stats.chi2_contingency([clusterGene, clusterCTRL])
        results.loc[gene, 'allPvalue'] = pvalue
        results.loc[gene, 'allScore'] = score
    reject, p_adjusted, _, _ = multitest.multipletests(results['allPvalue'], method='fdr_bh')
    results['allPvalue_adjust'] = p_adjusted
    results.fillna(1, inplace=True)
    results.sort_values('allScore', ascending=False, inplace=True)
    results = results.round(5)
    results.to_csv(fileout, sep='\t', header=True, index=True)
    return results

    
if __name__ == '__main__':
    path = sys.argv[1]
    species = sys.argv[2]
    # path = '/home/sirm/graduation_design/snakemake_folder/Demo_data/ATAC/fragment_file'
    # species = 'Hs'
    
    os.chdir(path)
    os.mkdir('corBetweenPerturb')
    
    doPCA()
    print ('Edist')
    doEdist()
    print('CorE')
    CorE()
    print('ClusterDistribution')
    ClusterDistribution('clusterGene_Distribution.tsv', 'clusterDistribution_pvalue.tsv')