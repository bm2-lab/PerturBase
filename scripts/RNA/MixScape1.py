from myUtil import *
import scanpy as sc
import warnings
warnings.filterwarnings('ignore')
from denoising import classGene
import sys
import pertpy as pt
import scanpy as sc


def calDistance(dirName):
    os.chdir(dirName)
    if not os.path.isfile('mixscape_hvg.h5ad'): return
    adata1 = sc.read_h5ad('mixscape_hvg.h5ad')
    strongP, weakP, *_ = classGene(adata1, threshold=.2)  

    mydict = {}
    adata= sc.read_h5ad('mixscape_hvg_filter.h5ad')
    distance = pt.tl.Distance('edistance', 'X_pca')
    pertGenes = adata.obs['gene'].unique()
    X = adata.obsm['X_pca'][adata.obs['gene'] == 'CTRL']
    for pertGene in tqdm(pertGenes):
        if pertGene == 'CTRL': continue
        Y = adata.obsm['X_pca'][adata.obs['gene'] == pertGene]
        result = distance(X, Y)
        mydict[pertGene] = result
    dat = pd.DataFrame({'perturbation': mydict.keys(), 'distance': mydict.values()})
    tmp = ['strongPerturbation' if i in strongP else 'weakPerturbation' for i in dat['perturbation']]
    dat['perturbationClass'] = tmp
    dat.sort_values(by='distance', ascending=False, inplace=True)

### count NP and SP
    tmp1 = []; tmp2 = []
    for pertGene in dat['perturbation']:
        tmp = adata1.obs[adata1.obs['gene'] == pertGene]
        NP = tmp['mixscape_class_global'].value_counts().get('NP', 0)
        KO = tmp['mixscape_class_global'].value_counts().get('KO', 0)
        tmp1.append(NP); tmp2.append(KO)
    dat['KO'] = tmp2
    dat['NP'] = tmp1
    dat['EscapeRatio'] = dat['NP'] / (dat['NP'] + dat['KO'])
    dat.to_csv('perturbationClass.tsv', sep='\t', index=False)


def pertCountPerClass(dirName):
    os.chdir(dirName)
    print (dirName)
    #if os.path.isfile('clusterDistribution_beforeMix.tsv'): continue
    adata1 = sc.read_h5ad('mixscape_hvg.h5ad')
    adata2 = sc.read_h5ad('mixscape_hvg_filter.h5ad')
    count1 = pd.crosstab(index=adata1.obs["leiden"], columns=adata1.obs["gene"])
    count2 = pd.crosstab(index=adata2.obs["leiden"], columns=adata2.obs["gene"])
    count1.to_csv('clusterDistribution_beforeMix.tsv', sep='\t', index=True)
    count2.to_csv('clusterDistribution_afterMix.tsv', sep='\t', index=True)


if __name__ == '__main__':
    path = sys.argv[1]
    species = sys.argv[2]
    calDistance(path)
    pertCountPerClass(path)
