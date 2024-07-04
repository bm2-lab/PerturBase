from myUtil import *
import scanpy as sc
import warnings
import anndata as ad
warnings.filterwarnings('ignore')
import matplotlib.pyplot as plt
import seaborn as sns
from MixScape2 import classGene

import pertpy as pt
import scanpy as sc
from seaborn import clustermap


def doPCA():
    for dirName in dirNames:
        print (dirName)
        os.chdir(dirName)
        if not os.path.isfile('rawPCA.h5ad'):
            adata = sc.read_h5ad('raw.h5ad')  ### 使用原始的数据直接进行PCA，保留最大的信息
            sc.pp.pca(adata)
            sc.pp.neighbors(adata)
            sc.tl.umap(adata)
            sc.tl.leiden(adata)
            adata.write_h5ad('rawPCA.h5ad')
        else:
            adata = sc.read_h5ad('rawPCA.h5ad')

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
    for dirName in dirNames:
        os.chdir(dirName)
        if not os.path.isdir('figures'): os.makedirs('figures')
        adata2 = sc.read_h5ad('rawPCA.h5ad')  ### 高变异基因

        perturbationClass = pd.read_csv('perturbationClass.tsv', sep='\t')
        count = pd.read_csv('clusterGene_Distribution.tsv', sep='\t', index_col=0)
        intersection = [i for i in  perturbationClass['perturbation'] if i in count.columns]
        if 'CTRL' not in intersection: intersection.append('CTRL')
        perturbationClass = perturbationClass[perturbationClass['perturbation'].isin(intersection)]
        strongP = list(perturbationClass['perturbation'][:25])
        strongP1 = list(perturbationClass['perturbation'][:100])
        if 'CTRL' not in strongP: strongP.append('CTRL')
        if 'CTRL' not in strongP1: strongP1.append('CTRL')

        adata3 = adata2[adata2.obs['gene'].isin(strongP)]
        adata4 = adata2[adata2.obs['gene'].isin(strongP1)]  ### csv表格保留的样品更多
        distance = pt.tl.Distance('edistance', 'X_pca')
        df = distance.pairwise(adata3, groupby="gene", verbose=True)
        df1 = distance.pairwise(adata4, groupby="gene", verbose=True)
        df = df / df.max().max()
        df1 = df1 / df1.max().max()
        df.columns.name = '';  df.index.name = ''
        df1.columns.name = '';  df1.index.name = ''
        df.to_csv('corBetweenPerturb/edistance_order.tsv', sep='\t', index=True, header=True)
        df.to_csv('corBetweenPerturb/corEdistance_order.tsv', sep='\t', index=True, header=True)
        df1.to_csv('corBetweenPerturb/edistance.tsv', sep='\t', index=True, header=True)
        df1.to_csv('corBetweenPerturb/corEdistance.tsv', sep='\t', index=True, header=True)




        plt.rcParams.update({'font.size': 15})   ### 改变字体大小
        g = clustermap(df, robust=False, figsize=(15,15), cbar_kws={"ticks": [0, 0.5, 1]})   ### 画edist距离，聚类热图
        g.ax_heatmap.set_title("Distance Heatmap of Perturbations", pad = 200)
        g.savefig('figures/clustermap.png', dpi=300, bbox_inches='tight', transparent=False)


'''
反应扰动之间的关系，比如在/home/wzt/project/HC-CrisprScreen/poolSC_data/8ECCITE-seq/PRJNA641353/ECCITE数据集中，
STAT1, IFNGR1, JAK2, IFNGR1四个扰动相似性很高，具有相同的扰动效果，而一些扰动和CTRL很相似，被定义为weak perturbation。
'''

dirNames = ['/home/wzt/project/HC-CrisprScreen/poolSC_data/7Perturb-ATAC/PRJNA658075/PRJNA658075', 
            '/home/wzt/project/HC-CrisprScreen/poolSC_data/7Perturb-ATAC/PRJNA714243/Day6', 
            '/home/wzt/project/HC-CrisprScreen/poolSC_data/7Perturb-ATAC/PRJNA714243/Day9',
            '/home/wzt/project/HC-CrisprScreen/poolSC_data/7Perturb-ATAC/PRJNA714243/Day21', 
            '/home/wzt/project/HC-CrisprScreen/poolSC_data/7Perturb-ATAC/PRJNA478043/Bcell_experiment1', 
            '/home/wzt/project/HC-CrisprScreen/poolSC_data/7Perturb-ATAC/PRJNA478043/Bcell_experiment2', 
            '/home/wzt/project/HC-CrisprScreen/poolSC_data/7Perturb-ATAC/PRJNA478043/Keratinocyte', 
            '/home/wzt/project/HC-CrisprScreen/poolSC_data/7Perturb-ATAC/PRJNA893678/PRJNA893678']



if __name__ == '__main__':
    print ('hello, world')
    #doPCA()
    doEdist()
