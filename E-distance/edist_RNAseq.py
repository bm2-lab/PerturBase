from myUtil import *
import scanpy as sc
import warnings
import anndata as ad
warnings.filterwarnings('ignore')
import matplotlib.pyplot as plt
import seaborn as sns
from MixScape2 import classGene

import pickle
import pertpy as pt # type: ignore
import scanpy as sc
from seaborn import clustermap

def doEdist(dirName):
    os.chdir(dirName)
    if os.path.isfile('mixscape_hvg_filter.h5ad') and os.path.isfile('perturbationClass.tsv'):
        adata = sc.read_h5ad('mixscape_hvg_filter_subset.h5ad')  ### 高变异基因
        perturbationClass = pd.read_csv('perturbationClass.tsv', sep='\t')
        strongP = list(perturbationClass['perturbation'][:25])
        if 'CTRL' not in strongP: strongP.append('CTRL')
        adata2 = adata[adata.obs['gene'].isin(strongP)]

        perturbationClass = pd.read_csv('perturbationClass.tsv', sep='\t')
        strong = list(perturbationClass[perturbationClass['perturbationClass'] == 'strongPerturbation']['perturbation'])
        weak = list(perturbationClass[perturbationClass['perturbationClass'] == 'weakPerturbation']['perturbation'])

        distance = pt.tl.Distance('edistance', 'X_pca')
        if not os.path.isfile('corBetweenPerturb/corEdistance_sub.XXX'):
            df = distance.pairwise(adata2, groupby="gene", verbose=True)
            df = df / df.max().max()
            df.to_csv('corBetweenPerturb/corEdistance_sub.tsv', sep='\t', header=True, index=True, index_label='')

            df1 = distance.pairwise(adata, groupby="gene", verbose=True)  ### 所有的数据提供下载
            df1 = df1 / df1.max().max()
            df1.to_csv('corBetweenPerturb/corEdistance.tsv', sep='\t', header=True, index=True, index_label='')
        else:
            df = pd.read_csv('corBetweenPerturb/corEdistance_sub.tsv', sep='\t', header=0, index_col=0)  
        plt.rcParams.update({'font.size': 15})   ### 改变字体大小
        row_colors = ['r' if i in strong else 'g' for i in df.index]
        g = clustermap(df, robust=False, figsize=(15,15), row_colors = row_colors, cbar_kws={"ticks": [0, 0.5, 1]})   ### 画edist距离，聚类热图
        
        reordered_rows = g.dendrogram_row.reordered_ind
        reordered_cols = g.dendrogram_col.reordered_ind
        cor_dat2 = df.iloc[reordered_rows, reordered_cols]
        cor_dat2.to_csv('corBetweenPerturb/corEdistance_order.tsv', sep='\t', header=True, index=True)
        with open('corBetweenPerturb/corEdistance_row.linkage.pkl', 'wb') as fout:
            pickle.dump(obj=g.dendrogram_row.linkage, file=fout)
        with open('corBetweenPerturb/corEdistance_col.linkage.pkl', 'wb') as fout:
            pickle.dump(obj=g.dendrogram_col.linkage, file=fout)    
        
        #legend_labels = ['strong', 'weak']; legend_colors = ['r', 'g']
        #legend_elements = [plt.Line2D([0], [0], marker='o', color='black', markerfacecolor=color, markersize=15, label=label,linewidth=0,) for label, color in zip(legend_labels, legend_colors)]
        #g.ax_heatmap.legend(legend_elements, legend_labels, loc='upper left', title='PerturbationClass', bbox_to_anchor=(-0.2, 1.2), frameon=False)
        #g.ax_heatmap.set_title("Distance Heatmap of Perturbations", pad = 200)
        #g.savefig('figures/clustermap.png', dpi=300, bbox_inches='tight', transparent=False)


def f_doEdist():
    mylist = []
    dat = pd.read_csv('/home/wzt/project/HC-CrisprScreen/results/dataInfo.tsv', sep='\t', header=0)
    for dirName, species in zip(dat["data"], dat["species"]):
        if dirNames and dirName not in dirNames:
            continue
        else:
            mylist.append(dirName)    
    myPool(doEdist, mylist, 5)

'''
反应扰动之间的关系，比如在/home/wzt/project/HC-CrisprScreen/poolSC_data/8ECCITE-seq/PRJNA641353/ECCITE数据集中，
STAT1, IFNGR1, JAK2, IFNGR1四个扰动相似性很高，具有相同的扰动效果，而一些扰动和CTRL很相似，被定义为weak perturbation。
'''
### 跑7个demo
dirNames  = ['/home/wzt/project/HC-CrisprScreen/poolSC_data/8ECCITE-seq/PRJNA641353/ECCITE', 
             #'/home/wzt/project/HC-CrisprScreen/poolSC_data/16TAP-seq/PRJNA559094/SCREEN_chr8', 
             #'/home/wzt/project/HC-CrisprScreen/poolSC_data/16TAP-seq/PRJNA559094/SCREEN_chr11', 
             #'/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb-seq/PRJNA587707/sciPlex3_A549', 
             #'/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb-seq/PRJNA587707/sciPlex3_K562', 
             #'/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb-seq/PRJNA628589/GSE149215',    ### failed
             #'/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb-seq/PRJNA831566/K562_essential',
             #'/NFS_home/NFS_home_2/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb-seq/PRJNA587707/sciPlex3_A549'
             ]


if __name__ == '__main__':
    print ('hello, world')
    f_doEdist()
