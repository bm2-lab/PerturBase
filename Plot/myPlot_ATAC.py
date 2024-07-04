import scanpy as sc
import warnings, os, subprocess
warnings.filterwarnings('ignore')
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
pd.set_option('display.max_colwidth', 100)
import matplotlib.colors as mcolors
from scipy import stats
from statsmodels.stats import multitest
from adjustText import adjust_text
import numpy as np


plt.rcParams.update({'font.size': 20})   ### 改变字体大小
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
plt.rc('font',family='Times New Roman')
plt.rcParams.update({'font.size': 10})   ### 改变字体大小


from sklearn.metrics.pairwise import cosine_similarity
def calCosine(Xtr, Xte):
    dat_cor = pd.DataFrame(cosine_similarity(Xte,Xtr))  ###行是Xte, 列是Xtr
    dat_cor.columns = Xtr.index
    dat_cor.index = Xte.index
    return dat_cor

#### 预处理ATAC数据降维，降低数据的维度，方便后续画图
def getHVG():
    for dirName in dirNames:
        os.chdir(dirName)
        adata = sc.read_h5ad('raw.h5ad')
        sc.pp.highly_variable_genes(adata, subset=True, n_top_genes=1000)
        adata.write_h5ad('rawHVG.h5ad')


### 画最基本的QC图
def plotQC():
    for dirName in dirNames:
        os.chdir(dirName)
        adata2 = sc.read_h5ad('rawHVG.h5ad')
        tmp = adata2.obs['gene'].value_counts()
        g = sns.displot(tmp, legend=False)
        g.fig.set_size_inches(8, 6)
        plt.xlabel("Cell Numbers")
        plt.ylabel("Count")
        plt.title("Number of cells per perturbation")
        g.savefig('figures/sgRNA_quality1.png', transparent=False, dpi=300, bbox_inches='tight')

### 生成画plotSgRNAPlot需要的barcode文件
def genBarcode():
    for dirName in dirNames:
        os.chdir(dirName)
        adata = sc.read_h5ad('rawHVG.h5ad')
        def fun1(adata):
            barcode = adata.obs[['gene']]
            mylist1  =[]; mylist2 = []
            for cell, gene in zip(barcode.index, barcode['gene']):
                genes = gene.split(',')
                [mylist1.append(cell) for _ in genes]
                [mylist2.append(i) for i in genes]
            return mylist1, mylist2
        mylist1, mylist2 = fun1(adata)
        barcode = pd.DataFrame({'cell':mylist1, 'barcode':mylist2, "gene": mylist2}, index=mylist1) ### 保留组合扰动， scmageck rra软件会自动去除这些细胞，lr会保留这些细胞
        barcode.to_csv('barcode.txt', sep='\t', index=False, header=True)

## 画sgRNA的质控图
def plotSgRNAPlot():
    for dirName in dirNames:
        os.chdir(dirName)
        cmd = 'Rscript   /home/wzt/project/HC-CrisprScreen/myPlot_ATAC.r'
        subprocess.call(cmd, shell=True)


def gen_mpl_labels(
    adata, groupby, exclude=(), ax=None, adjust_kwargs=None, text_kwargs=None
):
    if adjust_kwargs is None:
        adjust_kwargs = {"text_from_points": False}
    if text_kwargs is None:
        text_kwargs = {}

    medians = {}

    for g, g_idx in adata.obs.groupby(groupby).groups.items():
        if g in exclude:
            continue
        medians[g] = np.median(adata[g_idx].obsm["X_umap"], axis=0)

    if ax is None:
        texts = [
            plt.text(x=x, y=y, s=k, **text_kwargs) for k, (x, y) in medians.items()
        ]
    else:
        texts = [ax.text(x=x, y=y, s=k, **text_kwargs) for k, (x, y) in medians.items()]

    adjust_text(texts, **adjust_kwargs)




### 画聚类图
def plotMixscape():
    if not os.path.isdir('figures'): os.makedirs('figures')
    perturbationClass = pd.read_csv('perturbationClass.tsv', sep='\t')
    count = pd.read_csv('clusterGene_Distribution.tsv', sep='\t', index_col=0)
    intersection = [i for i in  perturbationClass['perturbation'] if i in count.columns]
    if 'CTRL' not in intersection: intersection.append('CTRL')
    perturbationClass = perturbationClass[perturbationClass['perturbation'].isin(intersection)]
    count = count[intersection]
    strongP = list(perturbationClass['perturbation'][:25])
    if 'CTRL' not in strongP: strongP.append('CTRL')
    count = count[strongP]
    if len(count.index.unique()) >= 50:  ### 最多保留50个cluster
        a = [int(i) for i in count.index.unique()]
        a = sorted(a)[:50]
        count = count.loc[a, :]
    count1 = count / count.sum(axis=0)
    count1.sort_index(axis=1, inplace=True)
    count1.sort_index(axis=0, inplace=True)

#### enrichmap图
    cmap_custom = mcolors.LinearSegmentedColormap.from_list("custom_cmap", ["white", "red"])
    fig, ax1 = plt.subplots(figsize=(12, 8))
    sns.heatmap(count1, vmin=0, vmax=1, cmap = cmap_custom, linecolor='black', linewidths=1, annot=False, center=True)
    ax1.set_xlabel('Perturbation'); ax1.set_ylabel('Cluster')
    fig.savefig('figures/enrichmap.png', transparent=False, dpi=300, bbox_inches='tight')

def f_plotMixscape():
    for dirName in dirNames:
        os.chdir(dirName)
        plotMixscape()
            




def ClusterDistribution(filein, fileout):
    dat = pd.read_csv(filein, sep='\t', index_col=0)
    pertGenes = [i for i in dat.columns]
    columns = [str(i) for i in dat.index] + ['allPvalue', 'allScore']
    results = pd.DataFrame(1, columns=columns, index=pertGenes)
    ### 统计在每个cluste的卡方差异
    for gene in pertGenes:
        for cluster in dat.index:
            sumGene = dat[gene].sum()   ### 扰动的总个数
            sumCTRL = dat['CTRL'].sum()  ### 对照的总个数
            sumCluster = dat.loc[cluster, :].sum()  ### cluster的总个数
            clusterGene = dat.loc[cluster, gene]  ### 扰动在cluster的个数
            clusterCTRL = dat.loc[cluster, 'CTRL']  ### ctrl在cluster的个数
            if  clusterGene / sumCluster <= 0.1: continue   ### 小于20%，则不统计
            if clusterGene == 0 and clusterCTRL == 0: continue  ### 如果观察值都是0，去除
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
    results.to_csv(fileout, sep='\t', header=True, index=True)
    return results


### 画聚类后，与ctrl对比的差异统计检验
def f_ClusterDistribution():
    for dirName in dirNames:
        os.chdir(dirName)
        try:
            ClusterDistribution('clusterGene_Distribution.tsv', 'clusterDistribution_pvalue.tsv')
        except Exception as e:
            print (e); print (dirName)


def PlotClusterDistribution():
    perturbationClass = pd.read_csv('perturbationClass.tsv', sep='\t')
    count = pd.read_csv('clusterGene_Distribution.tsv', sep='\t', index_col=0)
    intersection = [i for i in  perturbationClass['perturbation'] if i in count.columns]
    if 'CTRL' not in intersection: intersection.append('CTRL')
    perturbationClass = perturbationClass[perturbationClass['perturbation'].isin(intersection)]
    count = count[intersection]
    strongP = list(perturbationClass['perturbation'][:25])
    if 'CTRL' not in strongP: strongP.append('CTRL')
    count = count[strongP]
    if len(count.index.unique()) >= 50:  ### 最多保留50个cluster
        a = [int(i) for i in count.index.unique()]
        a = sorted(a)[:50]
        count = count.loc[a, :]
    count1 = count / count.sum(axis=0)
    count1.sort_index(axis=1, inplace=True)
    count1.sort_index(axis=0, inplace=True)

    dat = pd.read_csv('clusterDistribution_pvalue.tsv', sep='\t', index_col=0).T
    dat = dat.loc[[str(i) for i in count1.index] + ['allPvalue', 'allScore', 'allPvalue_adjust'], count1.columns]
#### enrichmap图
    cmap_custom = mcolors.LinearSegmentedColormap.from_list("custom_cmap", ["white", "red"])
    fig, ax1 = plt.subplots(figsize=(12, 8))
    ax = sns.heatmap(count1, vmin=0, vmax=1, cmap = cmap_custom, linecolor='black', linewidths=1, annot=False, center=True)
    ax1.set_ylabel('Cluster'); ax1.set_xlabel('Perturbation')

    # 在p-values小于0.01的单元格上添加星号
    for i in range(count1.shape[0]):  ### 行
        for j in range(count1.shape[1]):  ### 列
            if dat.iloc[i, j] < 0.01:
                ax.text(j + 0.5, i + 0.5, "*", ha='center', va='center', color='black', fontsize=12)

    # 根据差异表达的样品在xlable的基因名上添加星号
    for i in range(count1.shape[1]):
        if dat.loc['allPvalue_adjust', count1.columns[i]] < 0.01:
            ax.text(i + 0.5, count1.shape[0]+ .5, "*", ha='center', va='center', color='red', fontsize=12)   
    fig.savefig('figures/functionEnrichmap.png', transparent=False, dpi=300, bbox_inches='tight')

##
def f_PlotClusterDistribution():
    for dirName in dirNames:
        os.chdir(dirName)
        PlotClusterDistribution()



### 使用hvg表达谱计算扰动之间的相关性
def plotCor1():
    plt.rcParams.update({'font.size': 15})   ### 改变字体大小
    for dirName in dirNames:
        os.chdir(dirName)
        try:
            perturbationClass = pd.read_csv('perturbationClass.tsv', sep='\t')
            count = pd.read_csv('clusterGene_Distribution.tsv', sep='\t', index_col=0)
            intersection = [i for i in  perturbationClass['perturbation'] if i in count.columns]
            if 'CTRL' not in intersection: intersection.append('CTRL')
            perturbationClass = perturbationClass[perturbationClass['perturbation'].isin(intersection)]
            strongP = list(perturbationClass['perturbation'][:25])
            strongP1 = list(perturbationClass['perturbation'])

            if 'CTRL' not in strongP: strongP.append('CTRL')
            if 'CTRL' not in strongP1: strongP1.append('CTRL')
            adata = sc.read_h5ad('rawPCA.h5ad')
            dat = pd.DataFrame(adata.obsm['X_pca'])
            df = dat.groupby(list(adata.obs['gene']), axis=0).mean()
            strongP = [i for i in strongP if i in list(df.index)]
            strongP1 = [i for i in strongP1 if i in list(df.index)]
            df1 = df.loc[strongP, :]
            df = df.loc[strongP1, :]
            cor_dat1 = calCosine(df1, df1)
            cor_dat = calCosine(df, df)
            if not os.path.isdir('corBetweenPerturb'): os.makedirs('corBetweenPerturb')
            cor_dat1.to_csv('corBetweenPerturb/corExp_order.tsv', sep='\t', header=True, index=True) ### 用过滤之后的数据
            cor_dat.to_csv('corBetweenPerturb/corExp.tsv', sep='\t', header=True, index=True) ### 没过滤30个细胞数量的要求


            #g = sns.clustermap(cor_dat1, cmap="vlag", vmin=-1, vmax=1, cbar_kws={"ticks": [-1, 0, 1]}, linewidths=0.5, linecolor='black', cbar=True, figsize=(12, 12))
            #g.ax_heatmap.set_title("Correlation Heatmap of Perturbations", pad = 200)
            #g.savefig("figures/perturbation_cor_exp.png", dpi=300, bbox_inches='tight', transparent=False)
        except Exception as e:
            print (e); print (dirName)



dirNames = ['/home/wzt/project/HC-CrisprScreen/poolSC_data/7Perturb-ATAC/PRJNA658075/PRJNA658075',   ###可以作为demo数据
            '/home/wzt/project/HC-CrisprScreen/poolSC_data/7Perturb-ATAC/PRJNA714243/Day6', 
            '/home/wzt/project/HC-CrisprScreen/poolSC_data/7Perturb-ATAC/PRJNA714243/Day9',
            '/home/wzt/project/HC-CrisprScreen/poolSC_data/7Perturb-ATAC/PRJNA714243/Day21', 
            '/home/wzt/project/HC-CrisprScreen/poolSC_data/7Perturb-ATAC/PRJNA478043/Bcell_experiment1', 
            '/home/wzt/project/HC-CrisprScreen/poolSC_data/7Perturb-ATAC/PRJNA478043/Bcell_experiment2', 
            '/home/wzt/project/HC-CrisprScreen/poolSC_data/7Perturb-ATAC/PRJNA478043/Keratinocyte', 
            '/home/wzt/project/HC-CrisprScreen/poolSC_data/7Perturb-ATAC/PRJNA893678/PRJNA893678']



'''
/home/wzt/anaconda3/bin/python   myPlot_ATAC.py
'''

    
if __name__ == '__main__':
    print ('hello, world')
    #getHVG()
    #genBarcode()
    
    #plotQC()   ### 第一部分图
    #plotSgRNAPlot()
    #f_plotMixscape()   ###画聚类热图

    #f_ClusterDistribution()  ### 计算热图的pvalue  ### 第二部分图
    #f_PlotClusterDistribution()  ### 第二部分图

    plotCor1()