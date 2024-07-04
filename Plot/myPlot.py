import scanpy as sc, numpy as np
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
from myUtil import *
import pickle

try:
    import pertpy as pt
except:
    pass

import warnings
warnings.filterwarnings('ignore')

plt.rcParams.update({'font.size': 20})   ### 改变字体大小
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
plt.rc('font',family='Times New Roman')


from sklearn.metrics.pairwise import cosine_similarity
def calCosine(Xtr, Xte):
    dat_cor = pd.DataFrame(cosine_similarity(Xte,Xtr))  ###行是Xte, 列是Xtr
    dat_cor.columns = Xtr.index
    dat_cor.index = Xte.index
    return dat_cor


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





### 画最基本的QC图
def plotQC(dirName):
    os.chdir(dirName)
    if os.path.isfile('mixscape_hvg_filter.h5ad'):
        if not os.path.isdir('figures/violinfigures'): os.makedirs('figures/violinfigures')
        adata2 = sc.read_h5ad('mixscape_hvg.h5ad')

        # tmp = adata2.obs['gene'].value_counts()
        # g = sns.displot(tmp, legend=False)
        # g.figure.set_size_inches(8, 6)
        # plt.xlabel("Cell Number")
        # plt.ylabel("Count")
        # plt.title("Number of cells per perturbation")
        # g.savefig('figures/sgRNA_quality1.png', transparent=False,  dpi=300, bbox_inches='tight')

        if adata2.shape[0] >= 10000:        #### 第一张图
            sc.pp.subsample(adata2, n_obs=10000, random_state=42)
        # sc.pl.violin(adata2, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
        #     jitter=0.4, multi_panel=True, save='figures/afterQC1.png', stripplot=False, dpi = 300)
            
        g = sc.pl.violin(adata2, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                 jitter=0.4, multi_panel=True, save=False, stripplot=False, show=False)
        g.figure.set_size_inches(20, 5)
        g.savefig('figures/afterQC1.png', transparent=False,  dpi=800, bbox_inches='tight')
                

def f_plotQC():
    dat = pd.read_csv('/home/wzt/project/HC-CrisprScreen/results/dataInfo.tsv', sep='\t', header=0)
    for dirName, species in zip(dat["data"], dat["species"]):
        if dirNames and dirName not in dirNames: continue
        print (dirName)
        plotQC(dirName)



def plotSgRNAPlot(dirName):
    os.chdir(dirName)
    if os.path.isfile('mixscape_hvg.h5ad'):
        cmd = 'Rscript   /home/wzt/project/HC-CrisprScreen/myPlot.r'
        subprocess.call(cmd, shell=True)
    else: return


## 画sgRNA的质控图
def f_plotSgRNAPlot():
    mylist = []
    dat = pd.read_csv('/home/wzt/project/HC-CrisprScreen/results/dataInfo.tsv', sep='\t', header=0)
    for dirName, species in zip(dat["data"], dat["species"]):
        if dirNames and dirName not in dirNames: continue
        plotSgRNAPlot(dirName)
        #mylist.append(dirName)
    #myPool(plotSgRNAPlot, mylist, processes=1)


def enrichmentMap1(adata):
    count = pd.crosstab(index=adata.obs["leiden"], columns=adata.obs["gene"])
    count1 = count / count.sum(axis=0)
    cmap_custom = mcolors.LinearSegmentedColormap.from_list("custom_cmap", ["white", "red"])
    fig, ax = plt.subplots(figsize=(12, 8))
    sns.heatmap(count1, vmin=0, vmax=1, cmap = cmap_custom, linecolor='black', linewidths=1, annot=False, center=True)
    ax.set_xlabel('Perturbation'); ax.set_ylabel('Cluster')
    fig.savefig('figures/enrichmap_beforeMix.png', transparent=False, dpi=300, bbox_inches='tight')

def enrichmentMap2(adata):
    count = pd.crosstab(index=adata.obs["leiden"], columns=adata.obs["gene"])
    count1 = count / count.sum(axis=0)

    cmap_custom = mcolors.LinearSegmentedColormap.from_list("custom_cmap", ["white", "red"])
    fig, ax = plt.subplots(figsize=(12, 8))
    sns.heatmap(count1, vmin=0, vmax=1, cmap = cmap_custom, linecolor='black', linewidths=1, annot=False, center=True)
    ax.set_xlabel('Perturbation'); ax.set_ylabel('Cluster')
    fig.savefig('figures/enrichmap_afterMix.png', transparent=False, dpi=300, bbox_inches='tight')


###画 mixscape数据清洗前后的对比图
def plotMixscape():
    if not os.path.isdir('figures'): os.makedirs('figures')
    perturbationClass = pd.read_csv('perturbationClass.tsv', sep='\t')
    strongP = list(perturbationClass['perturbation'][:25])
    if 'CTRL' not in strongP: strongP.append('CTRL')
    
    adata = sc.read_h5ad('mixscape_hvg.h5ad')
    adata_tmp1 = adata[adata.obs['gene'].isin(strongP)]

#### 画mixscape pca聚类图
    adata.obs['tmp1'] = adata.obs['leiden'].apply(lambda x : x if int(x) < 25 else '26')
    g = sc.pl.umap(adata, color=['tmp1'], legend_fontsize='medium', legend_loc='on data', 
        frameon=True, title='',save=False, return_fig=True, groups=[str(i) for i in range(25)], na_in_legend=False)
    g.savefig('figures/pcaBeforeMixScape_leiden.png', transparent=False, dpi=300, bbox_inches='tight')
    
    adata.obs['tmp2'] = adata.obs['gene'].apply(lambda x: x if x in list(strongP) else 'Other')
    ax = sc.pl.umap(adata, color=['tmp2'], show=False, legend_loc=None, frameon=True, title='')
    gen_mpl_labels(adata, "tmp2",
    exclude=("None",),  # This was before we had the `nan` behaviour
    ax=ax,
    adjust_kwargs=dict(arrowprops=dict(arrowstyle='-', color='black')),
    text_kwargs=dict(fontsize=8))
    fig = ax.get_figure()
    fig.tight_layout()
    plt.show()
    fig.savefig('figures/pcaBeforeMixScape_gene.png', transparent=False, dpi=300, bbox_inches='tight')

    ### 过滤之后
    adata = sc.read_h5ad('mixscape_hvg_filter.h5ad')     
    adata_tmp2 = adata[adata.obs['gene'].isin(strongP)]       
    adata.obs['tmp1'] = adata.obs['leiden'].apply(lambda x : x if int(x) < 25 else '26')
    g = sc.pl.umap(adata, color=['tmp1'], legend_fontsize='medium', legend_loc='on data', 
        frameon=True, title='',save=False, return_fig=True, groups=[str(i) for i in range(25)], na_in_legend=False)
    g.savefig('figures/pcaMixScape_leiden.png', transparent=False, dpi=300, bbox_inches='tight')

    adata.obs['tmp2'] = adata.obs['gene'].apply(lambda x: x if x in list(strongP) else 'Other')
    ax = sc.pl.umap(adata, color=['tmp2'], show=False, legend_loc=None, frameon=True, title='')
    gen_mpl_labels(adata, "tmp2",
    exclude=("None",),  # This was before we had the `nan` behaviour
    ax=ax,
    adjust_kwargs=dict(arrowprops=dict(arrowstyle='-', color='black')),
    text_kwargs=dict(fontsize=8))
    fig = ax.get_figure()
    fig.tight_layout()
    plt.show()
    fig.savefig('figures/pcaMixScape_gene.png', transparent=False, dpi=300, bbox_inches='tight')


#### enrichmap图
    if len(adata_tmp1.obs['leiden'].unique()) >= 50:  ### 最多保留50个cluster
        a = [int(i) for i in adata_tmp1.obs['leiden'].unique()]
        a = sorted(a)[:50]
        adata_tmp1 = adata_tmp1[adata_tmp1.obs['leiden'].isin([str(i) for i in a])]
    enrichmentMap1(adata_tmp1)  ### 过滤之前

    if len(adata_tmp2.obs['leiden'].unique()) >= 50:
        a = [int(i) for i in adata_tmp2.obs['leiden'].unique()]
        a = sorted(a)[:50]
        adata_tmp2 = adata_tmp2[adata_tmp2.obs['leiden'].isin([str(i) for i in a])]
    enrichmentMap2(adata_tmp2)  ### 过滤之后

def f_plotMixscape():
    dat = pd.read_csv('/home/wzt/project/HC-CrisprScreen/results/dataInfo.tsv', sep='\t', header=0)
    for dirName, species in zip(dat["data"], dat["species"]):
        if dirNames and dirName not in dirNames: continue
        os.chdir(dirName)
        if os.path.isfile('mixscape_hvg.h5ad'):
            print (dirName)
            plotMixscape()

### 画scMageck的 phenotype的 气泡图
def doPhenotype():
    perturbationClass = pd.read_csv('perturbationClass.tsv', sep='\t')
    strongP = list(perturbationClass['perturbation'][:25])
    dat = pd.read_csv('scMageCK/RRA.txt', sep='\t')
    dat = dat[dat['Perturbation'].isin(strongP)]
    if dat.shape[0] == 0:
        fig, ax = plt.subplots(figsize=(14, 8))
        plt.text(0.2, 0.5, 'NO Enrichment Result for all of perturbations', fontsize=22)
        plt.show()
        fig.savefig("figures/Phenotype.png", transparent=False, dpi=300, bbox_inches='tight')
    else:
        fig, ax = plt.subplots(figsize=(14, 8))
        sns.scatterplot(data=dat, x='GeneSet', y='Perturbation', size='-log10(FDR)', hue='Selection')
        plt.xticks(rotation=90)
        plt.legend(loc=2, bbox_to_anchor=(1, 1))
        plt.show()
        fig.savefig("figures/Phenotype.png", transparent=False, dpi=300, bbox_inches='tight')



def f_doPhenotype():
    dat = pd.read_csv('/home/wzt/project/HC-CrisprScreen/results/dataInfo.tsv', sep='\t', header=0)
    for dirName, species in zip(dat["data"], dat["species"]):
        if dirNames and dirName not in dirNames: continue
        os.chdir(dirName)
        if os.path.isfile('mixscape_hvg_filter.h5ad') and os.path.isfile('scMageCK/RRA.txt'):
            if not os.path.isdir('figures'):  os.makedirs('figures')
            try:
                doPhenotype()
            except Exception as e:
                print (e); print (dirName)


def ClusterDistribution(filein, fileout):
    if os.path.isfile(fileout): return
    dat = pd.read_csv(filein, sep='\t', index_col=0)
    pertGenes = [i for i in dat.columns]
    columns = [str(i) for i in dat.index] + ['allPvalue', 'allScore']
    results = pd.DataFrame(1, columns=columns, index=pertGenes)
    ### 统计在每个cluster的卡方差异
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
    dat = pd.read_csv('/home/wzt/project/HC-CrisprScreen/results/dataInfo.tsv', sep='\t', header=0)
    for dirName, species in zip(dat["data"], dat["species"]):
        if dirNames and dirName not in dirNames: continue
        os.chdir(dirName)
        if os.path.isfile('mixscape_hvg_filter.h5ad') and os.path.isfile('clusterDistribution_beforeMix.tsv'):
            try:
                ClusterDistribution('clusterDistribution_beforeMix.tsv', 'clusterDistribution_beforeMix_pvalue.tsv')
                ClusterDistribution('clusterDistribution_afterMix.tsv', 'clusterDistribution_afterMix_pvalue.tsv')
            except Exception as e:
                print (e); print (dirName)


def PlotClusterDistribution():
    perturbationClass = pd.read_csv('perturbationClass.tsv', sep='\t')
    strongP = list(perturbationClass['perturbation'][:25])
    if 'CTRL' not in strongP: strongP.append('CTRL')
    count = pd.read_csv('clusterDistribution_afterMix.tsv', sep='\t', index_col=0)
    count = count[strongP]
    if len(count.index.unique()) >= 50:  ### 最多保留50个cluster
        a = [int(i) for i in count.index.unique()]
        a = sorted(a)[:50]
        count = count.loc[a, :]
    count1 = count / count.sum(axis=0)
    dat = pd.read_csv('clusterDistribution_afterMix_pvalue.tsv', sep='\t', index_col=0).T
    dat = dat.loc[[str(i) for i in count1.index] + ['allPvalue', 'allScore', 'allPvalue_adjust'], count1.columns]
#### enrichmap图
    cmap_custom = mcolors.LinearSegmentedColormap.from_list("custom_cmap", ["white", "red"])
    fig, ax1 = plt.subplots(figsize=(12, 8))
    ax = sns.heatmap(count1, vmin=0, vmax=1, cmap = cmap_custom, linecolor='black', linewidths=1, annot=False, center=True)
    ax.set_ylabel('Cluster')
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


#### 处理文件，让余老师产生动态图片
def PlotClusterDistribution1():
    perturbationClass = pd.read_csv('perturbationClass.tsv', sep='\t')
    strongP = list(perturbationClass['perturbation'][:25])
    if 'CTRL' not in strongP: strongP.append('CTRL')
    count = pd.read_csv('clusterDistribution_afterMix.tsv', sep='\t', index_col=0)
    count = count[strongP]
    if len(count.index.unique()) >= 50:  ### 最多保留50个cluster
        a = [int(i) for i in count.index.unique()]
        a = sorted(a)[:50]
        count = count.loc[a, :]
    count1 = count / count.sum(axis=0)
    dat = pd.read_csv('clusterDistribution_afterMix_pvalue.tsv', sep='\t', index_col=0).T
    dat = dat.loc[[str(i) for i in count1.index], count1.columns]
    dat.index.name='leiden'
    dat.to_csv('interact_plot_data/clusterDistribution_pvalue_afterMix.csv', sep=',', index=True)

def f_PlotClusterDistribution1():
    dat = pd.read_csv('/home/wzt/project/HC-CrisprScreen/results/dataInfo.tsv', sep='\t', header=0)
    for dirName, species in zip(dat["data"], dat["species"]):
        if dirNames and dirName not in dirNames: continue
        os.chdir(dirName)
        if os.path.isfile('mixscape_hvg_filter.h5ad') and os.path.isfile('clusterDistribution_afterMix_pvalue.tsv'):
            try:
                PlotClusterDistribution1()
            except Exception as e:
                print (e); print (dirName)



### 使用hvg表达谱计算扰动之间的相关性
def plotCor1():
    plt.rcParams.update({'font.size': 15})   ### 改变字体大小
    dat = pd.read_csv('/home/wzt/project/HC-CrisprScreen/results/dataInfo.tsv', sep='\t', header=0)
    for dirName, species in zip(dat["data"], dat["species"]):
        if dirNames and dirName not in dirNames: continue
        os.chdir(dirName)
        if not os.path.isfile('mixscape_hvg_filter.h5ad'): print (dirName); continue
        if os.path.isfile('figures/perturbation_cor_exp.png'): pass
        try:
            perturbationClass = pd.read_csv('perturbationClass.tsv', sep='\t')
            strongP = list(perturbationClass['perturbation'][:25])
            if 'CTRL' not in strongP: strongP.append('CTRL')
            adata = sc.read_h5ad('mixscape_hvg_filter_subset.h5ad')
            dat = pd.DataFrame(adata.obsm['X_pca'])
            df = dat.groupby(list(adata.obs['gene']), axis=0).mean()
            strongP = [i for i in strongP if i in list(df.index)]
            df1 = df.loc[strongP, :]
            cor_dat1 = calCosine(df1, df1)
            cor_dat = calCosine(df, df)
            if not os.path.isdir('corBetweenPerturb'): os.makedirs('corBetweenPerturb')
            cor_dat.to_csv('corBetweenPerturb/corExp.tsv', sep='\t', header=True, index=True)
            cor_dat1.to_csv('corBetweenPerturb/corExp_sub.tsv', sep='\t', header=True, index=True)

            perturbationClass = pd.read_csv('perturbationClass.tsv', sep='\t')
            strong = list(perturbationClass[perturbationClass['perturbationClass'] == 'strongPerturbation']['perturbation'])
            weak = list(perturbationClass[perturbationClass['perturbationClass'] == 'weakPerturbation']['perturbation'])
            row_colors = ['r' if i in strong else 'g' for i in cor_dat1.index]
            
            g = sns.clustermap(cor_dat1, cmap="vlag", vmin=-1, vmax=1, robust=False, figsize=(15,15), row_colors = row_colors, cbar_kws={"ticks": [-1, 0, 1]}, linecolor='black')   ###
            reordered_rows = g.dendrogram_row.reordered_ind
            reordered_cols = g.dendrogram_col.reordered_ind
            cor_dat2 = cor_dat1.iloc[reordered_rows, reordered_cols]
            cor_dat2.to_csv('corBetweenPerturb/corExp_order.tsv', sep='\t', header=True, index=True)
            with open('corBetweenPerturb/corExp_row.linkage.pkl', 'wb') as fout:
                pickle.dump(obj=g.dendrogram_row.linkage, file=fout)
            with open('corBetweenPerturb/corExp_col.linkage.pkl', 'wb') as fout:
                pickle.dump(obj=g.dendrogram_col.linkage, file=fout)     
            
            
            # legend_labels = ['strong', 'weak']; legend_colors = ['r', 'g']
            # legend_elements = [plt.Line2D([0], [0], marker='o', color='black', markerfacecolor=color, markersize=15, label=label,linewidth=0,) for label, color in zip(legend_labels, legend_colors)]
            # g.ax_heatmap.legend(legend_elements, legend_labels, loc='upper left', title='PerturbationClass', bbox_to_anchor=(-0.2, 1.2), frameon=False)
            # g.ax_heatmap.set_title("Correlation Heatmap of Perturbations", pad = 200)
            # g.savefig('figures/perturbation_cor_exp.png', dpi=300, bbox_inches='tight', transparent=False)
        except Exception as e:
            print (e); print (dirName)


def plotCor3():
    plt.rcParams.update({'font.size': 15})   ### 改变字体大小
    dat = pd.read_csv('/home/wzt/project/HC-CrisprScreen/results/dataInfo.tsv', sep='\t', header=0)
    for dirName, species in zip(dat["data"], dat["species"]):
        if dirNames and dirName not in dirNames: continue
        os.chdir(dirName)
        if not os.path.isfile('GSFA/fit.rds'): print (dirName); continue
        if os.path.isfile('figures/perturbation_cor.png'): pass
        try:
            perturbationClass = pd.read_csv('perturbationClass.tsv', sep='\t')
            perturbationClass = perturbationClass[[True if ',' not in i else False for i in perturbationClass['perturbation']]]  ### 没有组合扰动
            strongP = list(perturbationClass['perturbation'][:25])
            if 'CTRL' not in strongP: strongP.append('CTRL')
            dat = pd.read_csv('GSFA/beta_pm.tsv', sep='\t', index_col=0)
            tmp = [i for i in dat.columns if i != 'offset']
            dat = dat.loc[:, tmp]

            strongP = [i for i in strongP if i in list(dat.columns)]
            dat1 = dat[strongP]
            
            dat += 0.000001
            dat1 += 0.000001
            
            cor_dat = calCosine(dat.T, dat.T)
            cor_dat1 = calCosine(dat1.T, dat1.T)

            cor_dat.to_csv('corBetweenPerturb/corGSFA.tsv', sep='\t', header=True, index=True)
            cor_dat1.to_csv('corBetweenPerturb/corGSFA_sub.tsv', sep='\t', header=True, index=True)
            
            perturbationClass = pd.read_csv('perturbationClass.tsv', sep='\t')
            strong = list(perturbationClass[perturbationClass['perturbationClass'] == 'strongPerturbation']['perturbation'])
            weak = list(perturbationClass[perturbationClass['perturbationClass'] == 'weakPerturbation']['perturbation'])
            row_colors = ['r' if i in strong else 'g' for i in cor_dat1.index]
            
            g = sns.clustermap(cor_dat1, cmap="vlag", vmin=-1, vmax=1, robust=False, figsize=(15,15), row_colors = row_colors, cbar_kws={"ticks": [-1, 0, 1]}, linecolor='black')   ### 画edist距离，聚类热图
            reordered_rows = g.dendrogram_row.reordered_ind
            reordered_cols = g.dendrogram_col.reordered_ind
            cor_dat2 = cor_dat1.iloc[reordered_rows, reordered_cols]
            cor_dat2.to_csv('corBetweenPerturb/corGSFA_order.tsv', sep='\t', header=True, index=True)
            with open('corBetweenPerturb/corGSFA_row.linkage.pkl', 'wb') as fout:
                pickle.dump(obj=g.dendrogram_row.linkage, file=fout)
            with open('corBetweenPerturb/corGSFA_col.linkage.pkl', 'wb') as fout:
                pickle.dump(obj=g.dendrogram_col.linkage, file=fout)
            
            # legend_labels = ['strong', 'weak']; legend_colors = ['r', 'g']
            # legend_elements = [plt.Line2D([0], [0], marker='o', color='black', markerfacecolor=color, markersize=15, label=label,linewidth=0,) for label, color in zip(legend_labels, legend_colors)]
            # g.ax_heatmap.legend(legend_elements, legend_labels, loc='upper left', title='PerturbationClass', bbox_to_anchor=(-0.2, 1.2), frameon=False)
            # g.ax_heatmap.set_title("Correlation Heatmap of Perturbations", pad = 200)
            # g.savefig('figures/perturbation_cor.png', dpi=300, bbox_inches='tight', transparent=False)
        except Exception as e:
            print (e); print (dirName)



def perturbHeatMap_fun(pert):
    sc.set_figure_params(dpi=600)
    try:
        if pert != 'CTRL':
            pert1 = pert.replace('/', '')
            #if os.path.isfile('figures/heatmapfigures/{}.png'.format(pert1)):
            pt.pl.ms.heatmap(
            adata=adata,
            labels="gene",
            target_gene=pert,
            layer="X_pert",
            control="CTRL", save='figures/{}.png'.format(pert1), swap_axes=True, figsize=(10, 7))
    except:
        print ('*********wrong*********')

def perturbHeatMap(dirName):
    def fun1(x):
        if x[-2:] == 'KO':
            return x[:-2] + 'SP'
        else:
            return x
    os.chdir(dirName)
    print (dirName)
    if not os.path.isdir('figures/heatmapfigures'): os.makedirs('figures/heatmapfigures')
    global adata
    if os.path.isfile('mixscape_hvg_subset.h5ad'):
        adata = sc.read_h5ad('mixscape_hvg_subset.h5ad')
    else:
        adata = sc.read_h5ad('mixscape_hvg.h5ad')
    adata.obs['mixscape_class'] = list(adata.obs['mixscape_class'].apply(fun1))
    perts = adata.obs['gene'].unique()[:]
    myPool(perturbHeatMap_fun, perts, processes=20)


def f_perturbHeatMap():
    dat = pd.read_csv('/home/wzt/project/HC-CrisprScreen/results/dataInfo.tsv', sep='\t', header=0)
    for dirName, species in zip(dat["data"], dat["species"]):
        if dirNames and dirName not in dirNames: continue
        os.chdir(dirName)
        if os.path.isfile('mixscape_hvg.h5ad'):
            perturbHeatMap(dirName)




def MixscapeRatio():
    dat = pd.read_csv('/home/wzt/project/HC-CrisprScreen/results/dataInfo.tsv', sep='\t', header=0)
    for dirName, species in zip(dat["data"], dat["species"]):
        if dirNames and dirName not in dirNames: continue
        os.chdir(dirName)
        print (dirName)
        if not os.path.isfile('mixscape_hvg.h5ad'): continue
        perturbationClass = pd.read_csv('perturbationClass.tsv', sep='\t', index_col=0)
        perturbationClass.columns = ['distance', 'perturbationClass', 'SP', 'NP', 'EscapeRatio']
        perturbationClass.to_csv('perturbationClass.tsv', sep='\t', header=True, index=True)


        perturbationClass = pd.read_csv('perturbationClass.tsv', sep='\t')
        perturbationClass = perturbationClass.iloc[:25]   ### 选择基因

        fig, ax = plt.subplots(figsize = (10, 8))
        perturbationClass['SP'] = perturbationClass['SP'] / (perturbationClass['SP'] + perturbationClass['NP'])
        perturbationClass['NP'] = 1 - perturbationClass['SP']
        perturbationClass.sort_values('SP', ascending=False, inplace=True)
        ax.bar(perturbationClass['perturbation'], perturbationClass['SP'], label='SP')
        ax.bar(perturbationClass['perturbation'], perturbationClass['NP'], bottom= perturbationClass['SP'], label='NP', color='gray')
        ax.set_ylabel('Percentage of cells')
        ax.set_title('Mixscape Escape Ratio')
        plt.xticks(rotation=90)
        ax.legend()
        plt.show()
        fig.savefig('figures/Mixscape_escapeRatio.png',  dpi=300, bbox_inches='tight', transparent=False)

'''
/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb-seq/PRJNA679579/TP53
'''

### 跑7个demo
dirNames  = [#'/home/wzt/project/HC-CrisprScreen/poolSC_data/8ECCITE-seq/PRJNA641353/ECCITE', 
             #'/home/wzt/project/HC-CrisprScreen/poolSC_data/16TAP-seq/PRJNA559094/SCREEN_chr8', 
             #'/home/wzt/project/HC-CrisprScreen/poolSC_data/16TAP-seq/PRJNA559094/SCREEN_chr11', 
             #'/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb-seq/PRJNA587707/sciPlex3_A549', 
             #'/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb-seq/PRJNA587707/sciPlex3_K562', 
             #'/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb-seq/PRJNA628589/GSE149215',    ### failed
             #'/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb-seq/PRJNA831566/K562_essential'
             ]

### 大于200的数据集
# dirNames  = [#'/home/wzt/project/HC-CrisprScreen/poolSC_data/16TAP-seq/PRJNA559094/SCREEN_chr8', 
#              #'/home/wzt/project/HC-CrisprScreen/poolSC_data/16TAP-seq/PRJNA559094/SCREEN_chr11', 
#              #'/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb-seq/PRJNA551220/PRJNA551220',
#              #'/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb-seq/PRJNA587707/sciPlex3_A549', 
#              #'/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb-seq/PRJNA587707/sciPlex3_K562', 
#              #'/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb-seq/PRJNA587707/sciPlex3_MCF7',
#              #'/home/wzt/project/HC-CrisprScreen/poolSC_data/8ECCITE-seq/Perturb-CITE-seq/Perturb-CITE-seq', 
#              #'/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb-seq/PRJNA831566/K562_essential',
#              #'/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb-seq/PRJNA831566/RPE1_essential', 
#              #'/home/wzt/project/HC-CrisprScreen/poolSC_data/17SHARE-seq/PRJNA893678/210322_TFAtlas', 
#              #'/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb-seq/PRJNA831566/K562_GW',
#              ]


dirNames = [
'/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb-seq/PRJNA625951/PRJMA625951',
'/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb-seq/E-MTAB-10698/E-MTAB-10698',
'/home/wzt/project/HC-CrisprScreen/poolSC_data/13arrayBased/MIX-Seq/expt1',
'/home/wzt/project/HC-CrisprScreen/poolSC_data/13arrayBased/MIX-Seq/expt10',
'/home/wzt/project/HC-CrisprScreen/poolSC_data/13arrayBased/MIX-Seq/expt3',
'/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb-seq/PRJNA628589/GSE134836',
'/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb-seq/PRJNA628589/GSE134838',
'/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb-seq/PRJNA628589/GSE134839',
'/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb-seq/PRJNA628589/GSE134841',
'/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb-seq/PRJNA628589/GSE149214',
'/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb-seq/PRJNA628589/GSE160244',
'/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb-seq/PRJNA488861/PRJNA488861',
'/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb-seq/PRJNA505909/PRJNA505909',
'/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb-seq/PRJNA587707/sciPlex2',
'/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb-seq/PRJNA587707/sciPlex3_A549',
'/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb-seq/PRJNA587707/sciPlex3_K562',
'/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb-seq/PRJNA587707/sciPlex3_MCF7',
'/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb-seq/PRJNA587707/sciPlex4_A549',
'/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb-seq/PRJNA587707/sciPlex4_MCF7'
]



'''
/home/wzt/anaconda3/bin/python   myPlot.py
'''

if __name__ == '__main__':
    print ('hello, world')
    #f_plotQC()     
    #f_plotSgRNAPlot()    ### 第一部分图

    #MixscapeRatio()      ### 第二部分图
    f_perturbHeatMap()
    #f_plotMixscape()

    #f_doPhenotype()     ####   第三部分图
    #f_ClusterDistribution()  ### 与ctrl的差异的卡方统计检验
    #f_PlotClusterDistribution()  ### 与ctrl的差异的enrichmentmap图
    #f_PlotClusterDistribution1()  ### 与ctrl的差异的enrichmentmap图

    
    #plotCor1()
    #plotCor3()
