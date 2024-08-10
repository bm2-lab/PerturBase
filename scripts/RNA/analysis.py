from myUtil import *
import scanpy as sc, numpy as np
import pandas as pd
import anndata as ad
import warnings
warnings.filterwarnings('ignore')
import multiprocessing
from scipy import stats
import seaborn as sns
from statsmodels.stats import multitest
import pickle
try:
    from upsetplot import from_contents
    from upsetplot import UpSet
    from matplotlib import pyplot as plt
except:
    pass
from collections import defaultdict
try:
    import pertpy as pt
except:
    pass
import os
from tqdm import tqdm


pd.set_option('display.float_format', lambda x: '%.4f' % x)
import seaborn as sns

import matplotlib
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 20})
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
# plt.rc('font',family='Times New Roman')
import sys
import os

def doEdist(dirName):
    os.chdir(dirName)
    if os.path.isfile('mixscape_hvg_filter.h5ad') and os.path.isfile('perturbationClass.tsv'):
        adata = sc.read_h5ad('mixscape_hvg_filter.h5ad') 
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

            df1 = distance.pairwise(adata, groupby="gene", verbose=True) 
            df1 = df1 / df1.max().max()
            df1.to_csv('corBetweenPerturb/corEdistance.tsv', sep='\t', header=True, index=True, index_label='')
        else:
            df = pd.read_csv('corBetweenPerturb/corEdistance_sub.tsv', sep='\t', header=0, index_col=0)    

def ClusterDistribution(filein, fileout):
    adata1 = sc.read_h5ad('mixscape_hvg.h5ad')
    adata2 = sc.read_h5ad('mixscape_hvg_filter.h5ad')
    count1 = pd.crosstab(index=adata1.obs["leiden"], columns=adata1.obs["gene"])
    count2 = pd.crosstab(index=adata2.obs["leiden"], columns=adata2.obs["gene"])
    count1.to_csv('clusterDistribution_beforeMix.tsv', sep='\t', index=True)
    count2.to_csv('clusterDistribution_afterMix.tsv', sep='\t', index=True)

    dat = pd.read_csv(filein, sep='\t', index_col=0)
    pertGenes = [i for i in dat.columns]
    columns = [str(i) for i in dat.index] + ['allPvalue', 'allScore']
    results = pd.DataFrame(1, columns=columns, index=pertGenes)
    for gene in pertGenes:
        for cluster in dat.index:
            sumGene = dat[gene].sum()   ### pertb total nums
            sumCTRL = dat['CTRL'].sum()  ###  ctrl total nums
            sumCluster = dat.loc[cluster, :].sum()  ### cluster nums
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
    results.to_csv(fileout, sep='\t', header=True, index=True)
    return results



###  prepare barcode file for identifing differentially expressed genes
def f_preData1():
    try:
        adata = sc.read_h5ad('mixscape_hvg_filter.h5ad')
        strongP = getStrongPerturb(Nums=200)
        adata = adata[adata.obs['gene'].isin(strongP)]   ### only keep 200 perturbations
        def fun1(adata):
            barcode = adata.obs[['gene']]
            mylist1  =[]; mylist2 = []
            for cell, gene in zip(barcode.index, barcode['gene']):
                genes = gene.split(',')
                [mylist1.append(cell) for _ in genes]
                [mylist2.append(i) for i in genes]
            return mylist1, mylist2
        
        ### scMAGeCK file
        mylist1, mylist2 = fun1(adata)
        barcode = pd.DataFrame({'cell':mylist1, 'barcode':mylist2, "gene": mylist2}, index=mylist1)
        barcode.to_csv('barcode.txt', sep='\t', index=False, header=True)
        
        ###SCEPTRE barcode file
        barcode = adata.obs[['gene']]
        tmp = [True if ',' not in x else False for x in barcode['gene']]
        barcode = adata.obs[['gene']][tmp]
        
        a1 = barcode[barcode['gene'] != 'CTRL']
        a2 = barcode[barcode['gene'] == 'CTRL']
        from itertools import cycle
        mycycle = cycle(['CTRL1', 'CTRL2'])
        a2['gene'] =  [i for i, j in zip (mycycle, range(a2.shape[0]))]
        barcode = pd.concat([a2, a1], axis=0)

        barcode1 = pd.DataFrame(0, columns=barcode.index, index=barcode['gene'].unique())
        for cell, gene in zip(barcode.index, barcode['gene']):
            barcode1.loc[gene, cell] = 100
        barcode1.to_csv('barcode_sceptre.txt', sep='\t', index=True, header=True)

        tmp = ['non-targeting' if x.startswith('CTRL') else x for x in barcode['gene'].unique()]
        grna_group_data_frame = pd.DataFrame({'grna_id': barcode['gene'].unique(), 'grna_target': tmp})
        grna_group_data_frame.to_csv('grnaGroup.tsv', sep='\t', header=True, index=False)

        ### GSFA barcode file
        barcode = adata.obs[['gene']]
        geneUnique = list(set([x for i in barcode['gene'] for x in i.split(',')]))
        barcode1 = pd.DataFrame(0, columns=barcode.index, index= geneUnique)
        for cell, genes in zip(barcode.index, barcode['gene']):
            for gene in genes.split(','):
                barcode1.loc[gene, cell] = 100
        barcode1.to_csv('barcode_GSFA.txt', sep='\t', index=True, header=True)
    except Exception as e:
        print (e)


def runMagecK_rra(species,R_env,Script_base_path,Python_env,msigdb_Signature_path):
    try:
        if not os.path.isdir('scMageCK'): os.makedirs('scMageCK')
        path_to_python = f"{Python_env}/bin/python"
        myUtil_r = f'{Script_base_path}/scripts/RNA/myUtil.r'

        if species == 'Homo sapiens':
            cmd = '{}/bin/Rscript   {}/scripts/RNA/doScMAGeCK2.r    hsa {} {} {} {}'.format(R_env,Script_base_path,path_to_python,myUtil_r,msigdb_Signature_path,Script_base_path)
        else:
            cmd = '{}/bin/Rscript    {}/scripts/RNA/doScMAGeCK2.r    mmu {} {} {} {} '.format(R_env,Script_base_path,path_to_python,myUtil_r,msigdb_Signature_path,Script_base_path)
        subprocess.call(cmd, shell=True)
        files = glob.glob('scMageCK/GENE_SET/hallmark*')
        hallmarkName = [os.path.basename(i)[:-8] for i in files]
        dat = pd.read_table(files[0], index_col=0)
        with open('scMageCK/RRA.txt', 'w') as fout:
            fout.write('GeneSet\tPerturbation\tSelection\t-log10(FDR)\n')
            for file, hallmark in zip(files, hallmarkName):
                dat = pd.read_table(file, index_col=0)
                for gene in dat.index:
                    if dat.loc[gene, 'FDR.low'] <= 0.01:
                        fout.write('{}\t{}\tNegative\t{:.3f}\n'.format(hallmark, gene, -np.log10(dat.loc[gene, 'FDR.low'])))
                    elif dat.loc[gene, 'FDR.high'] <= 0.01:
                        fout.write('{}\t{}\tPositive\t{:.3f}\n'.format(hallmark, gene, -np.log10(dat.loc[gene, 'FDR.high'])))
        dat = pd.read_csv('scMageCK/RRA.txt', sep='\t')
        if dat.shape[0] == 0:  return
        dat.sort_values(['Perturbation', 'GeneSet'], inplace=True)
        fileout = 'scMageCK/RRA.txt'
        dat.to_csv(fileout, sep='\t', header=True, index=False)
    except Exception as e:
        print (e)


def runMagecK_lr(species,R_env,Script_base_path,Python_env):
    try:
        if not os.path.isdir('scMageCK'): os.makedirs('scMageCK')
        if not os.path.isfile('mixscape_hvg_filter.h5ad'): return
        path_to_python = f"{Python_env}/bin/python"
        myUtil_r = f'{Script_base_path}/scripts/RNA/myUtil.r'
        if species == 'Homo sapiens':
            cmd = '{}/bin/Rscript   {}/scripts/RNA/doScMAGeCK3.r    hsa {} {}'.format(R_env,Script_base_path,path_to_python,myUtil_r)
        else:
            cmd = '{}/bin/Rscript   {}/scripts/RNA/doScMAGeCK3.r    mmu {} {}'.format(R_env,Script_base_path,path_to_python,myUtil_r)
        subprocess.call(cmd, shell=True)
    except Exception as e:
        print(f"failed: {e}")


### second deg tool
def runsceptre(species,R_env,Script_base_path,Python_env):
    try:
        if not os.path.isfile('mixscape_hvg_filter.h5ad'): return
        if not os.path.isdir('sceptre'): os.makedirs('sceptre')
        path_to_python = f"{Python_env}/bin/python"
        myUtil_r = f'{Script_base_path}/scripts/RNA/myUtil.r'
        if species == 'Homo sapiens':
            cmd = '{}/bin/Rscript   {}/scripts/RNA/doSceptre.r    hsa {} {} '.format(R_env,Script_base_path,path_to_python,myUtil_r)
        else:
            cmd = '{}/bin/Rscript    {}/scripts/RNA/doSceptre.r    mmu {} {}'.format(R_env,Script_base_path,path_to_python,myUtil_r)
        subprocess.call(cmd, shell=True)
    except Exception as e:
        print(f"failed: {e}");


def runScanpy_ttest():
    try:
        if not os.path.isfile('mixscape_hvg_filter.h5ad'): return
        if not os.path.isdir('ttest'): os.makedirs('ttest')
        adata= sc.read_h5ad('mixscape_hvg_filter.h5ad')
        try:
            adata.uns['log1p']["base"] = None
        except:
            pass
        genes = [i for i in set(adata.obs['gene']) if i != 'CTRL']

        sc.tl.rank_genes_groups(adata, 'gene', groups=genes[:], reference='CTRL')
        result = adata.uns['rank_genes_groups']
        perturbations = result['names'].dtype.names
        genes = sorted(result['names'][perturbations[0]])
        final_result1 = pd.DataFrame(index=genes, columns=perturbations)
        for perturbation in perturbations:
            tmp1 = result['names'][perturbation]
            tmp2 = result['pvals_adj'][perturbation]
            sorted_lists = sorted(zip(tmp1, tmp2), key=lambda x: x[0])
            sorted_list1, sorted_list2 = zip(*sorted_lists)
            final_result1[perturbation] = sorted_list2
            final_result1.fillna(1, inplace=True)
        final_result1.to_csv('ttest/pvalue.tsv', sep='\t', header=True, index=True)

        final_result2 = pd.DataFrame(index=genes, columns=perturbations)
        for perturbation in perturbations:
            tmp1 = result['names'][perturbation]
            tmp2 = result['logfoldchanges'][perturbation]
            sorted_lists = sorted(zip(tmp1, tmp2), key=lambda x: x[0])
            sorted_list1, sorted_list2 = zip(*sorted_lists)
            sorted_list2 = [2 ** i for i in sorted_list2]
            final_result2[perturbation] = sorted_list2
            final_result2.fillna(1, inplace=True)
        final_result2.to_csv('ttest/foldchange.tsv', sep='\t', header=True, index=True)
    except Exception as e:
        print (e)



from sklearn.metrics.pairwise import cosine_similarity
def calCosine(Xtr, Xte):
    dat_cor = pd.DataFrame(cosine_similarity(Xte,Xtr))  ###row Xte, column Xtr
    dat_cor.columns = Xtr.index
    dat_cor.index = Xte.index
    return dat_cor


#### GSFA 
def runGSFA(species,R_env,Script_base_path,Python_env):
    try:
        if not os.path.isdir('GSFA'): os.makedirs('GSFA')
        if not os.path.isfile('mixscape_hvg_filter.h5ad'): return
        path_to_python = f"{Python_env}/bin/python"
        myUtil_r = f'{Script_base_path}/scripts/RNA/myUtil.r'
        if species == 'Homo sapiens':
            cmd = '{}/bin/Rscript   {}/scripts/RNA/doGSFA.r    hsa {} {} '.format(R_env,Script_base_path,path_to_python,myUtil_r)
        else:
            cmd = '{}/bin/Rscript    {}/scripts/RNA/doGSFA.r    mmu {} {}'.format(R_env,Script_base_path,path_to_python,myUtil_r)
        subprocess.call(cmd, shell=True)
    except Exception as e:
        print(f"failed: {e}")

def runScanpy_wilcoxon():
    try:
        if not os.path.isfile('mixscape_hvg_filter.h5ad'): return
        if not os.path.isdir('Wilcoxon'): os.makedirs('Wilcoxon')
        adata= sc.read_h5ad('mixscape_hvg_filter.h5ad')
        try:
            adata.uns['log1p']["base"] = None
        except:
            pass
        genes = [i for i in set(adata.obs['gene']) if i != 'CTRL']

        sc.tl.rank_genes_groups(adata, 'gene', groups=genes[:], reference='CTRL', method='wilcoxon')
        result = adata.uns['rank_genes_groups']
        perturbations = result['names'].dtype.names
        genes = sorted(result['names'][perturbations[0]])
        final_result1 = pd.DataFrame(index=genes, columns=perturbations)
        for perturbation in perturbations:
            tmp1 = result['names'][perturbation]
            tmp2 = result['pvals_adj'][perturbation]
            sorted_lists = sorted(zip(tmp1, tmp2), key=lambda x: x[0])
            sorted_list1, sorted_list2 = zip(*sorted_lists)
            final_result1[perturbation] = sorted_list2
            final_result1.fillna(1, inplace=True)
        final_result1.to_csv('Wilcoxon/pvalue.tsv', sep='\t', header=True, index=True)

        final_result2 = pd.DataFrame(index=genes, columns=perturbations)
        for perturbation in perturbations:
            tmp1 = result['names'][perturbation]
            tmp2 = result['logfoldchanges'][perturbation]
            sorted_lists = sorted(zip(tmp1, tmp2), key=lambda x: x[0])
            sorted_list1, sorted_list2 = zip(*sorted_lists)
            sorted_list2 = [2 ** i for i in sorted_list2]
            final_result2[perturbation] = sorted_list2
            final_result2.fillna(1, inplace=True)
        final_result2.to_csv('Wilcoxon/foldchange.tsv', sep='\t', header=True, index=True)
    except Exception as e:
        print (e)

def getDEG():
    ### wilcoxon
    dat1 = pd.read_csv('Wilcoxon/pvalue.tsv', sep='\t', index_col=0)
    dat1['gene'] = dat1.index
    melted_df1 = pd.melt(dat1, id_vars=['gene'], var_name='perturbation', value_name='Pvalue')

    dat2 = pd.read_csv('Wilcoxon/foldchange.tsv', sep='\t', index_col=0)
    dat2['gene'] = dat2.index
    melted_df2 = pd.melt(dat2, id_vars=['gene'], var_name='perturbation', value_name='foldchange')
    
    dat = pd.merge(left=melted_df1, right=melted_df2)
    dat = dat[['perturbation', 'gene', 'Pvalue', 'foldchange']]
    dat.to_csv('Wilcoxon/deg_result.tsv', sep='\t', header=True, index=False)


    ### ttest
    dat1 = pd.read_csv('ttest/pvalue.tsv', sep='\t', index_col=0)
    dat1['gene'] = dat1.index
    melted_df1 = pd.melt(dat1, id_vars=['gene'], var_name='perturbation', value_name='Pvalue')

    dat2 = pd.read_csv('ttest/foldchange.tsv', sep='\t', index_col=0)
    dat2['gene'] = dat2.index
    melted_df2 = pd.melt(dat2, id_vars=['gene'], var_name='perturbation', value_name='foldchange')
    
    dat = pd.merge(left=melted_df1, right=melted_df2)
    dat = dat[['perturbation', 'gene', 'Pvalue', 'foldchange']]
    dat.to_csv('ttest/deg_result.tsv', sep='\t', header=True, index=False)

### scMageck
    dat1 = pd.read_csv('scMageCK/scMageckLR_score_pval.txt', sep='\t', index_col=0)
    dat1['gene'] = dat1.index
    melted_df1 = pd.melt(dat1, id_vars=['gene'], var_name='perturbation', value_name='Pvalue')

    dat2 = pd.read_csv('scMageCK/scMageckLR_score.txt', sep='\t', index_col=0)
    dat2['gene'] = dat2.index
    melted_df2 = pd.melt(dat2, id_vars=['gene'], var_name='perturbation', value_name='scMageCK_LRscore')
    
    dat = pd.merge(left=melted_df1, right=melted_df2)
    dat = dat[['perturbation', 'gene', 'Pvalue', 'scMageCK_LRscore']]
    dat.to_csv('scMageCK/deg_result.tsv', sep='\t', header=True, index=False)

#### sceptre
    def tmp_fun(x):
        if x != 1 : return 2** x
        else: return 1
    dat1 = pd.read_csv('sceptre/rawResult.tsv', sep='\t', index_col=0)
    dat1['gene'] = dat1.index
    dat1 = dat1[['grna_target', 'gene', 'p_value', 'log_2_fold_change']]
    dat1.columns = ['perturbation', 'gene', 'Pvalue', 'foldchange']
    dat1.fillna(1, inplace=True)
    dat1['foldchange'] = dat1['foldchange'].apply(tmp_fun)
    dat1.to_csv('sceptre/deg_result.tsv', sep='\t', header=True, index=False)

### GSFA
    dat1 = pd.read_csv('GSFA/lfsr.tsv', sep='\t', index_col=0)
    dat1 = dat1.iloc[:, :-1]
    dat1['gene'] = dat1.index
    melted_df1 = pd.melt(dat1, id_vars=['gene'], var_name='perturbation', value_name='Pvalue')

    dat2 = pd.read_csv('GSFA/effect.tsv', sep='\t', index_col=0)
    dat2['gene'] = dat2.index
    melted_df2 = pd.melt(dat2, id_vars=['gene'], var_name='perturbation', value_name='GSFA_effect_score')
    
    dat = pd.merge(left=melted_df1, right=melted_df2)
    dat = dat[['perturbation', 'gene', 'Pvalue', 'GSFA_effect_score']]
    dat.to_csv('GSFA/deg_result.tsv', sep='\t', header=True, index=False)
def filterDEG(dirName):
    os.chdir(dirName)
    foldchange = 2; pvalue = 0.01; score_cut = 0.2

    ### wilcoxon   
    mydict = {}
    a = pd.read_csv('Wilcoxon/pvalue.tsv', sep='\t', index_col=0)
    b = pd.read_csv('Wilcoxon/foldchange.tsv', sep='\t', index_col=0)
    for pertGene in b.columns:
        tmp = (a[pertGene] <= pvalue) & ((b[pertGene] <= 1/foldchange) | (b[pertGene] >= foldchange))
        deg = list(tmp[tmp].index)
        mydict[pertGene] = deg
    with open('Wilcoxon/deg.pkl', 'wb') as fout:
        pickle.dump(obj=mydict, file=fout)


    ### ttest
    mydict = {}
    a = pd.read_csv('ttest/pvalue.tsv', sep='\t', index_col=0)
    b = pd.read_csv('ttest/foldchange.tsv', sep='\t', index_col=0)
    for pertGene in b.columns:
        tmp = (a[pertGene] <= pvalue) & ((b[pertGene] <= 1/foldchange) | (b[pertGene] >= foldchange))
        deg = list(tmp[tmp].index)
        mydict[pertGene] = deg
    with open('ttest/deg.pkl', 'wb') as fout:
        pickle.dump(obj=mydict, file=fout)

### scMageck
    mydict = {}
    a = pd.read_csv('scMageCK/scMageckLR_score_pval.txt', sep='\t', index_col=0)
    b = pd.read_csv('scMageCK/scMageckLR_score.txt', sep='\t', index_col=0)
    for pertGene in b.columns:
        tmp = (a[pertGene] <= pvalue) & ((b[pertGene] <= -score_cut) | (b[pertGene] >= score_cut))
        deg = list(tmp[tmp].index)
        mydict[pertGene] = deg
    with open('scMageCK/deg.pkl', 'wb') as fout:
        pickle.dump(obj=mydict, file=fout)

#### sceptre
    mydict = defaultdict(list)
    a = pd.read_csv('sceptre/rawResult.tsv', sep='\t', index_col=0)
    b = a[(a['p_value'] <= pvalue) & ((a["log_2_fold_change"] <= -np.log2(foldchange)) | (a["log_2_fold_change"] >= np.log2(foldchange)))]
    for i, j in zip(b['grna_target'], b.index):
        mydict[i].append(j)
    with open('sceptre/deg.pkl', 'wb') as fout:
        pickle.dump(obj=mydict, file=fout)

### GSFA
    mydict = defaultdict(list)
    dat = pd.read_csv('GSFA/lfsr.tsv', sep='\t', index_col=0)
    dat = dat.iloc[:, :-1]  ### offset
    for pertGene in dat.columns:
        tmp = dat[pertGene] <= pvalue
        deg = list(tmp[tmp].index)
        mydict[pertGene] = deg
    with open('GSFA/deg.pkl', 'wb') as fout:
        pickle.dump(obj=mydict, file=fout)

def plotCor1(dirName):
    os.chdir(dirName)
    try:
        perturbationClass = pd.read_csv('perturbationClass.tsv', sep='\t')
        strongP = list(perturbationClass['perturbation'][:25])
        if 'CTRL' not in strongP: strongP.append('CTRL')
        adata = sc.read_h5ad('mixscape_hvg_filter.h5ad')
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


def plotCor3(dirName):
    os.chdir(dirName)
    if not os.path.isfile('GSFA/fit.rds'): print (dirName); return
    if os.path.isfile('figures/perturbation_cor.png'): pass
    try:
        perturbationClass = pd.read_csv('perturbationClass.tsv', sep='\t')
        perturbationClass = perturbationClass[[True if ',' not in i else False for i in perturbationClass['perturbation']]]  
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



if __name__ == '__main__':
    path = sys.argv[1]
    os.chdir(path)
    species = sys.argv[2]
    R_env = sys.argv[3]
    Script_base_path = sys.argv[4]
    Python_env = sys.argv[5]
    msigdb_Signature_path = sys.argv[6]
    if species == 'Hs':
        species = 'Homo sapiens'
    else:
        species = 'Mus musculus'
    os.makedirs('corBetweenPerturb', exist_ok=True)
    print(species)
    ##  ClusterDistribution test
    print('ClusterDistribution')
    ClusterDistribution('clusterDistribution_beforeMix.tsv', 'clusterDistribution_beforeMix_pvalue.tsv')
    ClusterDistribution('clusterDistribution_afterMix.tsv',  'clusterDistribution_afterMix_pvalue.tsv')

    ### calculate edistance between perturbation
    print('ClusterDistribution')
    doEdist(path)
    print('prepare barcode')
    f_preData1()  ### prepare barcode file for deg methods

    #### deg detection
    print('DEG detection')
    print('runMagecK_rra')
    runMagecK_rra(species,R_env,Script_base_path,Python_env,msigdb_Signature_path)
    print('runMagecK_lr')
    runMagecK_lr(species,R_env,Script_base_path,Python_env) ### first deg method
    print('runsceptre')
    runsceptre(species,R_env,Script_base_path,Python_env)  ### second deg method
    print('runttest')
    runScanpy_ttest()  ### third and fourth deg method
    print('runwilcox')
    runScanpy_wilcoxon()
    print('runGSFA')
    runGSFA(species,R_env,Script_base_path,Python_env)
    
    # other correlation
    plotCor1(path)
    plotCor3(path)
    
    

    ### process deg file
    print('process deg file')
    getDEG()
