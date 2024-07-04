from myUtil import *
import scanpy as sc
import warnings
import anndata as ad
warnings.filterwarnings('ignore')
from scipy import sparse
from scipy import stats
import signal
import multiprocessing
import pickle
try:
    from upsetplot import from_contents
    from upsetplot import UpSet
    from matplotlib import pyplot as plt
except:
    pass

from tqdm import tqdm


pd.set_option('display.float_format', lambda x: '%.4f' % x)
import seaborn as sns

import matplotlib.pyplot as plt
import matplotlib
plt.rcParams.update({'font.size': 20})   ### 改变字体大小
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
plt.rc('font',family='Times New Roman')
plt.rcParams.update({'font.size': 10})   ### 改变字体大小



###  准备差异表达软件的barcode等
def f_preData1(X):
    try:
        dirName, _ = X
        os.chdir(dirName)
        if os.path.isfile('barcode_GSFA.txt') and os.path.isfile('barcode.txt') and os.path.isfile('mixscape_hvg_filter_subset.h5ad'): return
        if not os.path.isfile('mixscape_hvg_filter.h5ad'): 
            print ('mixscape_hvg_filter.h5ad  not exist! **{}**\n'.format(dirName))
            return
        if os.path.isfile('mixscape_hvg_filter_subset.h5ad'):
            adata = sc.read_h5ad('mixscape_hvg_filter_subset.h5ad')
        else:
            adata = sc.read_h5ad('mixscape_hvg_filter.h5ad')
            strongP = getStrongPerturb(Nums=200)
            adata = adata[adata.obs['gene'].isin(strongP)]
            adata.write_h5ad('mixscape_hvg_filter_subset.h5ad')   #### 只保留前200个扰动, 对于全基因扰动，CTRL保留10000个细胞
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
        
        ### 准备sceptre的barcode, 去除组合扰动
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
        grna_group_data_frame = pd.DataFrame({'grna_id': barcode['gene'].unique(), 'grna_group': tmp})
        grna_group_data_frame.to_csv('grnaGroup.tsv', sep='\t', header=True, index=False)

        ### 准备GSFA的barcode, 保留组合扰动
        barcode = adata.obs[['gene']]
        geneUnique = list(set([x for i in barcode['gene'] for x in i.split(',')]))
        barcode1 = pd.DataFrame(0, columns=barcode.index, index= geneUnique)
        for cell, genes in zip(barcode.index, barcode['gene']):
            for gene in genes.split(','):
                barcode1.loc[gene, cell] = 100
        barcode1.to_csv('barcode_GSFA.txt', sep='\t', index=True, header=True)
    except Exception as e:
        print (e); print (dirName)


def ff_preData1():
    mylist = []
    dat = pd.read_csv('/home/wzt/project/HC-CrisprScreen/results/dataInfo.tsv', sep='\t', header=0)
    for dirName, species in zip(dat["data"], dat["species"]):
        if dirNames and dirName not in dirNames: continue
        mylist.append([dirName, species])
    myPool(f_preData1, mylist, processes=2)




def runMagecK_rra(X):
    try:
        dirName, species = X
        os.chdir(dirName)
        if not os.path.isfile('mixscape_hvg_filter_subset.h5ad'): return
        if not os.path.isdir('scMageCK'): os.makedirs('scMageCK')
        if os.path.isfile('scMageCK/RRA.txt'): return     ### 第二次跑，避免重复
        if species == 'Homo sapiens':
            cmd = 'Rscript   /home/wzt/project/HC-CrisprScreen/doScMAGeCK2.r   {}  hsa'.format(dirName)
        else:
            cmd = 'Rscript   /home/wzt/project/HC-CrisprScreen/doScMAGeCK2.r   {}  mmu'.format(dirName)
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
        if dat.shape[0] == 0:  return ### 如果都没有差异或者没有文件
        dat.sort_values(['Perturbation', 'GeneSet'], inplace=True)
        fileout = 'scMageCK/RRA.txt'
        dat.to_csv(fileout, sep='\t', header=True, index=False)
    except Exception as e:
        print (e)



def f_runMagecK_rra():
    mylist = []
    dat = pd.read_csv('/home/wzt/project/HC-CrisprScreen/results/dataInfo.tsv', sep='\t', header=0)
    for dirName, species in zip(dat["data"], dat["species"]):
        if dirNames and dirName not in dirNames: continue
        mylist.append([dirName, species])
    myPool(runMagecK_rra, mylist, processes=3)




def runMagecK_lr(X):
    try:
        dirName, species = X
        os.chdir(dirName)
        if not os.path.isdir('scMageCK'): os.makedirs('scMageCK')
        if not os.path.isfile('mixscape_hvg_filter.h5ad'): return
        if os.path.isfile('scMageCK/deg.tsv'): return   ### 第二次跑，避免重复
        if species == 'Homo sapiens':
            cmd = 'Rscript   /home/wzt/project/HC-CrisprScreen/doScMAGeCK3.r   {}  hsa'.format(dirName)
        else:
            cmd = 'Rscript   /home/wzt/project/HC-CrisprScreen/doScMAGeCK3.r   {}  mmu'.format(dirName)
        subprocess.call(cmd, shell=True)
    except Exception as e:
        print(f"命令执行失败: {e}"); print (dirName)


def f_runMagecK_lr():
    mylist = []
    dat = pd.read_csv('/home/wzt/project/HC-CrisprScreen/results/dataInfo.tsv', sep='\t', header=0)
    for dirName, species in zip(dat["data"], dat["species"]):
        if dirNames and dirName not in dirNames: continue
        mylist.append([dirName, species])
    pool = multiprocessing.Pool(processes=3)
    pool.map(runMagecK_lr, mylist)
    pool.close()
    pool.join()


### 第二个差异基因软件
def runsceptre(X):
    try:
        dirName, species = X
        os.chdir(dirName)
        if not os.path.isfile('mixscape_hvg_filter.h5ad'): return
        if not os.path.isdir('sceptre'): os.makedirs('sceptre')
        if os.path.isfile('sceptre/rawResult.tsv'): return   ### 第二次跑，避免重复
        if species == 'Homo sapiens':
            cmd = 'Rscript   /home/wzt/project/HC-CrisprScreen/doSceptre.r   {}  hsa'.format(dirName)
        else:
            cmd = 'Rscript   /home/wzt/project/HC-CrisprScreen/doSceptre.r   {}  mmu'.format(dirName)
        subprocess.call(cmd, shell=True)
    except Exception as e:
        print(f"命令执行失败: {e}"); print (dirName)

def f_runsceptre():
    mylist = []
    dat = pd.read_csv('/home/wzt/project/HC-CrisprScreen/results/dataInfo.tsv', sep='\t', header=0)
    for dirName, species in zip(dat["data"], dat["species"]):
        if dirNames and dirName not in dirNames: continue
        mylist.append([dirName, species])
    pool = multiprocessing.Pool(processes=2)
    pool.map(runsceptre, mylist)
    pool.close()
    pool.join()


#### scanpy 内置的差异基因软件
def runScanpy(X):
    try:
        dirName, _ = X
        os.chdir(dirName)
        if not os.path.isfile('mixscape_hvg_filter.h5ad'): return
        if not os.path.isdir('ttest'): os.makedirs('ttest')
        if os.path.isfile('ttest/foldchange.tsv'): return
        adata= sc.read_h5ad('mixscape_hvg_filter_subset.h5ad')
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
        print (e); print (dirName)

### 第三个差异基因软件
def f_runScanpy():
    mylist = []
    dat = pd.read_csv('/home/wzt/project/HC-CrisprScreen/results/dataInfo.tsv', sep='\t', header=0)
    for dirName, species in zip(dat["data"], dat["species"]):
        if dirNames and dirName not in dirNames: continue
        mylist.append([dirName, species])
    myPool(runScanpy, mylist, processes=8)

from sklearn.metrics.pairwise import cosine_similarity
def calCosine(Xtr, Xte):
    dat_cor = pd.DataFrame(cosine_similarity(Xte,Xtr))  ###行是Xte, 列是Xtr
    dat_cor.columns = Xtr.index
    dat_cor.index = Xte.index
    return dat_cor


#### GSFA  差异表达基因 并画图
def runGSFA(X):
    dirName, species = X
    os.chdir(dirName)
    log_print = open('runGSFA.log', 'w'); sys.stdout = log_print; sys.stderr = log_print
    try:
        if not os.path.isdir('GSFA'): os.makedirs('GSFA')
        if not os.path.isfile('mixscape_hvg_filter.h5ad'): return
        if os.path.isfile('GSFA/fit.rds'): return   ### 第二次跑，避免重复
        if species == 'Homo sapiens':
            cmd = 'Rscript   /home/wzt/project/HC-CrisprScreen/doGSFA.r   {}  hsa'.format(dirName)
        else:
            cmd = 'Rscript   /home/wzt/project/HC-CrisprScreen/doGSFA.r   {}  mmu'.format(dirName)
        subprocess.call(cmd, shell=True)
    except Exception as e:
        print(f"命令执行失败: {e}"); print (dirName)
    log_print.close()


### 第四个差异基因软件
def f_runGSFA():
    mylist = []
    dat = pd.read_csv('/home/wzt/project/HC-CrisprScreen/results/dataInfo.tsv', sep='\t', header=0)
    for dirName, species in zip(dat["data"], dat["species"]):
        if dirNames and dirName not in dirNames: continue
        mylist.append([dirName, species])
    pool = multiprocessing.Pool(processes=3)
    pool.map(runGSFA, mylist)
    pool.close()
    pool.join()


#### scanpy 内置的差异基因软件
def runScanpy1(X):
    try:
        dirName, _ = X
        os.chdir(dirName)
        if not os.path.isfile('mixscape_hvg_filter_subset.h5ad'): return
        if not os.path.isdir('Wilcoxon'): os.makedirs('Wilcoxon')
        if os.path.isfile('Wilcoxon/foldchange.tsv'): return
        adata= sc.read_h5ad('mixscape_hvg_filter_subset.h5ad')
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
        print (e); print (dirName)

### 第五个差异基因软件
def f_runScanpy1():
    mylist = []
    dat = pd.read_csv('/home/wzt/project/HC-CrisprScreen/results/dataInfo.tsv', sep='\t', header=0)
    for dirName, species in zip(dat["data"], dat["species"]):
        if dirNames and dirName not in dirNames: continue
        mylist.append([dirName, species])
    myPool(runScanpy1, mylist, processes=8)


def getDEG(dirName):
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

# ### scMageck
#     mydict = {}
#     a = pd.read_csv('scMageCK/scMageckLR_score_pval.txt', sep='\t', index_col=0)
#     b = pd.read_csv('scMageCK/scMageckLR_score.txt', sep='\t', index_col=0)
#     for pertGene in b.columns:
#         tmp = (a[pertGene] <= pvalue) & ((b[pertGene] <= -score_cut) | (b[pertGene] >= score_cut))
#         deg = list(tmp[tmp].index)
#         mydict[pertGene] = deg
#     with open('scMageCK/deg.pkl', 'wb') as fout:
#         pickle.dump(obj=mydict, file=fout)

# #### sceptre
#     mydict = defaultdict(list)
#     a = pd.read_csv('sceptre/rawResult.tsv', sep='\t', index_col=0)
#     b = a[(a['p_value'] <= pvalue) & ((a["log_2_fold_change"] <= -np.log2(foldchange)) | (a["log_2_fold_change"] >= np.log2(foldchange)))]
#     for i, j in zip(b['grna_group'], b.index):
#         mydict[i].append(j)
#     with open('sceptre/deg.pkl', 'wb') as fout:
#         pickle.dump(obj=mydict, file=fout)

# ### GSFA
#     mydict = defaultdict(list)
#     dat = pd.read_csv('GSFA/lfsr.tsv', sep='\t', index_col=0)
#     dat = dat.iloc[:, :-1]  ### 最后一列为offset
#     for pertGene in dat.columns:
#         tmp = dat[pertGene] <= pvalue
#         deg = list(tmp[tmp].index)
#         mydict[pertGene] = deg
#     with open('GSFA/deg.pkl', 'wb') as fout:
#         pickle.dump(obj=mydict, file=fout)


def f_getDEG():
    dat = pd.read_csv('/home/wzt/project/HC-CrisprScreen/results/dataInfo.tsv', sep='\t', header=0)
    for dirName in tqdm(dat["data"]):
        if dirNames and dirName not in dirNames: continue
        os.chdir(dirName)
        if os.path.isfile('mixscape_hvg_filter.h5ad') and os.path.isfile('perturbationClass.tsv'):
            try:
                getDEG(dirName)
            except Exception as e:
                print (e); print (dirName)


#### 统计每个数据集的差异基因列表个数，以及画upset图片
def degPlot(geneName):
    def fun1(software, geneName):
        filein = '{}/deg.pkl'.format(software)
        if os.path.isfile(filein):
            with open(filein, 'rb') as fin:
                tmp_dict = pickle.load(fin)
                return tmp_dict.get(geneName, [])
        else:
            return []

    #if os.path.isfile('degSet/{}.png'.format(geneName)):return

    Wilcoxon = fun1('Wilcoxon', geneName)
    ttest = fun1('ttest', geneName)
    sceptre = fun1('sceptre', geneName)
    scMageCK = fun1('scMageCK', geneName)
    GSFA = fun1('GSFA', geneName)
    if  len(Wilcoxon + ttest + sceptre + scMageCK + GSFA) == 0: 
        with plt.style.context('seaborn-white'):
            fig, ax = plt.subplots(figsize=(10, 8))
            plt.text(0.3, 0.5, 'NO DEGs in all of the methods', fontsize=22)
            plt.savefig('degSet/{}.png'.format(geneName),  dpi=300, bbox_inches='tight', transparent=False)
            plt.show()
    else:
        countTable = from_contents({'Wilcoxon': Wilcoxon, 't-test': ttest, 'SCEPTRE': sceptre, 'scMAGeCK':scMageCK, 'GSFA':GSFA})
        with plt.style.context('seaborn-white'):
            fig, ax = plt.subplots(figsize=(10, 8))
            plt.xticks([])
            plt.yticks([])
            ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False); ax.spines['left'].set_visible(False)
            UpSet(countTable, subset_size='count', with_lines=True).plot(fig)
            #plt.suptitle('The intersection of DEGs')
            ax.set_xlabel('The intersection of DEGs')
            geneName = geneName.replace('/', '')
            plt.savefig('degSet/{}.png'.format(geneName),  dpi=300, bbox_inches='tight', transparent=False)
            plt.show()


def f_degPlot():
    dat = pd.read_csv('/home/wzt/project/HC-CrisprScreen/results/dataInfo.tsv', sep='\t', header=0)
    for dirName in dat["data"]:
        if dirNames and dirName not in dirNames: continue
        os.chdir(dirName)
        print (dirName)
        if os.path.isfile('mixscape_hvg_filter.h5ad') and os.path.isfile('perturbationClass.tsv'):
            dat = pd.read_csv('perturbationClass.tsv', sep='\t')
            geneNames = getStrongPerturb(Nums=2000000)
            if not os.path.isdir("degSet"): os.makedirs("degSet")
            for geneName in tqdm(geneNames):
                try:
                    degPlot(geneName)
                except Exception as e:
                    print (e); print (dirName)


### 差异基因的barplot
def degBarplot():
    dat = pd.read_csv('perturbationClass.tsv', sep='\t')
    geneNames = getStrongPerturb(Nums=25)   ### 一张图只展示前25个基因
    def fun1(software, geneNames):
        filein = '{}/deg.pkl'.format(software)
        if os.path.isfile(filein):
            with open(filein, 'rb') as fin:
                tmp =  pickle.load(fin)
                return [len(tmp.get(i, [])) for i in geneNames]
        else:
            return [0] * len(geneNames)
    mylen1 = fun1('Wilcoxon', geneNames)
    dat1 =  pd.DataFrame({'software':'Wilcoxon', 'Counts': mylen1, 'Perturbation':geneNames})
    
    mylen2 = fun1('ttest', geneNames)
    dat2 =  pd.DataFrame({'software':'t-test', 'Counts': mylen2, 'Perturbation':geneNames})
    
    mylen3 = fun1('sceptre', geneNames)
    dat3 =  pd.DataFrame({'software':'SCEPTRE', 'Counts': mylen3, 'Perturbation':geneNames})
    
    mylen4 = fun1('scMageCK', geneNames)
    dat4 =  pd.DataFrame({'software':'scMAGeCK', 'Counts': mylen4, 'Perturbation':geneNames})
    
    mylen5 = fun1('GSFA', geneNames)
    dat5 =  pd.DataFrame({'software':'GSFA', 'Counts': mylen5, 'Perturbation':geneNames})

    dat = pd.concat([dat1, dat2, dat3, dat4, dat5])
    fig, ax = plt.subplots(figsize=(14, 8))
    sns.barplot(dat, x='Perturbation', y='Counts', hue='software', ax=ax)
    for p in ax.patches:
        ax.annotate(f'{p.get_height():.0f}', (p.get_x() + p.get_width() / 2., p.get_height()), ha='center', va='center', fontsize=10, color='black', xytext=(0, 5), textcoords='offset points')
    plt.xticks(rotation=90)
    # 显示图形
    plt.show()
    fig.savefig('figures/degBarplot.png', transparent=False, dpi=300, bbox_inches='tight')

def f_degBarplot():
    dat = pd.read_csv('/home/wzt/project/HC-CrisprScreen/results/dataInfo.tsv', sep='\t', header=0)
    for dirName in tqdm(dat["data"]):
        if dirNames and dirName not in dirNames: continue
        os.chdir(dirName)
        if os.path.isfile('mixscape_hvg_filter.h5ad') and os.path.isfile('perturbationClass.tsv'):
            if not os.path.isdir('figures'):  os.makedirs('figures')
            try:
                degBarplot()
            except Exception as e:
                print (e); print (dirName)

#### 把每个软件的差异基因制作成表格方便展示
def getDEG1(dirName):
    os.chdir(dirName)

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

# ### scMageck
#     dat1 = pd.read_csv('scMageCK/scMageckLR_score_pval.txt', sep='\t', index_col=0)
#     dat1['gene'] = dat1.index
#     melted_df1 = pd.melt(dat1, id_vars=['gene'], var_name='perturbation', value_name='Pvalue')

#     dat2 = pd.read_csv('scMageCK/scMageckLR_score.txt', sep='\t', index_col=0)
#     dat2['gene'] = dat2.index
#     melted_df2 = pd.melt(dat2, id_vars=['gene'], var_name='perturbation', value_name='scMageCK_LRscore')
    
#     dat = pd.merge(left=melted_df1, right=melted_df2)
#     dat = dat[['perturbation', 'gene', 'Pvalue', 'scMageCK_LRscore']]
#     dat.to_csv('scMageCK/deg_result.tsv', sep='\t', header=True, index=False)

# #### sceptre
#     def tmp_fun(x):
#         if x != 1 : return 2** x
#         else: return 1
#     dat1 = pd.read_csv('sceptre/rawResult.tsv', sep='\t', index_col=0)
#     dat1['gene'] = dat1.index
#     dat1 = dat1[['grna_group', 'gene', 'p_value', 'log_2_fold_change']]
#     dat1.columns = ['perturbation', 'gene', 'Pvalue', 'foldchange']
#     dat1.fillna(1, inplace=True)
#     dat1['foldchange'] = dat1['foldchange'].apply(tmp_fun)
#     dat1.to_csv('sceptre/deg_result.tsv', sep='\t', header=True, index=False)

# ### GSFA
#     dat1 = pd.read_csv('GSFA/lfsr.tsv', sep='\t', index_col=0)
#     dat1 = dat1.iloc[:, :-1]
#     dat1['gene'] = dat1.index
#     melted_df1 = pd.melt(dat1, id_vars=['gene'], var_name='perturbation', value_name='Pvalue')

#     dat2 = pd.read_csv('GSFA/effect.tsv', sep='\t', index_col=0)
#     dat2['gene'] = dat2.index
#     melted_df2 = pd.melt(dat2, id_vars=['gene'], var_name='perturbation', value_name='GSFA_effect_score')
    
#     dat = pd.merge(left=melted_df1, right=melted_df2)
#     dat = dat[['perturbation', 'gene', 'Pvalue', 'GSFA_effect_score']]
#     dat.to_csv('GSFA/deg_result.tsv', sep='\t', header=True, index=False)

def f_getDEG1():
    dat = pd.read_csv('/home/wzt/project/HC-CrisprScreen/results/dataInfo.tsv', sep='\t', header=0)
    for dirName in tqdm(dat["data"]):
        if dirNames and dirName not in dirNames: continue
        os.chdir(dirName)
        if os.path.isfile('mixscape_hvg_filter.h5ad') and os.path.isdir('GSFA'):
            try:
                getDEG1(dirName)
            except Exception as e:
                print (e); print (dirName)


def volcanoPlot(dirName):
    ### 五个软件
    os.chdir(dirName)
    print (dirName)
    ### ttest
    try:
        if not os.path.isdir('figures/volcano/ttest'): os.makedirs('figures/volcano/ttest')
        dat1 = pd.read_csv('ttest/deg_result.tsv', sep='\t')
        strongP = getStrongPerturb(Nums=20000)
        for gene in tqdm(strongP):
            fileout = 'figures/volcano/ttest/{}.png'.format(gene.replace('/', ''))
            if os.path.isfile(fileout): continue                                                     
            if  gene == 'CTRL': continue
            dat = dat1[dat1['perturbation'] == gene]
            with open('ttest/deg.pkl', 'rb') as fin:
                mydict = pickle.load(fin)[gene]
            dat.loc[dat['foldchange'] == 1, 'Pvalue'] = 1
            dat['-log10(Pvalue)']  = -np.log10(dat['Pvalue'])
            dat['-log10(Pvalue)'] = dat['-log10(Pvalue)'].apply(lambda x: x if x <= 10 else 10)
            dat['-log2(foldchange)']  = np.log2(dat['foldchange'])
            dat['-log2(foldchange)']  =  dat['-log2(foldchange)'].apply(lambda x: x if abs(x) <= 5 else 5)
            dat['deg'] = dat['gene'].apply(lambda x: x in mydict)
            color_palette = {True: 'red', False: 'grey'}
            fig, ax = plt.subplots(figsize=(8, 8))
            sns.scatterplot(data=dat, x = '-log2(foldchange)', y = '-log10(Pvalue)', hue='deg', legend=False, marker='+', palette = color_palette)
            plt.axhline(y=2, color='black', linestyle='--')
            plt.axvline(x=1, color='black', linestyle='--')
            plt.axvline(x=-1, color='black', linestyle='--')
            fig.savefig(fileout, transparent=False, dpi=300, bbox_inches='tight')

        # wilcoxon
        if not os.path.isdir('figures/volcano/Wilcoxon'): os.makedirs('figures/volcano/Wilcoxon')
        dat1 = pd.read_csv('Wilcoxon/deg_result.tsv', sep='\t')
        strongP = getStrongPerturb(Nums=20000)
        for gene in tqdm(strongP):
            fileout = 'figures/volcano/Wilcoxon/{}.png'.format(gene.replace('/', ''))
            if os.path.isfile(fileout): continue  
            if  gene == 'CTRL': continue
            dat = dat1[dat1['perturbation'] == gene]
            with open('Wilcoxon/deg.pkl', 'rb') as fin:
                mydict = pickle.load(fin)[gene]
            dat.loc[dat['foldchange'] == 1, 'Pvalue'] = 1
            dat['-log10(Pvalue)']  = -np.log10(dat['Pvalue'])
            dat['-log10(Pvalue)'] = dat['-log10(Pvalue)'].apply(lambda x: x if x <= 10 else 10)
            dat['-log2(foldchange)']  = np.log2(dat['foldchange'])
            dat['-log2(foldchange)']  =  dat['-log2(foldchange)'].apply(lambda x: x if abs(x) <= 5 else 5)
            dat['deg'] = dat['gene'].apply(lambda x: x in mydict)
            color_palette = {True: 'red', False: 'grey'}
            fig, ax = plt.subplots(figsize=(8, 8))
            sns.scatterplot(data=dat, x = '-log2(foldchange)', y = '-log10(Pvalue)', hue='deg', legend=False, marker='+', palette = color_palette)
            plt.axhline(y=2, color='black', linestyle='--')
            plt.axvline(x=1, color='black', linestyle='--')
            plt.axvline(x=-1, color='black', linestyle='--')
            fig.savefig(fileout, transparent=False, dpi=300, bbox_inches='tight')


    #     ### sceptre
    #     if not os.path.isdir('figures/volcano/sceptre'): os.makedirs('figures/volcano/sceptre')
    #     dat1 = pd.read_csv('sceptre/deg_result.tsv', sep='\t')
    #     strongP = getStrongPerturb(Nums=200)
    #     for gene in strongP:
    #         fileout = 'figures/volcano/sceptre/{}.png'.format(gene.replace('/', ''))
    #         if os.path.isfile(fileout): continue 
    #         if  gene == 'CTRL': continue
    #         dat = dat1[dat1['perturbation'] == gene]
    #         with open('sceptre/deg.pkl', 'rb') as fin:
    #             mydict = pickle.load(fin)
    #         #if gene not in mydict: continue   ### 只会考虑单个扰动
    #         if ',' in gene: continue
    #         mydict = mydict[gene]   ####
    #         dat.loc[dat['foldchange'] == 1, 'Pvalue'] = 1
    #         dat['-log10(Pvalue)']  = -np.log10(dat['Pvalue'])
    #         dat['-log10(Pvalue)'] = dat['-log10(Pvalue)'].apply(lambda x: x if x <= 10 else 10)
    #         dat['-log2(foldchange)']  = np.log2(dat['foldchange'])
    #         dat['-log2(foldchange)']  =  dat['-log2(foldchange)'].apply(lambda x: x if abs(x) <= 5 else 5)
    #         dat['deg'] = dat['gene'].apply(lambda x: x in mydict)
    #         color_palette = {True: 'red', False: 'grey'}
    #         fig, ax = plt.subplots(figsize=(8, 8))
    #         sns.scatterplot(data=dat, x = '-log2(foldchange)', y = '-log10(Pvalue)', hue='deg', legend=False, marker='+', palette = color_palette)
    #         plt.axhline(y=2, color='black', linestyle='--')
    #         plt.axvline(x=1, color='black', linestyle='--')
    #         plt.axvline(x=-1, color='black', linestyle='--')
    #         fig.savefig(fileout, transparent=False, dpi=300, bbox_inches='tight')

    # ## scMageck
    #     if not os.path.isdir('figures/volcano/scMageck'): os.makedirs('figures/volcano/scMageck')
    #     dat1 = pd.read_csv('scMageCK/deg_result.tsv', sep='\t')
    #     dat1.columns = ['perturbation', 'gene', 'Pvalue', 'scMAGeCK_LRscore']
    #     strongP = getStrongPerturb(Nums=200)
    #     for gene in strongP:
    #         fileout = 'figures/volcano/scMageck/{}.png'.format(gene.replace('/', ''))
    #         if os.path.isfile(fileout): continue
    #         if  gene == 'CTRL': continue
    #         dat = dat1[dat1['perturbation'] == gene]
    #         with open('scMageCK/deg.pkl', 'rb') as fin:
    #             mydict = pickle.load(fin)
    #         #if gene not in mydict: continue   ### scMageck只会考虑单个扰动
    #         if ',' in gene: continue
    #         mydict = mydict[gene]   ####
    #         dat['-log10(Pvalue)']  = -np.log10(dat['Pvalue'])
    #         dat['-log10(Pvalue)'] = dat['-log10(Pvalue)'].apply(lambda x: x if x <= 10 else 10)
    #         dat['deg'] = dat['gene'].apply(lambda x: x in mydict)
    #         color_palette = {True: 'red', False: 'grey'}
    #         fig, ax = plt.subplots(figsize=(8, 8))
    #         sns.scatterplot(data=dat, x = 'scMAGeCK_LRscore', y = '-log10(Pvalue)', hue='deg', legend=False, marker='+', palette = color_palette)
    #         plt.axhline(y=2, color='black', linestyle='--')
    #         plt.axvline(x=0.2, color='black', linestyle='--')
    #         plt.axvline(x=-0.2, color='black', linestyle='--')
    #         fileout = 'figures/volcano/scMageck/{}.png'.format(gene.replace('/', ''))
    #         fig.savefig(fileout, transparent=False, dpi=300, bbox_inches='tight')

    #     ## GSFA
    #     if not os.path.isdir('figures/volcano/GSFA'): os.makedirs('figures/volcano/GSFA')
    #     dat1 = pd.read_csv('GSFA/deg_result.tsv', sep='\t')
    #     strongP = getStrongPerturb(Nums=200)
    #     for gene in strongP:
    #         fileout = 'figures/volcano/GSFA/{}.png'.format(gene.replace('/', ''))
    #         if os.path.isfile(fileout): continue 
    #         if  gene == 'CTRL': continue
    #         dat = dat1[dat1['perturbation'] == gene]
    #         with open('GSFA/deg.pkl', 'rb') as fin:
    #             mydict = pickle.load(fin)
    #         #if gene not in mydict: continue   ### GSFA 不会输出组合扰动
    #         if ',' in gene: continue
    #         mydict = mydict[gene]   ####
    #         dat['-log10(Pvalue)']  = -np.log10(dat['Pvalue'])
    #         dat['-log10(Pvalue)'] = dat['-log10(Pvalue)'].apply(lambda x: x if x <= 10 else 10)
    #         dat['deg'] = dat['gene'].apply(lambda x: x in mydict)
    #         color_palette = {True: 'red', False: 'grey'}
    #         fig, ax = plt.subplots(figsize=(8, 8))
    #         sns.scatterplot(data=dat, x = 'GSFA_effect_score', y = '-log10(Pvalue)', hue='deg', legend=False, marker='+', palette = color_palette)
    #         plt.axhline(y=2, color='black', linestyle='--')
    #         fig.savefig(fileout, transparent=False, dpi=300, bbox_inches='tight')
    except Exception as e:
        print (e); print (dirName)


def f_volcanoPlot():
    dat = pd.read_csv('/home/wzt/project/HC-CrisprScreen/results/dataInfo.tsv', sep='\t', header=0)
    for dirName in tqdm(dat["data"]):
        if dirNames and dirName not in dirNames: continue
        os.chdir(dirName)
        if os.path.isfile('mixscape_hvg_filter.h5ad'):
            volcanoPlot(dirName)

    # mylist = []
    # for dirName in tqdm(dat["data"]):
    #     if dirNames and dirName not in dirNames: continue
    #     mylist.append(dirName)
    # myPool(volcanoPlot, mylist, processes=4)

### 查看哪些文件已经完成
if False:
    dat = pd.read_csv('/home/wzt/project/HC-CrisprScreen/results/dataInfo.tsv', sep='\t', header=0)
    for dirName, species in zip(dat["data"], dat["species"]):
        os.chdir(dirName)
        if os.path.isfile('mixscape_hvg_filter.h5ad'):
            #for fileName in ['scMageCK/RRA.txt','scMageCK/deg.tsv', 'sceptre/rawResult.tsv', 'GSFA/fit.rds', 'ttest/pvalue.tsv', 'Wilcoxon/foldchange.tsv']:
            for fileName in ['perturbationClass.tsv']:
                if not os.path.isfile(fileName):
                    print (fileName, dirName)

'''
/home/wzt/anaconda3/bin/python   calDEG.py
'''

#### 对RRA进行subsample
def subSample_RRA():
    dat = pd.read_csv('/home/wzt/project/HC-CrisprScreen/results/dataInfo.tsv', sep='\t', header=0)
    for dirName, species in zip(dat["data"], dat["species"]):
        if dirNames and dirName not in dirNames:
            continue
        else:
            os.chdir(dirName)
            perturbationClass = pd.read_csv('perturbationClass.tsv', sep='\t')
            strongP = list(perturbationClass['perturbation'][:25])
            tmp = [True if i in strongP else False for i in dat['Perturbation']]
            dat = pd.read_csv('scMageCK/RRA.txt', sep='\t')
            dat = dat[tmp]
            dat.to_csv('scMageCK/RRA_order.txt', index=False)


dirNames  = [#'/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb-seq/PRJNA831566/RPE1_essential', 
             #'/home/wzt/project/HC-CrisprScreen/poolSC_data/17SHARE-seq/PRJNA893678/210322_TFAtlas', 
             #'/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb-seq/PRJNA831566/K562_essential',
             '/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb-seq/PRJNA831566/K562_GW',
             ]


if __name__ == '__main__':
    print ('hello, world')
    #ff_preData1()  ### 产生barcode文件
    #f_runMagecK_rra()
   
    #f_runMagecK_lr()
    #f_runsceptre()
    #f_runScanpy()
    #f_runScanpy1()  ### wilcoxon
    #f_runGSFA()
    
    #f_getDEG1()    ### 产生deg_result.tsv文件，方便展示
    #f_getDEG()     ### 获得差异基因pkl文件
    
    ####f_degBarplot()  ### 差异基因柱状图
    f_volcanoPlot()   ### 火山图
    f_degPlot()  ### 五组维恩图
    