from myUtil import *
import scanpy as sc
import warnings
warnings.filterwarnings('ignore')

import pertpy as pt #type: ignore

pd.set_option('display.max_colwidth', 100)

def fun1(adata, isWrite = True):
    if isWrite: adata.write_h5ad('raw.h5ad')
    a = adata.shape
    adata = adata[~adata.obs['gene'].isin(['None', 'CTRL'])]

    tmp = [i.split(',') for i in adata.obs['gene']]
    tmp = [i for i in tmp if i != 'CTRL' and i != 'None']
    tmp = np.unique([i for j in tmp for i in j])
    return a, len(tmp)

def classGene(adata, threshold = 0.2):  #### 20%以上产生了性状，
    def fun(x):
        if x in strongP: return 'strong'
        elif x == 'CTRL': return 'CTRL'
        else: return 'weak'

    count = pd.crosstab(index=adata.obs["mixscape_class_global"], columns=adata.obs["gene"])
    all_cells_percentage = pd.melt(count / count.sum(), ignore_index=False).reset_index()
    KO_cells_percentage = all_cells_percentage[all_cells_percentage["mixscape_class_global"] == "KO"]
    KO_cells_percentage = KO_cells_percentage.sort_values("value", ascending=False)
    strongP = list(KO_cells_percentage[KO_cells_percentage['value'] >= threshold].gene)
    weakP = [i for i in list(adata.obs['gene'].unique()) if i not in strongP]
    weakP = [i for i in weakP if i != 'CTRL']
    adata.obs['pertclass'] = adata.obs['gene'].apply(fun)
    adata1 = adata[~((adata.obs['mixscape_class_global'] == 'NP') & (adata.obs['pertclass'] == 'strong'))]  #### 过滤掉扰动为strong，但是没有产生性状的细胞
    adata2 = adata[((adata.obs['mixscape_class_global'] == 'NP') & (adata.obs['pertclass'] == 'strong'))]   #### 扰动为strong，但是没有产生性状的细胞
    return strongP, weakP, adata1, adata2


def doMixScape1(adata, filterNone=True, minNums = 30, shuffle=True, filterCom = False, mtpercent=10, min_genes=200, keepall = False):  #### 得到hvg基因的local_perturbation
    if adata.shape[1] <= 1500: min_genes = 0
    if keepall: minNums = 1
    filterNoneNums, filterCells, filterMT, filterMinNums, adata = preData(adata, filterNone, minNums, shuffle, filterCom,  seed=42, mtpercent=mtpercent, min_genes=min_genes)
    if len(adata.obs['gene'].unique()) <=1:
        return filterNoneNums, filterCells, filterMT, filterMinNums, adata, 0
    if adata.shape[1] <= 1500:  ### 针对TAP-seq, 保留所有基因
        sc.pp.highly_variable_genes(adata, subset=True, n_top_genes=adata.shape[1])
    else:
        sc.pp.highly_variable_genes(adata, subset=True, n_top_genes=4000)
    sc.pp.pca(adata)
    if not isinstance(adata.X, np.ndarray):  ###
        adata.X = adata.X.toarray()
    adata = doPCA1(adata)
    adata.write_h5ad('RawNormal.h5ad')   #### 进行了基本的过滤 和log矫正
    allGeneAffilter = adata.var_names

    mixscape_identifier = pt.tl.Mixscape()
    mixscape_identifier.perturbation_signature(adata, 'gene', 'CTRL')
    if keepall:
        adata.obs["mixscape_class_global"] = 'KO'
    else:
        mixscape_identifier.mixscape(adata = adata, control = 'CTRL', labels='gene', layer='X_pert')

    if not isinstance(adata.layers['X_pert'], np.ndarray):
        adata.layers['X_pert'] = adata.layers['X_pert'].toarray()
    return filterNoneNums, filterCells, filterMT, filterMinNums, adata, allGeneAffilter


def doPCA(adata):    ###对过滤后的数据进行PCA等处理, 使用copy保存原始的X数据
    adata.layers['normal_counts'] = adata.X.copy()
    adata.X = adata.layers['X_pert']
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata)
    return adata

def doPCA1(adata):    ###对过滤后的数据进行PCA等处理, 使用copy保存原始的X数据
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata)
    return adata

def doMixScape2(adata, filterNone=True, minNums = 30, shuffle=True, filterCom = False, mtpercent=10, min_genes=200, keepall = False):  ###  #### 得到所有基因的local_perturbation, 后续方便计算
    if adata.shape[1] <= 1500: min_genes = 0
    if keepall: minNums = 1
    *_, adata = preData(adata, filterNone, minNums, shuffle, filterCom,  seed=42, mtpercent=mtpercent, min_genes=min_genes)
    sc.pp.highly_variable_genes(adata, subset=False)
    sc.pp.pca(adata)
    
    perturbation_signature(adata, 'gene', 'CTRL')
    return adata


def doTest1():
    dat = pd.read_csv('/home/wzt/project/HC-CrisprScreen/poolSC_data/info.tsv', sep='\t', header=None)
    for tech, PRJ, i in zip(dat[1], dat[6], dat[3]):
        dirName = '/home/wzt/project/HC-CrisprScreen/poolSC_data/{}/{}/{}'.format(tech, PRJ, i)
        os.chdir(dirName)
        if not os.path.isfile('mixscape.h5ad'):
            print (dirName)   

'''
没有ctrl 与 None
/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb-seq/PRJNA549582/TCell1_nonorm
/home/wzt/project/HC-CrisprScreen/poolSC_data/8ECCITE-seq/PRJNA556195/EGR2_KnockOut
/home/wzt/project/HC-CrisprScreen/poolSC_data/17SHARE-seq/PRJNA893678/180124_perturb
/home/wzt/project/HC-CrisprScreen/poolSC_data/17SHARE-seq/PRJNA893678/210715_combinatorial
/home/wzt/project/HC-CrisprScreen/poolSC_data/13arrayBased/CaltechDATA/GPCTP_2019   #### 药物处理



有CTRL但数据量太大, 后续处理   *******
/home/wzt/project/HC-CrisprScreen/poolSC_data/17SHARE-seq/PRJNA893678/210322_TFAtlas

enhancer
/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb-seq/PRJNA494734/pilot_lowmoi
/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb-seq/PRJNA494734/pilot_highmoi

/home/wzt/project/HC-CrisprScreen/poolSC_data/5Mosaic-seq/PRJNA322853/TAD2_SE2
/home/wzt/project/HC-CrisprScreen/poolSC_data/5Mosaic-seq/PRJNA322853/15SE_71HS

/home/wzt/project/HC-CrisprScreen/poolSC_data/16TAP-seq/PRJNA559094/SCREEN_chr8
/home/wzt/project/HC-CrisprScreen/poolSC_data/16TAP-seq/PRJNA559094/SCREEN_chr11
'''

def f_doMixScape(dirName, species):
    os.chdir(dirName)
    adata = sc.read_h5ad('raw.h5ad')
    rawTranscriptNums = adata.shape[1]
    adata = transID(adata, species)
    TranscriptNums = adata.shape[1]
    highMOI = ['pilot_lowmoi', 'pilot_highmoi', 'at_scale', 'TAD2_SE2', '15SE_71HS', 'K562-dCas9-KRAB_5K', 'SCREEN_chr8', 'SCREEN_chr11']
    keepall = False if os.path.basename(dirName) in highMOI else False
    filterNoneNums, filterCells, filterMT, filterMinNums, adata1, allGeneAffilter = doMixScape1(adata, mtpercent=mtpercent, keepall=keepall)
    if len(adata1.obs['gene'].unique()) <= 1:
        RawCellNums, RawGeneNums, CellNums, GeneNums = fun1(adata, False)[0][0], fun1(adata, False)[1], 0, 0
        return RawCellNums, RawGeneNums, CellNums, GeneNums, 0, 0, 0, 0, filterNoneNums, filterCells, filterMT, filterMinNums, rawTranscriptNums, TranscriptNums
    strongP, weakP, adata_pert, adata_pert1 = classGene(adata1, threshold=.2)  ### 对strong perturbation 把np过滤掉
    adata_filter = adata[adata_pert.obs_names, allGeneAffilter]   ### 没经过矫正，只是过滤的原始数据
    adata1.write_h5ad('mixscape_hvg.h5ad')

    #adata2 = doMixScape2(adata, mtpercent= mtpercent, keepall=keepall)
    #adata2.write_h5ad('mixscape.h5ad')
    #adata3 = adata2[adata_pert.obs_names, ]   #### 过滤掉无用的数据
    #adata3.write_h5ad('mixscape_filter.h5ad')

    adata_pert = doPCA(adata_pert)
    adata_pert.write_h5ad('mixscape_hvg_filter.h5ad')
    adata_filter.write_h5ad('filterRaw.h5ad')  ### 没经过矫正，只是过滤强扰动的未扰动信息，数据为count计数

    RawCellNums, RawGeneNums, CellNums, GeneNums = fun1(adata, False)[0][0], fun1(adata, False)[1], fun1(adata1, False)[0][0], fun1(adata1, False)[1]
    return RawCellNums, RawGeneNums, CellNums, GeneNums, len(strongP), len(weakP), adata_pert.shape[0], adata_pert1.shape[0], filterNoneNums, filterCells, filterMT, filterMinNums, rawTranscriptNums, TranscriptNums


def ff_doMixScape():
    fileout = '/home/wzt/project/HC-CrisprScreen/results/dataInfo.tsv'
    with open(fileout, 'w') as fout:
        fout.write('data\tspecies\tRawCellNums\tRawGeneNums\tCellNums\tGeneNums\tstrongP\tweakP\tCellNums_strongP\tCellNums_weakP\tfilterNoneNums\tfilterCells\tfilterMT\tfilterMinNums\trawTranscriptNums\tTranscriptNums\n')
        dirNameList = []; speciesList = []
        dat = pd.read_csv('/home/wzt/project/HC-CrisprScreen/datInfo.tsv', sep='\t', header=None)
        for tech, PRJ, i, species in zip(dat[1], dat[6], dat[3], dat[11]):
            dirName = '/home/wzt/project/HC-CrisprScreen/poolSC_data/{}/{}/{}'.format(tech, PRJ, i)
            dirNameList.append(dirName)
            speciesList.append(species)
        dat = pd.read_csv('/home/wzt/project/HC-CrisprScreen/datInfo1.tsv', sep='\t', header=None)
        dirNameList.extend(list(dat[0]))
        speciesList.extend(list(dat[1]))

        for dirName, species in zip(dirNameList, speciesList):
            if specific_dataset and dirName != specific_dataset: continue
            os.chdir(dirName)
            if dirName not in unuseList:
                print (dirName, species)
                result = f_doMixScape(dirName, species)
                result = [str(i) for i in result]
                fout.write('{}\t{}\t{}\n'.format(dirName, species, '\t'.join(result)))


### 追加模式，已经处理过的不需要处理了
def ff_doMixScape1():
    fileout = '/home/wzt/project/HC-CrisprScreen/results/dataInfo.tsv'
    with open(fileout, 'a+') as fout:
        complete_list = list(pd.read_csv(fileout, sep='\t')['data'])
        dirNameList = []; speciesList = []
        dat = pd.read_csv('/home/wzt/project/HC-CrisprScreen/datInfo.tsv', sep='\t', header=None)
        for tech, PRJ, i, species in zip(dat[1], dat[6], dat[3], dat[11]):
            dirName = '/home/wzt/project/HC-CrisprScreen/poolSC_data/{}/{}/{}'.format(tech, PRJ, i)
            dirNameList.append(dirName)
            speciesList.append(species)
        dat = pd.read_csv('/home/wzt/project/HC-CrisprScreen/datInfo1.tsv', sep='\t', header=None)
        dirNameList.extend(list(dat[0]))
        speciesList.extend(list(dat[1]))
        for dirName, species in zip(dirNameList, speciesList):
            if specific_dataset and dirName != specific_dataset: continue
            os.chdir(dirName)
            if dirName not in unuseList + complete_list:
                print (dirName, species)
                result = f_doMixScape(dirName, species)
                result = [str(i) for i in result]
                fout.write('{}\t{}\t{}\n'.format(dirName, species, '\t'.join(result)))


### 处理数据, 但是不输出日志文件
def ff_doMixScape2():
    dirNameList = []; speciesList = []
    dat = pd.read_csv('/home/wzt/project/HC-CrisprScreen/datInfo.tsv', sep='\t', header=None)
    for tech, PRJ, i, species in zip(dat[1], dat[6], dat[3], dat[11]):
        dirName = '/home/wzt/project/HC-CrisprScreen/poolSC_data/{}/{}/{}'.format(tech, PRJ, i)
        dirNameList.append(dirName)
        speciesList.append(species)
    dat = pd.read_csv('/home/wzt/project/HC-CrisprScreen/datInfo1.tsv', sep='\t', header=None)
    dirNameList.extend(list(dat[0]))
    speciesList.extend(list(dat[1]))
    for dirName, species in zip(dirNameList, speciesList):
        if specific_dataset and dirName != specific_dataset: continue
        os.chdir(dirName)
        if dirName not in unuseList:
            print (dirName, species)
            result = f_doMixScape(dirName, species)
            result = [str(i) for i in result]




### 追加模式，处理还未处理的数据集
def ff_doMixScape3():
    fileout = '/home/wzt/project/HC-CrisprScreen/results/dataInfo.tsv'
    with open(fileout, 'a+') as fout:
            result = f_doMixScape(specific_dataset, 'Homo sapiens')
            result = [str(i) for i in result]
            fout.write('{}\t{}\t{}\n'.format(specific_dataset, "Homo sapiens", '\t'.join(result)))


### 处理特定的数据集，不输出文件
def ff_doMixScape4():
    f_doMixScape(specific_dataset, 'Homo sapiens')




unuseList = ["/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb-seq/PRJNA549582/TCell1_nonorm",
"/home/wzt/project/HC-CrisprScreen/poolSC_data/8ECCITE-seq/PRJNA556195/EGR2_KnockOut",
"/home/wzt/project/HC-CrisprScreen/poolSC_data/17SHARE-seq/PRJNA893678/180124_perturb",
"/home/wzt/project/HC-CrisprScreen/poolSC_data/17SHARE-seq/PRJNA893678/210715_combinatorial",
"/home/wzt/project/HC-CrisprScreen/poolSC_data/13arrayBased/CaltechDATA/GPCTP_2019"]



'''
A549  K562   MCF7
/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb-seq/PRJNA831566/RPE1_essential
/home/wzt/project/HC-CrisprScreen/poolSC_data/17SHARE-seq/PRJNA893678/210322_TFAtlas
/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb-seq/PRJNA831566/K562_GW
'''


mtpercent = 10
if __name__ == '__main__':
    print ('hello, world')
    ff_doMixScape()

