from myUtil import *
import scanpy as sc
import warnings
import anndata as ad
import re
warnings.filterwarnings('ignore')
from scipy import sparse
from mudata import MuData

def fun1(adata, isWrite = True, fileName = 'raw.h5ad'):
    if isWrite: adata.write_h5ad(fileName)
    print (adata.shape)
    adata = adata[~adata.obs['gene'].isin(['None', 'CTRL'])]

    tmp = [i.split(',') for i in adata.obs['gene']]
    tmp = [i for i in tmp if i != 'CTRL' and i != 'None']
    tmp = np.unique([i for j in tmp for i in j])
    print (tmp)
    print (len(tmp))

    #print (len(adata.obs['gene'].unique()))
    #print (len(adata.obs['sgRNA'].unique()))
    #tmp = [i.split(',') for i in adata.obs['sgRNA']]
    #tmp = np.unique([i for j in tmp for i in j])
    #print (len(tmp))

## Direct-capture perturbseq
'''
Combinatorial single-cell CRISPR screens by direct guide RNA capture and targeted sequencing
Direct-capture perturbseq, 主要优化了perturbseq技术，并做了一系列的探讨
exp1-5 探讨了敲除一些基因，对UPR通路的一些影响，其中使用了五种技术建库，并讨论了基因互作的一些分析。
exp6  Cholesterol/DNA repair;
exp7, exp8: CRISPRi 和 CRISPRa，

multiplexing，即在设计sgRNA的时候，同一个质粒载体，对一个基因
同时设计多个sgRNA提高编辑的效率。在cell_identities文件中，把 sgRNA + sgRNA和 sgRNA + control，即non_targeting，进行对比。 
exp9, exp10 对需要测序的基因先进行富集，比如只测定LINCS 978个基因，增加了对这些基因的检测性能。后续的分析方法是否不同

'''
def PRJNA609688():
    os.chdir('/home/wzt/project/HC-CrisprScreen/poolSC_data/9Perturb-seq/PRJNA609688/exp1-5')
    adata = sc.read_10x_mtx('./', cache=True)
    dat = pd.read_csv('cell_identities.csv', sep=',', index_col=0)
    df, adata_obs = dat.align(adata.obs, join="inner", axis=0)
    adata = adata[adata_obs.index, :]
    adata.obs = adata.obs.merge(df, left_index=True, right_index=True)

    adata = adata[:, ~adata.var['gene_ids'].str.startswith('sg')]
    adata.obs['sgRNA'] = adata.obs['guide_identity']
    adata.obs['gene'] = adata.obs['guide_identity'].apply(lambda x: x[2:])
    adata.obs['gene'] = adata.obs['gene'].apply(lambda x: 'CTRL' if x in ('NegCtrl3', 'NegCtrl2') else x)
    fun1(adata)

    tmpPath = '/NFS_home/NFS_home_2/wzt/project/HC-CrisprScreen/poolSC_data/9Perturb-seq/PRJNA609688/'
    for dirName in ['exp6', 'exp9', 'exp10']:
        dirName = os.path.join(tmpPath, dirName)
        os.chdir(dirName)
        adata = sc.read_10x_mtx('./', cache=True)
        dat = pd.read_csv('cell_identities.csv', sep=',', index_col=0)
        df, adata_obs = dat.align(adata.obs, join="inner", axis=0)
        adata = adata[adata_obs.index, :]
        adata.obs = adata.obs.merge(df, left_index=True, right_index=True)
        adata = adata[:, ~adata.var['gene_ids'].str.startswith('sg')]

        adata.obs['gene_A1'] = adata.obs['gene_A'].apply(lambda x: 'CTRL' if x in ('NegCtrl02093', 'NegCtrl4', 'NegCtrl8', 'NegCtrl5', 'NegCtrl1', 'NegCtrl9') else x)
        adata.obs['gene_A1'] = adata.obs['gene_A1'].apply(lambda x: 'HUS1' if x in ('HUS1_2') else x)
        adata.obs['gene_B1'] = adata.obs['gene_B'].apply(lambda x: 'CTRL' if x in ('NegCtrl3') else x)
        adata.obs['gene_B1'] = adata.obs['gene_B1'].apply(lambda x: 'FDPS' if x in ('FDPS_2') else x)
        adata.obs['gene'] = [','.join(sorted((i, j)))for i, j in zip(adata.obs['gene_A1'], adata.obs['gene_B1'])]
        adata.obs = adata.obs[['gene_A', 'gene_B', 'gene']]
        adata.obs['gene'] = adata.obs['gene'].apply(preGene)
        fun1(adata)

    tmpPath = '/NFS_home/NFS_home_2/wzt/project/HC-CrisprScreen/poolSC_data/9Perturb-seq/PRJNA609688/'
    for dirName in ['exp7', 'exp8']:
        dirName = os.path.join(tmpPath, dirName)
        os.chdir(dirName)
        adata = sc.read_10x_mtx('./', cache=True)
        dat = pd.read_csv('cell_identities.csv', sep=',', index_col=0)
        df, adata_obs = dat.align(adata.obs, join="inner", axis=0)
        adata = adata[adata_obs.index, :]
        adata.obs = adata.obs.merge(df, left_index=True, right_index=True)
        adata = adata[:, ~adata.var['gene_ids'].str.startswith('sg')]
 
        adata.obs['gene'] = adata.obs['gene'].apply(lambda x: 'CTRL' if x in ('Non-Targeting') else x)
        adata.obs = adata.obs[['gene']]
        fun1(adata)


#### CROP-seq
'''
Pooled CRISPR screening with single-cell transcriptome readout
直接使用段斌处理好的数据
'''

def PRJNA358686():
    os.chdir('/home/wzt/project/HC-CrisprScreen/poolSC_data/4CROP-seq/PRJNA358686')
    dat = pd.read_csv('perturb_matrix.tsv', sep='\t', index_col=0)
    
    dat = dat.T
    adata = ad.AnnData(X=sparse.csr_matrix(dat.values), obs=pd.DataFrame(index=dat.index), var=pd.DataFrame(index=dat.columns))

    dat1 = pd.read_csv('perturb_lable.tsv', sep='\t', index_col=0)
    dat1 = dat1.T
    cell2gene = dat1.to_dict()['gene']

    dat2 = pd.read_csv('perturb_lable_condition.tsv', sep='\t', index_col=0)
    dat2 = dat2.T
    cell2condition = dat2.to_dict()['condition']

    dat3 = pd.read_csv('perturb_lable_sgRNA.tsv', sep='\t', index_col=0)
    dat3 = dat3.T
    cell2sgRNA = dat3.to_dict()['grna']

    adata.obs['gene'] = [cell2gene.get(i, 'None') for i in adata.obs_names]
    adata.obs['condition'] = [cell2condition.get(i, 'None') for i in adata.obs_names]
    adata.obs['sgRNA'] = [cell2sgRNA.get(i, 'None') for i in adata.obs_names]
    fun1(adata)


## improved crop-seq
'''
直接使用段斌处理好的数据
'''
def PRJNA428344():
    def myfun2(x):
        if x == 'BLANK_CTRL': return 'None'
        else:
            xs = x.split(',')
            if len(xs) == 0: return 'None'
            elif len(xs) == 1 and 'CTRL' in xs: return 'CTRL'
            elif len(xs) == 1 and 'CTRL' not in xs: 
                return ','.join(sorted(xs))
            else: 
                xs = [i for i in xs if i != 'CTRL']
                return ','.join(sorted(xs))
    os.chdir('/home/wzt/project/HC-CrisprScreen/poolSC_data/6CROP-seq/PRJNA428344/backup')
    dat = pd.read_csv('perturb_matrix.tsv', sep='\t', index_col=0)
    
    dat = dat.T
    adata = ad.AnnData(X=sparse.csr_matrix(dat.values), obs=pd.DataFrame(index=dat.index), var=pd.DataFrame(index=dat.columns))

    dat1 = pd.read_csv('perturb_lable.tsv', sep='\t', index_col=0, header=None)
    cell2gene = dat1.to_dict()[1]

    dat2 = pd.read_csv('perturb_lable_condition.tsv', sep='\t', index_col=0, header=None)
    cell2condition = dat2.to_dict()[1]

    adata.obs['gene'] = [cell2gene.get(i, 'None') for i in adata.obs_names]
    adata.obs['gene'] = adata.obs['gene'].apply(myfun2)
    adata.obs['condition'] = [cell2condition.get(i, 'None') for i in adata.obs_names]
    adata1 = adata[adata.obs['condition'] == 'mock']
    adata2 = adata[adata.obs['condition'] == 'dox_100nm']
    adata1.write_h5ad('mock.raw.h5ad')
    adata2.write_h5ad('dox.raw.h5ad')
    fun1(adata)



## perturb-ATAC
'''
Coupled Single-Cell CRISPR Screening and Epigenomic Profiling Reveals Causal Gene Regulatory Networks

'''
### 第一个数据集
def PRJNA478043():
    os.chdir('/home/wzt/project/HC-CrisprScreen/poolSC_data/7Perturb-ATAC/PRJNA478043/BCell')
    peak = pd.read_csv('GSE116249_Peak_counts_GM12878_experiment2.txt', sep='\t')
    peak.index = ['_'.join([i, str(j), str(k)]) for i,j,k in zip(peak['chr'], peak['start'], peak['stop'])]
    peak = peak.iloc[:, 3:]
    peak.columns = ['-'.join(i.split('-')[:5]) for i in peak.columns]

    gbc = pd.read_csv('GSE116285_GBC_counts_GM12878_experiment2.txt', sep='\t', index_col=0)
    gbc.index = ['-'.join(i.split('-')[:5]) for i in gbc.index]
    gbc.columns = ['CTRL' if i in ('NT1', 'NT2') else i for i in gbc.columns]
    myrows, mycols = np.where(gbc >= 1000)  ###大于1000，则赋值，沿用文献原文的方法。
    mydict = defaultdict(set); mydict1 = defaultdict()
    for i, j in zip(myrows, mycols):
        cell = gbc.index[i]; gene = gbc.columns[j]
        mydict[cell].add(gene)
    for cell in mydict:
        mydict1[cell] = ','.join(mydict[cell])

    peak = peak.T
    adata = ad.AnnData(X=sparse.csr_matrix(peak.values), obs=pd.DataFrame(index=peak.index), var=pd.DataFrame(index=peak.columns))
    adata.obs['gene'] = [mydict1.get(i, 'None') for i in adata.obs_names]
    adata.obs['sgRNA'] = adata.obs['gene']
    fun1(adata)

#### 第二个数据集
    os.chdir('/NFS_home/NFS_home_2/wzt/project/HC-CrisprScreen/poolSC_data/7Perturb-ATAC/PRJNA478043/Keratinocyte')
    peakKO = pd.read_csv('GSE116249_Peak_counts_Keratinocyte_KO.txt', sep='\t', index_col=0)
    peakWT = pd.read_csv('GSE116248_Peak_counts_Keratinocyte_WT.txt', sep='\t', index_col=0)
    peakKO.columns = ['-'.join(i.split('-')[:2]) for i in peakKO.columns]
    peakWT.columns = ['-'.join(i.split('-')[:2]) for i in peakWT.columns]
    
    
    sgRNA = pd.read_csv('GSE116285_sgRNA_counts_Keratinocyte.txt', sep='\t', index_col=False)
    sgRNA.set_index('Cell', inplace=True)
    sgRNA.index = ['-'.join(i.split('-')[:2]) for i in sgRNA.index]
    sgRNA.columns = ['CTRL' if i in ('NT1', 'NT2') else i for i in sgRNA.columns]
    myrows, mycols = np.where(sgRNA >= 1000)
    mydict = defaultdict(set); mydict1 = defaultdict()
    for i, j in zip(myrows, mycols):
        cell = sgRNA.index[i]; gene = sgRNA.columns[j]
        mydict[cell].add(gene)
    for cell in mydict:
        mydict1[cell] = ','.join(mydict[cell])

    peak = pd.concat([peakKO, peakWT], axis=1)

    peak = peak.T
    adata = ad.AnnData(X=sparse.csr_matrix(peak.values), obs=pd.DataFrame(index=peak.index), var=pd.DataFrame(index=peak.columns))
    adata.obs['gene'] = [mydict1.get(i, 'None') for i in adata.obs_names[:peakKO.shape[1]]] + ['WT'] * peakWT.shape[1]
    adata.obs['sgRNA'] = adata.obs['gene']
    fun1(adata)



### SHARE-seq
'''
张峰实验室最近的文章, 一共15个数据集，其中
GSE216595, GSE217066, GSE217215, GSE217460   有数据可以下载，但是不一定是扰动数据

GSE216457, GSE216463, GSE216602, GSE216601
GSE216477, GSE216479, GSE218506, GSE218562, GSE218789,
GSE219000, GSE219058数据不符合要求
'''

###  180124_perturb
def PRJNA893678():
    os.chdir('/home/wzt/project/HC-CrisprScreen/poolSC_data/17SHARE-seq/PRJNA893678/GSE216595')
    adata = sc.read_h5ad('GSE216595_180124_perturb.h5ad')
    adata.obs.drop(labels=['temp'], axis=1, inplace=True)
    adata.obs['gene'] = adata.obs['TF'].apply(lambda x: x.split('-')[0])
    fun1(adata)


####  210715_combinatorial 数据集
def PRJNA893678_2():
    ### c, 第二个数据集, 因为h5ad数据是矫正过后的数据，整理原始数据
    ### 第一步，获取保留的cells
    os.chdir('/home/wzt/project/HC-CrisprScreen/poolSC_data/17SHARE-seq/PRJNA893678/GSE217066')
    def myfun1(filein, adata, fileout):
        tmp = pd.read_csv(filein, sep='\t', chunksize=1000, index_col=0, compression='gzip')
        mylist = []
        for chunk in tmp:
            keepcells = [i for i in adata.obs_names if i in chunk.columns]
            chunk_sub = chunk[keepcells]
            mylist.append(chunk_sub)
        dat = pd.concat(mylist, axis=0)
        dat.to_csv(fileout, sep='\t', header=True, index=True)
    adata1 = sc.read_h5ad('GSE217066_210715_combinatorial_subsample.h5ad')
    adata1.obs.index = [i[:-2] for i in adata1.obs_names]
    for i in range(1, 10):
        filein = glob.glob('*210715_combinatorial_S{}_counts.csv.gz'.format(i))[0]
        fileout = 'S{}_subset_counts.tsv'.format(i)
        myfun1(filein, adata1, fileout)


    def myfun2(x):
        xs = x.split(',')
        xs = [i.split('-')[1] for i in xs]
        if len(xs) == 0: return 'None'
        elif len(xs) == 1 and 'GFP' in xs: return 'CTRL'
        elif len(xs) == 1 and 'GFP' not in xs:  return ','.join(sorted(xs))
        else: 
            xs = [i for i in xs if i != 'GFP']
            return ','.join(sorted(xs))
    ###第二步，整理成h5ad格式
    files = glob.glob('*subset_counts.tsv')
    dataList = []
    for file in files:
        dat_tmp = pd.read_csv(file, sep='\t', index_col=0)
        dataList.append(dat_tmp)
    dat = pd.concat(dataList, axis=1, join='inner')
    adata = ad.AnnData(X=sparse.csr_matrix(dat.T.values), obs=pd.DataFrame(index=dat.columns), var=pd.DataFrame(index=dat.index))
    
    df, adata_obs = adata1.obs.align(adata.obs, join="inner", axis=0)
    adata = adata[adata_obs.index, :]
    adata.obs = adata.obs.merge(df, left_index=True, right_index=True)
    adata.obs['gene'] = adata.obs['TF'].apply(myfun2)
    fun1(adata)

'''
GSE217215
'''
def PRJNA893678_3():
    os.chdir('/home/wzt/project/HC-CrisprScreen/poolSC_data/17SHARE-seq/PRJNA893678/GSE217215')
    #### 首先利用R包读取rds获取子集数据，要不然原始数据太大
    ### library('Matrix)
    ### dat = readRDS('GSE217215_201218_ATAC_subsample.rds')  
    ##  write.table(colnames(dat@assays$RNA@data), file = 'RNA_colnames.tsv', sep='\t', quote = F, row.names = F, col.names = F)
    ##  write.table(rownames(dat@assays$RNA@data), file = 'RNA_rownames.tsv', sep='\t', quote = F, row.names = F, col.names = F)
    ##  writeMM(dat@assays$RNA@counts, file = "RNA.mtx")
    
    #writeMM(dat@assays$ATAC@counts, file = "ATAC.mtx")
    #write.table(colnames(dat@assays$ATAC@data), file = 'ATAC_colnames.tsv', sep='\t', quote = F, row.names = F, col.names = F)
    #write.table(rownames(dat@assays$ATAC@data), file = 'ATAC_rownames.tsv', sep='\t', quote = F, row.names = F, col.names = F)
    import scipy.io as spio
    metaData = pd.read_csv('metaData.tsv', sep='\t', index_col=0)
    sparse_matrix = spio.mmread("RNA.mtx")
    tmp1 = pd.read_csv('RNA_colnames.tsv', header=None, sep='\t')
    tmp2 = pd.read_csv('RNA_rownames.tsv', header=None, sep='\t')
    adata = ad.AnnData(sparse.csr_matrix(sparse_matrix.T), obs=pd.DataFrame(index=tmp1[0].values), var=pd.DataFrame(index=tmp2[0].values))
    
    df, adata_obs = metaData.align(adata.obs, join="inner", axis=0)
    adata = adata[adata_obs.index, :]
    adata.obs = adata.obs.merge(df, left_index=True, right_index=True)

    adata.obs['gene'] = adata.obs['TF'].apply(lambda x: x.split('-')[1])
    adata.obs['gene'] = adata.obs['gene'].apply(lambda x: 'CTRL' if x in ('GFP', 'mCherry') else x)
    fun1(adata)

    sparse_matrix = spio.mmread("ATAC.mtx")
    tmp1 = pd.read_csv('ATAC_colnames.tsv', header=None, sep='\t')
    tmp2 = pd.read_csv('ATAC_rownames.tsv', header=None, sep='\t')
    adata = ad.AnnData(sparse.csr_matrix(sparse_matrix.T), obs=pd.DataFrame(index=tmp1[0].values), var=pd.DataFrame(index=tmp2[0].values))
    
    df, adata_obs = metaData.align(adata.obs, join="inner", axis=0)
    adata = adata[adata_obs.index, :]
    adata.obs = adata.obs.merge(df, left_index=True, right_index=True)

    adata.obs['gene'] = adata.obs['TF'].apply(lambda x: x.split('-')[1])
    adata.obs['gene'] = adata.obs['gene'].apply(lambda x: 'CTRL' if x in ('GFP', 'mCherry') else x)
    fun1(adata)

###  GSE217460
def PRJNA893678_4():
    os.chdir('/home/wzt/project/HC-CrisprScreen/poolSC_data/17SHARE-seq/PRJNA893678/210322_TFAtlas')
    adata = sc.read_h5ad('GSE217460_210322_TFAtlas_subsample.Raw.h5ad')
    adata.obs['gene'] = adata.obs['TF'].apply(lambda x: x.split('-')[1])
    adata.obs['gene'] = adata.obs['gene'].apply(lambda x: 'CTRL' if x in ('GFP', 'mCherry') else x)

    adata.obs = adata.obs[['batch', 'gene']]
    genes = adata.obs['gene'].unique()
    tmp = [adata[adata.obs['gene'] == i][:500] for i in genes]  ### 每个扰动最多保留500个细胞
    result = ad.concat(tmp)

    tmp = result.obs['gene'].value_counts()   ###去除细胞数量太少的扰动
    genes = list(tmp[tmp >= 100].index)  ###
    if 'CTRL' not in genes: genes += ['CTRL']
    result = result[result.obs['gene'].isin(genes), :]

    result.write_h5ad('raw.h5ad')
    fun1(adata, isWrite=True)

'''
Ultra-high throughput single-cell RNA sequencing by combinatorial fluidic indexing
优化了droplet的问题，本质上是array-based, 特点是对每个转录谱加上标记
'''
def PRJNA713314():
    os.chdir('/home/wzt/project/HC-CrisprScreen/poolSC_data/19scifi-RNA-seq')
    dat = pd.read_csv('GSM5151370_PD213_scifi_2_CRISPR-TCR_77300_MeOH-cells.csv', sep=',')
    dat = dat[['gRNA_ID', 'gRNA_seq', 'TCR_status', 'plate_well']]
    dat['gene'] = list(dat['gRNA_ID'].apply(lambda x: x.split('_')[0]))
    dat['gene'] = dat['gene'].apply(lambda x: 'CTRL' if x.startswith('CTRL') else x)
    adata = sc.read_h5ad('GSM5151370_PD213_scifi_2_CRISPR-TCR_77300_MeOH-cells.h5ad')
    adata.obs['index'] = adata.obs_names
    adata.obs = adata.obs.merge(dat, left_on='plate_well', right_on='plate_well')
    adata.obs.set_index('index', inplace=True)
    adata1 = adata[adata.obs['TCR_status'] == 'stimulated']
    adata2 = adata[adata.obs['TCR_status'] == 'unstimulated']
    fun1(adata1); fun1(adata2)

'''
Efficient combinatorial targeting of RNA transcripts in single cells with Cas13 RNA Perturb-seq
主要是优化了组合扰动
GSE213957
'''

def PRJNA883380():
    ### 数据集2
    ### 先进行预处理 ADT文件
    # os.chdir('/home/wzt/project/HC-CrisprScreen/poolSC_data/20CaRPool-seq/PRJNA883380/THP1/ADT/')
    # names = ['1', '2', '3', '4']
    # for i in names:
    #     filein = 'THP1-CaRPool-seq_and_HEK293FTstabRNA.ADT{}.features.tsv'.format(i)
    #     fileout = 'ADT{}.genes.tsv'.format(i)
    #     dat = pd.read_csv(filein, sep='\t', header=None)
    #     mylist = []
    #     for i in dat[0]:
    #         tmp = i.split('-')
    #         if len(tmp) == 1 or len(tmp) == 2: mylist.append(tmp[0])
    #         elif len(tmp) == 3: mylist.append('-'.join(tmp[:2]))
    #     dat[1] = mylist; dat[2] = dat[1]; dat = dat[[1, 2]]
    #     dat.to_csv(fileout, sep='\t', header= False, index=False)

    os.chdir('/home/wzt/project/HC-CrisprScreen/poolSC_data/20CaRPool-seq/PRJNA883380/THP1')
    def ADTfun(prefix):
        ADT = sc.read_10x_mtx(path = 'ADT/', prefix = 'ADT{}.'.format(prefix), cache=True)
        ADT = ADT[:, ADT.var['gene_ids'] != 'unmapped']
        ADT.obs_names = ['L{}_'.format(prefix) + i for i in ADT.obs_names]
        return ADT
    ADT1 = ADTfun('1'); ADT2 = ADTfun('2')
    ADT3 = ADTfun('3'); ADT4 = ADTfun('4')
    adata1 = ad.concat([ADT1, ADT2, ADT3, ADT4], join='outer', merge='same')

    def GEXfun(prefix):
        GEX = sc.read_10x_mtx(path = 'GEXGDO/', prefix = 'GEXGDO{}.'.format(prefix), cache=True)
        dat = pd.read_csv('GEXGDO/GEXGDO{}.genes.tsv'.format(prefix), sep='\t', header=None)
        GEX.var['annotation'] = dat[2].values
        GEX.obs_names = ['L{}_'.format(prefix) + i[:-2] for i in GEX.obs_names]
        return GEX
    GEX1 = GEXfun('1'); GEX2 = GEXfun('2')
    GEX3 = GEXfun('3'); GEX4 = GEXfun('4')
    adata2 = ad.concat([GEX1, GEX2, GEX3, GEX4], join='outer', merge='same')
    adata2 = adata2[:, adata2.var['annotation'] != 'CRISPR Guide Capture']

    dat = pd.read_csv('GSE213957_THP1-CaRPool-seq.metadata.tsv', sep='\t')
    dat.index = [i[:-2] for i in dat.index]
    dat = dat[['GenePair', 'Phase']]

    def mergeFun(dat, adata):
        df, adata_obs = dat.align(adata.obs, join="inner", axis=0)
        adata = adata[adata_obs.index, :]
        adata.obs = adata.obs.merge(df, left_index=True, right_index=True)
        adata.obs['gene1'] = adata.obs['GenePair'].apply(lambda x: ','.join(sorted(x.split('_'))))
        return adata
    adata1 = mergeFun(dat, adata1)
    adata2 = mergeFun(dat, adata2)


    def myfun2(x):
        xs = x.split(',')
        xs1 = [i for i in xs if i != 'NT']
        if len(xs1) == 0:
            return 'CTRL'
        else:
            return ','.join(sorted(xs1))
    adata1.obs['gene'] = adata1.obs['gene1'].apply(myfun2)
    adata1.obs = adata1.obs[['GenePair', 'Phase', 'gene']]

    adata2.obs['gene'] = adata2.obs['gene1'].apply(myfun2)
    adata2.obs = adata2.obs[['GenePair', 'Phase', 'gene']]

    mdata = MuData({'RNA': adata2, 'protein': adata1})  ### 多模态数据
    mdata.write("mudata.h5mu")

    fun1(adata2)


'''
TAP-seq  只测序一部分基因
'''

### 第一批数据集，扰动的是基因，
def PRJNA559094_1():
    def PRJNA559094_fun(dirName):
        os.chdir(dirName)
        count_files = sorted(glob.glob('*counts.csv'))
        pertStatus_files = sorted(glob.glob('*pertStatus.csv'))
        adata_list = []
        for count_file, pertStatus_file in zip(count_files, pertStatus_files):
            dat1 = pd.read_csv(count_file, sep=',', index_col=0)
            dat2 = pd.read_csv(pertStatus_file, sep=',', index_col=0)
            
            dat1.columns = [re.split(r'\.|-', i)[-1] for i in dat1.columns]            
            dat2.columns = [re.split(r'\.|-', i)[-1] for i in dat2.columns]
            dat2.index = ['CTRL' if i.startswith('non-targeting') else i for i in dat2.index]
            dat2.index = ['CTRL' if i.startswith('NT') else i for i in dat2.index]
            dat2.index = [re.split(r'\-|_', i)[0] for i in dat2.index]
            dat1 = dat1.T; dat2 = dat2.T
            adata = ad.AnnData(X=sparse.csr_matrix(dat1.values), obs=pd.DataFrame(index=dat1.index), var=pd.DataFrame(index=dat1.columns))

            myrows, mycols = np.where(dat2 == 1)  ## 对细胞的扰动标签进行分配
            mydict = defaultdict(set); mydict1 = defaultdict()
            for i, j in zip(myrows, mycols):
                cell = dat2.index[i]; gene = dat2.columns[j]
                mydict[cell].add(gene)
            for cell in mydict:
                tmp = mydict[cell]
                if len(tmp) == 0:
                    mydict1[cell] = 'None'
                elif len(tmp) == 1 and 'CTRL' in tmp:
                    mydict1[cell] = 'CTRL'
                elif len(tmp) == 1 and 'CTRL' not in tmp:
                    mydict1[cell] = ','.join(sorted(tmp))
                else:
                    tmp = [i for i in tmp if i != 'CTRL']
                    mydict1[cell] = ','.join(sorted(tmp))
            adata.obs['gene'] = [mydict1.get(i, 'None') for i in adata.obs.index]
            adata_list.append(adata)    
        adata = ad.concat(adata_list, join='inner', axis=0, merge='same')
        fun1(adata)
    
    #PRJNA559094_fun(dirName='/home/wzt/project/HC-CrisprScreen/poolSC_data/16TAP-seq/TAP_DIFFEX')
    #PRJNA559094_fun(dirName='/home/wzt/project/HC-CrisprScreen/poolSC_data/16TAP-seq/WTX_DIFFEX')
    #PRJNA559094_fun(dirName='/home/wzt/project/HC-CrisprScreen/poolSC_data/16TAP-seq/L1000')
    PRJNA559094_fun(dirName='/home/wzt/project/HC-CrisprScreen/poolSC_data/16TAP-seq/REDESIGN')


### 第二批数据集，扰动的是enhancer，
def PRJNA559094_2():
    def PRJNA559094_fun2(dirName):
        os.chdir(dirName)
        count_files = sorted(glob.glob('*counts.csv'))
        pertStatus_files = sorted(glob.glob('*pertStatus.csv'))
        adata_list = []
        for count_file, pertStatus_file in zip(count_files, pertStatus_files):
            dat1 = pd.read_csv(count_file, sep=',', index_col=0)
            dat2 = pd.read_csv(pertStatus_file, sep=',', index_col=0)
            
            dat1.columns = [re.split(r'\.|-', i)[-1] for i in dat1.columns]            
            dat2.columns = [re.split(r'\.|-', i)[-1] for i in dat2.columns]
            tmp = [False if i.startswith('CROPseq_dCas9_DS') else True for i in dat1.index]
            dat1 = dat1[tmp]
            dat2.index = ['CTRL' if i.startswith('non-targeting') else i for i in dat2.index]
            dat1 = dat1.T; dat2 = dat2.T
            adata = ad.AnnData(X=sparse.csr_matrix(dat1.values), obs=pd.DataFrame(index=dat1.index), var=pd.DataFrame(index=dat1.columns))

            myrows, mycols = np.where(dat2 == 1)  ## 对细胞的扰动标签进行分配
            mydict = defaultdict(set); mydict1 = defaultdict()
            for i, j in zip(myrows, mycols):
                cell = dat2.index[i]; gene = dat2.columns[j]
                mydict[cell].add(gene)
            for cell in mydict:
                tmp = mydict[cell]
                if len(tmp) == 0:
                    mydict1[cell] = 'None'
                elif len(tmp) == 1 and 'CTRL' in tmp:
                    mydict1[cell] = 'CTRL'
                elif len(tmp) == 1 and 'CTRL' not in tmp:
                    mydict1[cell] = ','.join(sorted(tmp))
                else:
                    tmp = [i for i in tmp if i != 'CTRL']
                    mydict1[cell] = ','.join(sorted(tmp))
            adata.obs['gene'] = [mydict1.get(i, 'None') for i in adata.obs.index]
            adata_list.append(adata)    
        adata = ad.concat(adata_list, join='inner', axis=0, merge='same')
        fun1(adata)
    
    PRJNA559094_fun2('/home/wzt/project/HC-CrisprScreen/poolSC_data/16TAP-seq/PRJNA559094/SCREEN_chr8')
    PRJNA559094_fun2('/home/wzt/project/HC-CrisprScreen/poolSC_data/16TAP-seq/PRJNA559094/SCREEN_chr11')


#### 处理全基因组数据 







###
if __name__ == '__main__':
    print ('hello, world')
    #PRJNA559094_1()
    #PRJNA559094_2()
    #PRJNA883380()
    PRJNA609688()

    