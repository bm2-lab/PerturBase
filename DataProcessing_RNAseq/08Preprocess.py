from myUtil import *
import scanpy as sc
import warnings
import anndata as ad
warnings.filterwarnings('ignore')
from scipy import sparse
from mudata import MuData # type: ignore

### ECCITE-seq
basepath =  '/home/wzt/project/HC-CrisprScreen/poolSC_data/8ECCITE-seq/'

def fun1(adata, isWrite = True):
    if isWrite: adata.write_h5ad('raw.h5ad')
    print (adata.shape)
    adata = adata[~adata.obs['gene'].isin(['None', 'CTRL'])]

    tmp = [i.split(',') for i in adata.obs['gene']]
    tmp = [i for i in tmp if i != 'CTRL' and i != 'None']
    tmp = np.unique([i for j in tmp for i in j])
    print (tmp)
    print (len(tmp))

'''
'''
def PRJNA797432():
    os.chdir('/home/wzt/project/HC-CrisprScreen/poolSC_data/8ECCITE-seq/PRJNA797432/ORF_screen')
    dat = pd.read_csv('GSM5819660_GEX_counts.csv', sep=',', index_col=0)
    dat = dat.T

    adata = ad.AnnData(X=sparse.csr_matrix(dat.values), obs=pd.DataFrame(index=dat.index), var=pd.DataFrame(index=dat.columns))

    dat1 = pd.read_csv('GSM5819659_ORF_counts.csv', sep=',', index_col=0)
    mydict = dat1.idxmax().to_dict()
    adata.obs['gene'] = [mydict[i] for i in adata.obs_names]
    adata.obs['gene'].replace({'NGFR':'CTRL'}, inplace=True)
    fun1(adata, True)

    dat2 = pd.read_csv('GSM5819658_ADT_counts.csv', sep=',', index_col=0).T
    adata1 = ad.AnnData(X=sparse.csr_matrix(dat2.values), obs=pd.DataFrame(index=dat2.index), var=pd.DataFrame(index=dat2.columns))
    mdata = MuData({'RNA': adata, 'protein': adata1})
    mdata.write("mudata.h5mu")



'''
Antigen-driven EGR2 expression is required for exhausted CD8+ T cell stability and maintenance
Knock out
'''
def PRJNA556195():
    os.chdir('/NFS_home/NFS_home_2/wzt/project/HC-CrisprScreen/poolSC_data/8ECCITE-seq/PRJNA556195/EGR2_KnockOut')
    adata = sc.read_10x_h5('GSM4664613_GEX_filtered_feature_bc_matrix.h5')
    adata.obs['gene'] = 'EGR2'
    fun1(adata)


'''
Multiplexed detection of proteins, transcriptomes, clonotypes and CRISPR perturbations in single cells
'''

### 第二篇数据
def PRJNA521522():   ### 把protein 和 RNA数据取 使用mudata函数弄成多模态数据
    os.chdir('/NFS_home/NFS_home_2/wzt/project/HC-CrisprScreen/poolSC_data/8ECCITE-seq/PRJNA521522/K')
    RNA= pd.read_csv('GSM3596090_K-cDNA.txt', sep=' ', index_col=0)
    protein = pd.read_csv('GSM3596091_K-ADT-count.csv', sep=',', index_col=0)
    protein = protein.iloc[[0, 1], :]; protein.index = ['CD29', 'CD46']
    cells = [i for i in protein.columns if i in RNA.columns]
    RNA = RNA[cells]; protein = protein[cells]
    RNA = RNA.T;  protein = protein.T

    adata_RNA = ad.AnnData(X=sparse.csr_matrix(RNA.values), obs=pd.DataFrame(index=RNA.index), var=pd.DataFrame(index=RNA.columns))
    adata_protein = ad.AnnData(X=sparse.csr_matrix(protein.values), obs=pd.DataFrame(index=protein.index), var=pd.DataFrame(index=protein.columns))

    GDO = pd.read_csv('GSM3596093_K-GDO-count.csv', sep=',', index_col=0)
    GDO_gene = GDO.iloc[:13, ]
    GDO_gene.index = ['CD29', 'CD29', 'CD46', 'CD46', 'CD46', 'JAK1', 'JAK1', 'JAK1', 
    'CTRL', 'CTRL', 'P53', 'P53', 'P53']

    mydict1 = GDO_gene.idxmax().to_dict()

    adata_RNA.obs['gene'] = [mydict1.get(i, 'None') for i in adata_RNA.obs_names]
    adata_protein.obs['gene'] = [mydict1.get(i, 'None') for i in adata_protein.obs_names]
    fun1(adata_RNA, isWrite= True);  fun1(adata_protein, isWrite= False)
    
    mdata = MuData({'RNA': adata_RNA, 'protein': adata_protein})
    mdata.write("mudata.h5mu")


'''
第三篇 
'''
def PRJNA641353():
    import pertpy as pt #type: ignore
    os.chdir('/home/wzt/project/HC-CrisprScreen/poolSC_data/8ECCITE-seq/PRJNA641353/ECCITE')
    mdata = pt.dt.papalexi_2021()
    adata_RNA = mdata['rna']
    adata_RNA.obs = adata_RNA.obs[['replicate', 'gene_target', 'Phase', 'NT']]
    adata_RNA.obs.columns = ['replicate', 'gene', 'Phase', 'sgRNA']
    adata_RNA.obs['gene'].replace({'NT': 'CTRL'}, inplace=True)
    adata_RNA.obs['sgRNA'].replace({'NT': 'CTRL'}, inplace=True)
    
    adata_protein = mdata['adt']
    adata_protein.obs = adata_protein.obs[['replicate', 'gene_target', 'Phase', 'NT']]
    adata_protein.obs.columns = ['replicate', 'gene', 'Phase', 'sgRNA']
    adata_protein.obs['gene'].replace({'NT': 'CTRL'}, inplace=True)
    adata_protein.obs['sgRNA'].replace({'NT': 'CTRL'}, inplace=True)
    fun1(adata_protein, False)

    mdata = MuData({'RNA': adata_RNA, 'protein': adata_protein})
    mdata.write("mudata.h5mu")



###
if __name__ == '__main__':
    print ('hello, world')
    PRJNA521522()
    