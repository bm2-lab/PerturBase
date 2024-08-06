from myUtil import *
import scanpy as sc
import warnings
import anndata as ad
warnings.filterwarnings('ignore')
from scipy import sparse
from itertools import compress

#### Mosaic-seq

def fun1(adata, isWrite = True):
    if isWrite: adata.write_h5ad('raw.h5ad')
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

'''
ADT, Antibody-derived tag sequencing
GDO, guide-derived oligonucleotide in eccite_seq is cell barcode
Multiplexed CRISPR technologies, in which numerous gRNAs or Cas enzymes are expressed at once, 
have facilitated powerful biological engineering applications
'''

'''
Multiplexed Engineering and Analysis of Combinatorial Enhancer Activity in Single Cells
'''

### first dataset
def PRJNA322853_1():
    def PRJNA322853_fun(GSM, MOI, batch):
        def PRJNA322853_funfun(sgRNA):
            if sgRNA == 'None': return 'None'
            else:
                sgRNAs = sgRNA.split(',')
                return ','.join(sgRNA2gene[tmp] for tmp in sgRNAs)          

        def PRJNA322853_funfun1(sgRNA, read):  #### min reads 5
            sgRNAs = sgRNA.split(',')
            reads = read.split(',')
            tmp = [int(i) >=5 for i in reads]
            result = list(compress(sgRNAs, tmp))
            if len(result) == 0: return 'None'
            else:
                return ','.join(sorted(result))

        cell2sgRNA = {}; sgRNA2gene = {}
        filein = '{}_K562_dCas9-KRAB_10_sgRNAs_targeting_TAD2_SE2_enhancers_{}_MOI_batch_{}.counts.txt'.format(GSM, MOI, batch)
        dat = pd.read_csv(filein, sep='\t', index_col=0)
        dat.index = [i.strip().split(';')[0] for i in dat.index]
        dat = dat.iloc[:, 5:]; dat.columns = [i.split('.')[0] for i in dat.columns]
        dat = dat.T
        adata = ad.AnnData(X=sparse.csr_matrix(dat.values),
        obs=pd.DataFrame(index=dat.index),
        var=pd.DataFrame(index=dat.columns))

        filein = '{}_K562_dCas9-KRAB_10_sgRNAs_targeting_TAD2_SE2_enhancers_{}_MOI_batch_{}.cell_barcode_to_sgRNA_barcode.txt'.format(GSM, MOI, batch)
        dat1 = pd.read_csv(filein, sep='\t')
        for cell, sgRNAs, read in zip(dat1['cell_barcode'], dat1['sgRNA_barcodes_detected'], dat1['#reads/sgRNA_barcode']):
            if sgRNAs is np.nan: sgRNA = 'None'
            else: 
                sgRNA = PRJNA322853_funfun1(sgRNAs, read)
            cell2sgRNA[cell]  = sgRNA
        dat2 = pd.read_excel('../Mosaic_mmc2.xlsx')
        dat2.iloc[-3, 2] = 'GFP'; dat2.iloc[-4, 2] = 'GFP'
        dat2['target'] = dat2['target'].apply(lambda x: 'CTRL' if x.startswith('GFP') else x)
        sgRNA2gene = dat2[['target', 'sgRNA_barcode_seq']].set_index('sgRNA_barcode_seq').to_dict()['target']

        adata.obs['sgRNA'] = [cell2sgRNA[cell] for cell in adata.obs_names]
        adata.obs['gene'] = adata.obs['sgRNA'].apply(PRJNA322853_funfun)
        adata.obs['gene'] = adata.obs['gene'].apply(preGene)  #### delete CTRL,CTRL
        adata.obs['MOI'] = MOI; adata.obs['batch'] = batch
        adata.var_names_make_unique()
        return adata
    

    os.chdir('/home/wzt/project/HC-CrisprScreen/poolSC_data/5Mosaic-seq/PRJNA322853/TAD2_SE2')
    GSMs = ['GSM2544738', 'GSM2544739', 'GSM2544740', 'GSM2544741', 'GSM2544742']
    MOIs = ['low', 'low', 'low', 'low', 'high']; batchs = [1, 2, 3, 4, 1]
    mylist = [PRJNA322853_fun(GSM, MOI, batch) for GSM, MOI, batch in zip(GSMs, MOIs, batchs)]
    adata = ad.concat(mylist)
    fun1(adata)

### second dataset
def PRJNA322853_2():
    def PRJNA322853_fun(counts_file, sgRNABarcode_file):
        def PRJNA322853_funfun(sgRNA):
            if sgRNA == 'None': return 'None'
            else:
                sgRNAs = sgRNA.split(',')
                return ','.join(sgRNA2gene[tmp] for tmp in sgRNAs)    

        def PRJNA322853_funfun1(sgRNA, read):  #### min reads 5
            sgRNAs = sgRNA.split(',')
            reads = read.split(',')
            tmp = [int(i) >=5 for i in reads]
            result = list(compress(sgRNAs, tmp))
            if len(result) == 0: return 'None'
            else:
                return ','.join(sorted(result))


        cell2sgRNA = {}; sgRNA2gene = {}
        dat = pd.read_csv(counts_file, sep='\t', index_col=0)
        dat.index = [i.strip().split(';')[0] for i in dat.index]
        dat = dat.iloc[:, 5:]; dat.columns = [i.split('.')[0] for i in dat.columns]
        dat = dat.T
        adata = ad.AnnData(X=sparse.csr_matrix(dat.values),
        obs=pd.DataFrame(index=dat.index),
        var=pd.DataFrame(index=dat.columns))

        dat1 = pd.read_csv(sgRNABarcode_file, sep='\t')
        for cell, sgRNAs, read, in zip(dat1['cell_barcode'], dat1['sgRNA_barcodes_detected'],  dat1['#reads/sgRNA_barcode']):
            if sgRNAs is np.nan: sgRNA = 'None'
            else:  sgRNA = PRJNA322853_funfun1(sgRNAs, read)
            cell2sgRNA[cell]  = sgRNA
        dat2 = pd.read_excel('../Mosaic_mmc2.xlsx')
        dat2.iloc[-3, 2] = 'GFP'; dat2.iloc[-4, 2] = 'GFP'
        dat2['target'] = dat2['target'].apply(lambda x: 'CTRL' if x.startswith('GFP') else x)
        sgRNA2gene = dat2[['target', 'sgRNA_barcode_seq']].set_index('sgRNA_barcode_seq').to_dict()['target']

        adata.obs['sgRNA'] = [cell2sgRNA[cell] for cell in adata.obs_names]
        adata.obs['gene'] = adata.obs['sgRNA'].apply(PRJNA322853_funfun)
        adata.obs['gene'] = adata.obs['gene'].apply(preGene)  #### delete  CTRL,CTRL
        adata.obs['set'] = counts_file.split('_')[5]; adata.obs['batch'] = counts_file.split('_')[-1][0]
        adata.var_names_make_unique()
        return adata

    os.chdir('/home/wzt/project/HC-CrisprScreen/poolSC_data/5Mosaic-seq/PRJNA322853/15SE_71HS')
    files1 = sorted(glob.glob('*counts.txt')); files2 = sorted(glob.glob('*sgRNA_barcode.txt'))
    mylist = [PRJNA322853_fun(counts_file, sgRNABarcode_file) for counts_file, sgRNABarcode_file in zip(files1, files2)]
    adata = ad.concat(mylist)
    fun1(adata)


'''
Global Analysis of Enhancer Targets Reveals Convergent Enhancer-Driven Regulatory Modules
'''
def PRJNA532921():
    def PRJNA532921_fun(h5_file, sgRNA_file):
        def PRJNA532921_funfun(x):
            xs = x.split(';')
            tmp = [sgRNA2gene.get(i, 'None') for i in xs]
            return ';'.join(tmp)
        
        def PRJNA532921_funfun1(x):   #### assign None to CTRL
            if x == 'None': return 'CTRL'
            xs = x.split(';')
            if np.all(np.array(xs) == 'None'): return 'CTRL'
            xs = [i for i in xs if i != 'None']
            return ','.join(xs)

        cell2sgRNA = {}; sgRNA2gene = {}
        adata = sc.read_10x_h5(h5_file)
        
        dat = pd.read_csv(sgRNA_file, sep='\t', header=None)
        tmp = []
        for _, x in dat.iterrows():
            UMIs, Counts = x[3].split(';'), x[4].split(';')
            UMIs = np.array(UMIs)
            Counts = np.array([int(i) for i in Counts])
            UMIs = UMIs[Counts >= 5]   ### min reads 5 
            tmp.append(';'.join(UMIs))
        dat[3] = tmp
        cell2sgRNA = dat[[0, 3]].set_index(0).to_dict()[3]
        
        
        dat1 = pd.read_excel('mmc2.xlsx')
        sgRNA2gene = dat1[['spacer sequence', 'region pos (hg38)']].set_index('spacer sequence').to_dict()['region pos (hg38)']
        adata.obs['sgRNA'] = [cell2sgRNA.get(i[:-2], 'None') for i in adata.obs_names]
        adata.obs['gene'] = adata.obs['sgRNA'].apply(PRJNA532921_funfun)
        adata.obs['batch'] = h5_file[43:46]
        adata.var_names_make_unique()
        adata.obs['gene'] = adata.obs['gene'].apply(PRJNA532921_funfun1)
        return adata

    os.chdir('/home/wzt/project/HC-CrisprScreen/poolSC_data/5Mosaic-seq/PRJNA532921')
    files1 = sorted(glob.glob('*bc_matrices_h5.h5')); files2 = sorted(glob.glob('*sgRNA_UMI.txt'))
    mylist = [PRJNA532921_fun(h5_file, sgRNA_file) for h5_file, sgRNA_file in zip(files1, files2)]
    adata = ad.concat(mylist)
    os.chdir('K562-dCas9-KRAB_5K')
    fun1(adata)


###
if __name__ == '__main__':
    print ('hello, world')
    #PRJNA322853_1()
    #PRJNA322853_2()
    PRJNA532921()
    
