from myUtil import *
import scanpy as sc
import warnings
import anndata as ad
warnings.filterwarnings('ignore')
from scipy import sparse

def fun1(adata, isWrite = True):
    if isWrite: adata.write_h5ad('raw.h5ad')
    print (adata.shape)
    adata = adata[~adata.obs['gene'].isin(['None', 'CTRL'])]

    tmp = [i.split(',') for i in adata.obs['gene']]
    tmp = [i for i in tmp if i != 'CTRL' and i != 'None']
    tmp = np.unique([i for j in tmp for i in j])
    print (tmp)
    print (len(tmp))

############ 01 perturb seq

'''
Massively parallel phenotyping of coding variants in cancer with Perturb-seq
'''

def PRJNA679579(gene='TP53'):
    os.chdir('PRJNA679579/{}'.format(gene))
    filein = 'A549_{}.rawcounts.matrix.mtx'.format(gene)
    fileout = 'matrix.mtx'
    with open(filein, 'r') as fin, open(fileout, 'w') as fout:
        fout.write(fin.readline()); fout.write(fin.readline());
        for line in fin:
            lines = line.strip().split(' ')
            fout.write('{} {} {}\n'.format(lines[1], lines[0], lines[2]))
    
    filein = 'A549_{}.rawcounts.genes.csv'.format(gene)
    fileout = 'genes.tsv'
    dat = pd.read_csv(filein, sep=' ', header=None)
    dat['temp'] = dat[0]
    dat.to_csv(fileout, sep='\t', header=False, index=False)

    filein = 'A549_{}.variants2cell.csv'.format(gene)
    dat = pd.read_csv(filein, sep='\t')
    dat.index = dat['cell']
    adata = sc.read_10x_mtx('./', cache=True)
    df, adata_obs = dat.align(adata.obs, join="inner", axis=0)
    adata = adata[adata_obs.index, :]
    adata.obs = adata.obs.merge(df, left_index=True, right_index=True)
    adata.obs['gene'] = adata.obs['variant.detailed_multi']
    adata.obs = adata.obs[['batch', 'gene']]
    adata.obs['gene'].replace({'unassigned': 'None'}, inplace=True)
    adata.obs['gene'].replace({'WT': 'CTRL'}, inplace=True)
    genes = list(adata.obs['gene'].value_counts()[adata.obs['gene'].value_counts() >= 100].index)
    adata = adata[adata.obs['gene'].isin(genes)]
    fun1(adata, isWrite = False)

def f_PRJNA679579():
    PRJNA679579("KRAS")
    PRJNA679579("TP53")
    
'''
In vivo Perturb-Seq reveals neuronal and glial abnormalities associated with autism risk genes
loss of function
'''

def PRJNA663571(dirName):
    os.chdir(dirName)
    sample = os.path.basename(dirName)
    mydict = {}; mydict1 = {}; mydict2 = {}
    dat2 = pd.read_csv('../NIHMS1666890-supplement-Table_S5.csv', sep=',', skiprows=1)
    mydict1 = dat2[['Perturbation barcode', 'gene']].set_index('Perturbation barcode').to_dict()['gene']
    dat1 = pd.read_csv('UMI.Counts.csv', sep=',', index_col=0)
    for index, i in dat1.iterrows():
        i.sort_values(ascending=False, inplace=True)
        if i[1] * 1.3 < i[0] and i[0] >=1:  ###
            mydict[index] = i.index[0]
    for cell in mydict:
        barcode = mydict[cell]
        if barcode in mydict1:
            gene = mydict1[barcode]
            mydict2[cell + '-1'] = gene
    adata=  sc.read_10x_h5('sample.h5')
    adata.obs['gene'] = [mydict2.get(i, 'None') for i in  adata.obs_names]
    adata.obs['gene'] = adata.obs['gene'].apply(lambda x: 'CTRL' if x == 'GFP (control)' else x)
    adata = adata[adata.obs['gene'] != 'None']
    adata.obs.index = adata.obs.index + '_{}'.format(sample)
    adata.var_names_make_unique()
    return adata

def f_PRJNA663571():
    os.chdir('PRJNA663571')
    dirNames = ['sample' + str(i) for i in range(1, 19)]
    dirNames = [os.path.join(basepath + '/PRJNA663571/', dirName) for dirName in dirNames]
    datList = [PRJNA663571(dirName) for dirName in dirNames]
    dat = ad.concat(datList, axis=0)
    os.chdir('PRJNA663571')
    dat.write_h5ad('raw.h5ad')

'''
Defining the Teratoma as a Model for Multi-lineage Human Development
'''

def PRJNA656958_fun():
    cell2gene = {}; cell2sgRNA = {}
    dat = pd.read_csv('cbc_gbc_dict.tsv', sep=',', header=None)
    for sgRNA, cells in zip(dat[0], dat[1]):
        gene = sgRNA.split('-')[0]
        cells = cells.split(',')
        for cell in cells:
            cell = cell.strip(); cell = cell[:-2] + '-1'
            cell2gene[cell] = gene; cell2sgRNA[cell] = sgRNA
    return cell2gene, cell2sgRNA

def PRJNA656958_1(dirName):
    os.chdir(dirName)
    def PRJNA656958_1fun(x):
        hg19= sc.read_10x_mtx('{}/hg19/'.format(x), cache=True)
        datfile = pd.read_csv('{}/species_class.csv'.format(x), sep=',')
        datfile = datfile[datfile['call'] != 'Multiplet']
        tmp = hg19.obs_names.isin(list(datfile['barcode'])); hg19 = hg19[tmp]; hg19.obs['species'] = 'hg19'
        hg19.obs['batch'] = x
        return hg19
    dat11_hg19 = PRJNA656958_1fun('m1-1'); dat12_hg19 = PRJNA656958_1fun('m1-2')
    dat21_hg19 = PRJNA656958_1fun('m2-1'); dat22_hg19 = PRJNA656958_1fun('m2-2')
    dat31_hg19 = PRJNA656958_1fun('m3-1'); dat32_hg19 = PRJNA656958_1fun('m3-2')
    dat_hg19 = ad.concat([dat11_hg19, dat12_hg19, dat21_hg19, dat22_hg19, dat31_hg19, dat32_hg19])
    cell2gene, cell2sgRNA = PRJNA656958_fun()

    dat_hg19.obs['gene'] = [cell2gene.get(i, 'None') for i in dat_hg19.obs_names]
    dat_hg19.obs['gene'] = dat_hg19.obs['gene'].apply(lambda x: 'CTRL' if x.startswith('NTC') else x)
    fun1(dat_hg19)

def PRJNA656958_2(dirName):
    os.chdir(dirName)
    def PRJNA656958_2fun(x):
        dat= sc.read_10x_mtx('{}/'.format(x), cache=True)
        datfile = pd.read_csv('{}/species_class.csv'.format(x), sep=',')
        datfile = datfile[datfile['call'] != 'Multiplet']
        mydict = datfile[['barcode', 'call']].set_index('barcode').to_dict()['call']
        tmp = dat.obs_names.isin(list(datfile['barcode'])); dat = dat[tmp]
        dat.obs['species'] = [mydict[i] for i in dat.obs_names]
        dat.obs['batch'] = x
        return dat
    dat11 = PRJNA656958_2fun('m1-1'); dat12 = PRJNA656958_2fun('m1-2')
    dat21 = PRJNA656958_2fun('m2-1'); dat22 = PRJNA656958_2fun('m2-2')
    dat31 = PRJNA656958_2fun('m3-1'); dat32 = PRJNA656958_2fun('m3-2')
    dat = ad.concat([dat11, dat12, dat21, dat22, dat31, dat32])
    cell2gene, cell2sgRNA = PRJNA656958_fun()
    dat.obs['gene'] = [cell2gene.get(i, 'None') for i in dat.obs_names]
    dat.obs['gene'] = dat.obs['gene'].apply(lambda x: 'CTRL' if x.startswith('NTC') else x)
    fun1(dat)
    
       
    dat = sc.read_h5ad('Raw.h5ad')
    tmp1 = dat.obs['species'] == 'hg19'
    tmp2 = dat.var_names.str.startswith('hg19')
    dat_hg19 = dat[tmp1, tmp2]
    dat_hg19.write_h5ad('hg19.raw.h5ad')

def PRJNA656958_3(dirName):
    os.chdir(dirName)
    def PRJNA656958_3fun(x):
        dat= sc.read_10x_mtx('{}/'.format(x), cache=True)
        datfile = pd.read_csv('{}/species_class.csv'.format(x), sep=',')
        datfile = datfile[datfile['call'] != 'Multiplet']
        mydict = datfile[['barcode', 'call']].set_index('barcode').to_dict()['call']
        tmp = dat.obs_names.isin(list(datfile['barcode'])); dat = dat[tmp]
        dat.obs['species'] = [mydict[i] for i in dat.obs_names]
        dat.obs['batch'] = x
        return dat
    dat11 = PRJNA656958_3fun('m1-1'); dat12 = PRJNA656958_3fun('m1-2')
    dat21 = PRJNA656958_3fun('m2-1'); dat22 = PRJNA656958_3fun('m2-2')
    dat = ad.concat([dat11, dat12, dat21, dat22])
    cell2gene, cell2sgRNA = PRJNA656958_fun()
    dat.obs['gene'] = [cell2gene.get(i, 'None') for i in dat.obs_names]
    dat.obs['gene'] = dat.obs['gene'].apply(lambda x: 'CTRL' if x.startswith('NTC') else x)
    fun1(dat)

    dat = sc.read_h5ad('Raw.h5ad')
    tmp1 = dat.obs['species'] == 'hg19'
    tmp2 = dat.var_names.str.startswith('hg19')
    dat_hg19 = dat[tmp1, tmp2]
    dat_hg19.write_h5ad('hg19.raw.h5ad')

def f_PRJNA656958():
    dirNames = ['screen-matrices', 'screen-repool-matrices', 'neural-matrices']
    dirNames = [os.path.join(basepath + '/PRJNA656958/', dirName) for dirName in dirNames]
    #PRJNA656958_1(dirNames[0])
    #PRJNA656958_2(dirNames[1])
    PRJNA656958_3(dirNames[2])


'''
A genome-wide framework for mapping gene regulation via cellular genetic screens
'''
def PRJNA494734(dirName, prefix):
    os.chdir(dirName)
    def fun(x):
        xs = x.split('_')
        genes = [sgRNA2gene[i] for i in xs]
        genes = ['CTRL' if i in NTC else i for i in genes]
        genes = sorted(genes)
        return ','.join(genes)
    tmp = pd.read_csv('GSE120861_grna_groups.txt', sep='\t', header=None)
    sgRNA2gene = tmp.set_index(1).to_dict()[0]
    tmp = pd.read_csv('GSE120861_gene_gRNAgroup_pair_table.txt', sep='\t')
    NTC = list(tmp[tmp['gRNAgroup.start'] == 'NTC']['gRNAgroup'].unique())

    filein = 'GSE120861_{}_screen.genes.txt'.format(prefix); fileout = 'genes.tsv'
    dat = pd.read_csv(filein, sep=' ', header=None); dat['temp'] = dat[0]
    dat.to_csv(fileout, sep='\t', header=False, index=False)
    adata = sc.read_10x_mtx('.', cache=True)
    filein = 'GSE120861_{}_screen.phenoData.txt'.format(prefix)
    metaData = pd.read_csv(filein, sep=' ', header=None)
    metaData = metaData[[1, 6]]
    metaData.columns = ['cell', 'sgRNA']
    metaData = metaData[~metaData['sgRNA'].isna()]
    metaData['gene'] = metaData['sgRNA'].apply(lambda x: fun(x))
    metaData['gene'] = metaData['gene'].apply(preGene)
    metaData.set_index('cell', inplace=True)
    df, adata_obs = metaData.align(adata.obs, join="inner", axis=0)
    adata = adata[adata_obs.index, :]
    adata.obs = adata.obs.merge(df, left_index=True, right_index=True)
    fun1(adata)

def f_PRJNA494734():
    prefixs = ['pilot_lowmoi', 'pilot_highmoi', 'at_scale'][2:]
    dirNames = [os.path.join(basepath + '/PRJNA494734/', prefix) for prefix in prefixs]
    for dirName, prefix in zip(dirNames, prefixs):
        PRJNA494734(dirName, prefix)

'''
Functional single-cell genomics of human cytomegalovirus infection
'''
def PRJNA693896(dirName):
    os.chdir(dirName)
    metaData = pd.read_excel('metaData.xlsx', index_col=0)
    cbc_gbc = pd.read_csv('cbc_gbc_dict.tsv', sep=',', index_col=0)
    dat = pd.merge(cbc_gbc, metaData, left_index=True, right_index=True, suffixes=('','_y'))
    adata = sc.read_10x_mtx('./', cache=True)
    df, adata_obs = dat.align(adata.obs, join="inner", axis=0)
    adata = adata[adata_obs.index, :]
    adata.obs = adata.obs.merge(df, left_index=True, right_index=True)
    adata.obs['sgRNA'] = adata.obs['guide_identity']
    adata.obs['gene'] = adata.obs['sgRNA'].apply(lambda x : x.split('_')[0])
    adata.obs['gene'] = adata.obs['gene'].apply(lambda x: 'CTRL' if x.startswith('GFP') else x)

    adata = adata[:, adata.var['gene_ids'].str.startswith('ENSG')]
    fun1(adata)

def f_PRJNA693896():
    dirName = 'CRISPRi_perturb_host'
    PRJNA693896(dirName)

'''
'''
def PRJNA549582():
    os.chdir('PRJNA549582')
    dat = pd.read_csv('GSE132959_TCell1_nonorm_matrix.rn.2.txt', sep='\t', index_col=0)
    ensemblID = [i.split(':')[1].strip() for i in dat.index]
    dat.index = [i.split(':')[0].strip() for i in dat.index]
    dat = dat.T
    adata = ad.AnnData(X=sparse.csr_matrix(dat.values), obs=pd.DataFrame(index=dat.index), var=pd.DataFrame(index=dat.columns))
    adata.var['ensembl'] = ensemblID
    adata.obs['gene'] = 'DHX37'
    fun1(adata)

'''
PTPN2 regulates the generation of exhausted CD8+ T cell subpopulations and restrains tumor immunity
'''
def PRJNA554074():
    os.chdir('PRJNA554074')
    adata1 = sc.read_10x_mtx('PTPN2', cache=True)
    adata1.obs['gene']= 'PTPN2'
    adata2 = sc.read_10x_mtx('CTRL', cache=True)
    adata2.obs['gene']= 'CTRL'
    adata = ad.concat(adatas=[adata1, adata2])
    fun1(adata)


'''
Perturb-Seq: Dissecting Molecular Circuits with Scalable Single-Cell RNA Profiling of Pooled Genetic Screens
'''
def PRJNA354362(dirName):
    os.chdir(dirName)
    adata = sc.read_10x_mtx('./', cache=True)
    dat = pd.read_csv('cbc_gbc_dict.tsv', sep='\t', header=None)
    mydict = defaultdict(list); mydict1 = {}
    for cell, gene in zip(dat[0], dat[1]):
        cell = cell.strip()
        mydict[cell].append(gene)
    for cell in mydict: mydict1[cell] = ','.join(mydict[cell])
    adata.obs['gene'] = [mydict1.get(i, 'None') for i in adata.obs_names]
    adata.obs['gene'] = adata.obs['gene'].apply(preGene)
    fun1(adata)


def f_PRJNA354362():
    dirNames = ['dc_0hr', 'dc_3hr', 'k562_tfs_7', 'k562_tfs_13', 'k562_tfs_highmoi', 'k562_ccycle']
    tmpPath = 'PRJNA354362'
    for dirName in dirNames:
        dirName = os.path.join(tmpPath, dirName)
        PRJNA354362(dirName)


'''
Genome-wide CRISPR Screens in Primary Human T Cells Reveal Key Regulators of Immune Function
'''

def PRJNA489369(dirName):
    os.chdir(dirName)
    adata = sc.read_10x_mtx('./', cache=True)
    dat1 = pd.read_csv('cbc_gbc_dict.tsv', sep='\t')
    cell2sgRNA = dat1[['Cell.BC', 'sgRNA']].set_index('Cell.BC').to_dict()['sgRNA']
    cell2gene = dat1[['Cell.BC', 'gene']].set_index('Cell.BC').to_dict()['gene']
    adata.obs['gene'] = [cell2gene.get(i[:-2], "None") for i in adata.obs_names]
    adata.obs['sgRNA'] = [cell2sgRNA.get(i[:-2], "None") for i in adata.obs_names]
    fun1(adata)

def f_PRJNA489369():
    dirNames = ['D1N', 'D1S', 'D2N', 'D2S']
    dirNames = [os.path.join(basepath + '/PRJNA489369/', dirName) for dirName in dirNames]
    for dirName in dirNames:
        PRJNA489369(dirName)


'''
Charting oncogenicity of genes and variants across lineages via multiplexed screens in teratomas
'''

def PRJNA715235(dirName):
    rawName = ['AR-V7', 'BCL-XL', 'BCL2', 'Beta-catenin_S33A-S37A-T41A-S45A',
       'Beta-catenin_S33Y', 'CCNB1', 'CCND1', 'CDK1', 'CDK4', 'CDK4_R24C',
       'CDK6', 'Caspase3_C163A', 'Caspase8_C360A',
       'FLAG-Ikkbeta_S177E-S181E', 'FLAG-MKK6_S207E-T211E',
       'FLAG-MKK7-JNK2', 'FLAG-Rheb_Q64L', 'FLAG-YAP2_8SA', 'GSK3b_K85A',
       'Gli2_Truncation', 'Hras_G12V', 'Hras_G12V-E37G',
       'Ikkalpha_S176E-S180E', 'JNK2-WT-O-E_MAPK9', 'KLF4', 'Kras_G12V',
       'MEK1_S218D-S222D', 'MEK5_DD_S311D-T315D', 'None', 'Notch1_ICD',
       'Notch3_ICD', 'POU5F1', 'PTEN_C124S', 'RHOJ', 'RalA_G23V_full',
       'RalA_G23V_mature-peptide', 'Rgl2-CAAX', 'SV40-Large-T-Antigen',
       'SmoM2_W535L', 'Stat3_A662C-N664C-V667L', 'TGFbetaR1_T204D',
       'c-MYC', 'hTERT', 'mCherry', 'myr-FLAG-AKT1', 'myr-FLAG-MEK5',
       'myr-FLAG-PIK3CA', 'p38-WT-O-E_MAPK14',
       'p53_Dominant-negative_R175H']
    
    usedName = ['AR-V7', 'BCL-XL', 'BCL2', 'CTNNB1',
       'CTNNB1', 'CCNB1', 'CCND1', 'CDK1', 'CDK4', 'CDK4',
       'CDK6', 'CASP3', 'CASP8',
       'IKBKB', 'MAP2K6',
       'MAP2K7_JNK2', 'RHEB', 'YAP2', 'GSK3B',
       'GLI2', 'HRAS', 'HRAS',
       'IKBKA', 'JNK2', 'KLF4', 'KRAS',
       'MEK1', 'MEK5', 'None', 'NOTCH1',
       'NOTCH3', 'POU5F1', 'PTEN', 'RHOJ', 'RALA',
       'RALA', 'RGL2', 'SV40',
       'SmoM2', 'STAT3', 'TGFBRI',
       'MYC', 'TERT', 'CTRL', 'AKT1', 'MEK5',
       'PIK3CA', 'MAPK14',
       'TP53']

    mydict2 = {}
    for i, j in zip(rawName, usedName):  mydict2[i] = j

    os.chdir(dirName)
    adata = sc.read_10x_mtx('./', cache=True)
    dat = pd.read_csv('genotype.tsv', sep=',', header=None)
    mydict = defaultdict(list); mydict1 = {}
    for gene, celljoin in zip(dat[0], dat[1]):
        cells = celljoin.strip().split(',')
        for cell in cells:
            mydict[cell].append(mydict2[gene])
    mydict1 = {}
    for cell in mydict:
        tmp = sorted(mydict[cell])
        if 'CTRL' in tmp and len(tmp) >= 2:
            tmp = [i for i in tmp if i != 'CTRL']
            mydict1[cell + '-1'] = ','.join(sorted(tmp))
        else:
            mydict1[cell + '-1'] = ','.join(sorted(tmp))
    adata.obs['gene'] = [mydict1.get(i, "None") for i in adata.obs_names]
    adata.obs['gene'] = adata.obs['gene'].apply(preGene)
    return adata

def f_PRJNA715235():
    dirName = 'PRJNA715235/Pre_inject'
    adata = PRJNA715235(dirName); fun1(adata)
    
    tmpPath = 'PRJNA715235/Driver_lib'
    dirNames = ['Driver_lib_round1_ter1', 'Driver_lib_round1_ter2', 'Driver_lib_round1_ter3',
    'Driver_lib_round1_ter4', 'Driver_lib_round2_ter1', 'Driver_lib_round2_ter2']
    dirNames = [os.path.join(tmpPath, dirName) for dirName in dirNames]
    adataList = [PRJNA715235(dirName) for dirName in dirNames]
    adata = ad.concat(adataList)
    tmp1 = adata.var_names.str.startswith('hg19')
    adata = adata[:, tmp1]
    adata.var.index = [i.split('_')[1] for i in adata.var.index]

    os.chdir(tmpPath)
    fun1(adata)
    

    tmpPath = 'PRJNA715235/Driver_sub_lib'
    dirNames =  ['Driver_sub_lib_round1_ter1', 'Driver_sub_lib_round1_ter2', 'Driver_sub_lib_round1_ter3', 
    'Driver_sub_lib_round2_ter1', 'Driver_sub_lib_round2_ter2']
    dirNames = [os.path.join(tmpPath, dirName) for dirName in dirNames]
    adataList = [PRJNA715235(dirName) for dirName in dirNames]
    adata = ad.concat(adataList)

    tmp1 = adata.var_names.str.startswith('hg19')
    adata = adata[:, tmp1]
    adata.var.index = [i.split('_')[1] for i in adata.var.index]

    os.chdir(tmpPath)
    fun1(adata)


### scPerturb
'''
Highly multiplexed single-cell RNA-seq by DNA oligonucleotide tagging of cellular proteins
'''
def GehringPachter2019():
    os.chdir('13arrayBased/96-plex_scRNA-seq')
    dat = sc.read_h5ad('GehringPachter2019.h5ad')
    tmp = []
    for _, x in dat.obs.iterrows():
        tmp1 = []
        if x['dose_value'] != 0: tmp1.append('BMP4')
        tmp1.append('EGF_bFGF')
        if x['dose_value_3'] != 0: tmp1.append('Scriptaid:decitabine')
        if len(tmp1) != x['nperts']: tmp1.append('retinoic acid')
        tmp.append(','.join(tmp1))
    dat.obs['gene'] = tmp
    fun1(dat)    

    

#### Mosaic-seq
'''
ADT, Antibody-derived tag sequencing
GDO, guide-derived oligonucleotidee

Multiplexed CRISPR technologies, in which numerous gRNAs or Cas enzymes are expressed at once, 
have facilitated powerful biological engineering applications
'''

'''
Multiplexed Engineering and Analysis of Combinatorial Enhancer Activity in Single Cells
'''

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
    

    os.chdir('5Mosaic-seq/PRJNA322853/TAD2_SE2')
    GSMs = ['GSM2544738', 'GSM2544739', 'GSM2544740', 'GSM2544741', 'GSM2544742']
    MOIs = ['low', 'low', 'low', 'low', 'high']; batchs = [1, 2, 3, 4, 1]
    mylist = [PRJNA322853_fun(GSM, MOI, batch) for GSM, MOI, batch in zip(GSMs, MOIs, batchs)]
    adata = ad.concat(mylist)
    fun1(adata)

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
        adata.obs['gene'] = adata.obs['gene'].apply(preGene)  #### delete CTRL,CTRL
        adata.obs['set'] = counts_file.split('_')[5]; adata.obs['batch'] = counts_file.split('_')[-1][0]
        adata.var_names_make_unique()
        return adata

    os.chdir('5Mosaic-seq/PRJNA322853/15SE_71HS')
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
        
        def PRJNA532921_funfun1(x):
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
            UMIs = UMIs[Counts >= 5]
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

    os.chdir('5Mosaic-seq/PRJNA532921')
    files1 = sorted(glob.glob('*bc_matrices_h5.h5')); files2 = sorted(glob.glob('*sgRNA_UMI.txt'))
    mylist = [PRJNA532921_fun(h5_file, sgRNA_file) for h5_file, sgRNA_file in zip(files1, files2)]
    adata = ad.concat(mylist)
    os.chdir('K562-dCas9-KRAB_5K')
    fun1(adata)



def PRJNA797432():
    os.chdir('8ECCITE-seq/PRJNA797432/ORF_screen')
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
    os.chdir('8ECCITE-seq/PRJNA556195/EGR2_KnockOut')
    adata = sc.read_10x_h5('GSM4664613_GEX_filtered_feature_bc_matrix.h5')
    adata.obs['gene'] = 'EGR2'
    fun1(adata)


'''
Multiplexed detection of proteins, transcriptomes, clonotypes and CRISPR perturbations in single cells
'''

def PRJNA521522():
    os.chdir('8ECCITE-seq/PRJNA521522/K')
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


def PRJNA641353():
    import pertpy as pt #type: ignore
    os.chdir('8ECCITE-seq/PRJNA641353/ECCITE')
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


## Direct-capture perturbseq
'''
Combinatorial single-cell CRISPR screens by direct guide RNA capture and targeted sequencing
Direct-capture perturbseq
'''
def PRJNA609688():
    os.chdir('9Perturb-seq/PRJNA609688/exp1-5')
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

    tmpPath = '9Perturb-seq/PRJNA609688/'
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

    tmpPath = '9Perturb-seq/PRJNA609688/'
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
'''

def PRJNA358686():
    os.chdir('4CROP-seq/PRJNA358686')
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
    os.chdir('6CROP-seq/PRJNA428344/backup')
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
def PRJNA478043():
    os.chdir('7Perturb-ATAC/PRJNA478043/BCell')
    peak = pd.read_csv('GSE116249_Peak_counts_GM12878_experiment2.txt', sep='\t')
    peak.index = ['_'.join([i, str(j), str(k)]) for i,j,k in zip(peak['chr'], peak['start'], peak['stop'])]
    peak = peak.iloc[:, 3:]
    peak.columns = ['-'.join(i.split('-')[:5]) for i in peak.columns]

    gbc = pd.read_csv('GSE116285_GBC_counts_GM12878_experiment2.txt', sep='\t', index_col=0)
    gbc.index = ['-'.join(i.split('-')[:5]) for i in gbc.index]
    gbc.columns = ['CTRL' if i in ('NT1', 'NT2') else i for i in gbc.columns]
    myrows, mycols = np.where(gbc >= 1000)  ### threshold in reference study
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

    os.chdir('7Perturb-ATAC/PRJNA478043/Keratinocyte')
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
GSE216595, GSE217066, GSE217215, GSE217460
'''

###  180124_perturb
def PRJNA893678():
    os.chdir('17SHARE-seq/PRJNA893678/GSE216595')
    adata = sc.read_h5ad('GSE216595_180124_perturb.h5ad')
    adata.obs.drop(labels=['temp'], axis=1, inplace=True)
    adata.obs['gene'] = adata.obs['TF'].apply(lambda x: x.split('-')[0])
    fun1(adata)


def PRJNA893678_2():
    os.chdir('17SHARE-seq/PRJNA893678/GSE217066')
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
    os.chdir('17SHARE-seq/PRJNA893678/GSE217215')
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
    os.chdir('17SHARE-seq/PRJNA893678/210322_TFAtlas')
    adata = sc.read_h5ad('GSE217460_210322_TFAtlas_subsample.Raw.h5ad')
    adata.obs['gene'] = adata.obs['TF'].apply(lambda x: x.split('-')[1])
    adata.obs['gene'] = adata.obs['gene'].apply(lambda x: 'CTRL' if x in ('GFP', 'mCherry') else x)

    adata.obs = adata.obs[['batch', 'gene']]
    genes = adata.obs['gene'].unique()
    tmp = [adata[adata.obs['gene'] == i][:500] for i in genes] 
    result = ad.concat(tmp)

    tmp = result.obs['gene'].value_counts()
    genes = list(tmp[tmp >= 100].index)  ###
    if 'CTRL' not in genes: genes += ['CTRL']
    result = result[result.obs['gene'].isin(genes), :]

    result.write_h5ad('raw.h5ad')
    fun1(adata, isWrite=True)

'''
Ultra-high throughput single-cell RNA sequencing by combinatorial fluidic indexing
'''
def PRJNA713314():
    os.chdir('poolSC_data/19scifi-RNA-seq')
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
GSE213957
'''

def PRJNA883380():
    # os.chdir('20CaRPool-seq/PRJNA883380/THP1/ADT/')
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

    os.chdir('20CaRPool-seq/PRJNA883380/THP1')
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

    mdata = MuData({'RNA': adata2, 'protein': adata1})
    mdata.write("mudata.h5mu")

    fun1(adata2)


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

            myrows, mycols = np.where(dat2 == 1)  #
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
    
    #PRJNA559094_fun(dirName='16TAP-seq/TAP_DIFFEX')
    #PRJNA559094_fun(dirName='16TAP-seq/WTX_DIFFEX')
    #PRJNA559094_fun(dirName='16TAP-seq/L1000')
    PRJNA559094_fun(dirName='16TAP-seq/REDESIGN')


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

            myrows, mycols = np.where(dat2 == 1)  ## assignment perturb
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
    
    PRJNA559094_fun2('PRJNA559094/SCREEN_chr8')
    PRJNA559094_fun2('PRJNA559094/SCREEN_chr11')
