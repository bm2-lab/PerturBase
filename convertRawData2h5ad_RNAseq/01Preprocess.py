from myUtil import *
import scanpy as sc
import warnings
import anndata as ad
warnings.filterwarnings('ignore')
from scipy import sparse

basepath =  '/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb_seq/'

def fun1(adata, isWrite = True):
    if isWrite: adata.write_h5ad('raw.h5ad')
    print (adata.shape)
    adata = adata[~adata.obs['gene'].isin(['None', 'CTRL'])]

    tmp = [i.split(',') for i in adata.obs['gene']]
    tmp = [i for i in tmp if i != 'CTRL' and i != 'None']
    tmp = np.unique([i for j in tmp for i in j])
    print (tmp)
    print (len(tmp))

#mydata = scanpy.read_10x_mtx(path='./', prefix='GSM3564450_CSTARVE_', var_names='gene_ids', cache=True)

############ 01 perturb seq

'''
Massively parallel phenotyping of coding variants in cancer with Perturb-seq
'''

def PRJNA679579(gene='TP53'):
    os.chdir('/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb_seq/PRJNA679579/{}'.format(gene))
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

###  预处理文件。
# def f_PRJNA663571():
#     dirNames = ['sample' + str(i) for i in range(1, 19)]
#     dirNames = [os.path.join(basepath + '/PRJNA663571/', dirName) for dirName in dirNames]
#     for dirName in dirNames:
#         os.chdir(dirName)
#         tmp = glob.glob('*UMI*sample*.h5')[0]
#         cmd = 'ln -sf {} sample.h5'.format(tmp); subprocess.call(cmd, shell=True)
        
#         tmp = glob.glob('*sample*UMI.Counts.csv')[0]
#         cmd = 'ln -sf {}  UMI.Counts.csv'.format(tmp); subprocess.call(cmd, shell=True)
#         print (os.path.basename(dirName))
#         PRJNA663571()

def f_PRJNA663571():
    os.chdir('/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb_seq/PRJNA663571')
    dirNames = ['sample' + str(i) for i in range(1, 19)]
    dirNames = [os.path.join(basepath + '/PRJNA663571/', dirName) for dirName in dirNames]
    datList = [PRJNA663571(dirName) for dirName in dirNames]
    dat = ad.concat(datList, axis=0)
    os.chdir('/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb_seq/PRJNA663571')
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
    metaData['gene'] = metaData['gene'].apply(preGene)  ###
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
    dirName = '/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb-seq/PRJNA693896/CRISPRi_perturb_host'
    PRJNA693896(dirName)

'''
'''
def PRJNA549582():
    os.chdir('/NFS_home/NFS_home_2/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb_seq/PRJNA549582')
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
    os.chdir('/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb_seq/PRJNA554074')
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
    tmpPath = '/home/wzt/project/HC-CrisprScreen/poolSC_data/01Perturb_seq/PRJNA354362'
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
    dirName = '/home/wzt/project/HC-CrisprScreen/poolSC_data/6CROP-seq/PRJNA715235/Pre_inject'
    adata = PRJNA715235(dirName); fun1(adata)
    
    tmpPath = '/home/wzt/project/HC-CrisprScreen/poolSC_data/6CROP-seq/PRJNA715235/Driver_lib'
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
    

    tmpPath = '/home/wzt/project/HC-CrisprScreen/poolSC_data/6CROP-seq/PRJNA715235/Driver_sub_lib'
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
    os.chdir('/home/wzt/project/HC-CrisprScreen/poolSC_data/13arrayBased/96-plex_scRNA-seq')
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

    