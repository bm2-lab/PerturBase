from myUtil import *
import scanpy as sc #type: ignore
import warnings
warnings.filterwarnings('ignore')
import pertpy as pt #type: ignore
import os, sys


#### basic quality control and denosing with Mixscape

def fun1(adata, isWrite = True):
    if isWrite: adata.write_h5ad('raw.h5ad')
    a = adata.shape
    adata = adata[~adata.obs['gene'].isin(['None', 'CTRL'])]

    tmp = [i.split(',') for i in adata.obs['gene']]
    tmp = [i for i in tmp if i != 'CTRL' and i != 'None']
    tmp = np.unique([i for j in tmp for i in j])
    return a, len(tmp)

def classGene(adata, threshold = 0.2):  ### classify a perturbation to weak or strong perturbation accroding to the ratio of non-perturbed cells
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
    adata1 = adata[~((adata.obs['mixscape_class_global'] == 'NP') & (adata.obs['pertclass'] == 'strong'))]  #### filter non-perturbed cells
    adata2 = adata[((adata.obs['mixscape_class_global'] == 'NP') & (adata.obs['pertclass'] == 'strong'))]  ### non-perturbed cells
    return strongP, weakP, adata1, adata2


def doMixScape1(adata, filterNone=True, minNums = 30, shuffle=True, filterCom = False, mtpercent=10, min_genes=200, keepall = False):  #### denosing using Mixscape
    if adata.shape[1] <= 1500: min_genes = 0
    if keepall: minNums = 1  ### do not filter cells
    filterNoneNums, filterCells, filterMT, filterMinNums, adata = preData(adata, filterNone, minNums, shuffle, filterCom,  seed=42, mtpercent=mtpercent, min_genes=min_genes)  ### basic quality control using scanpy standard pipeline
    if len(adata.obs['gene'].unique()) <=1:
        return filterNoneNums, filterCells, filterMT, filterMinNums, adata, 0
    if adata.shape[1] <= 1500:  ### min gene 1500
        sc.pp.highly_variable_genes(adata, subset=True, n_top_genes=adata.shape[1])
    else:
        sc.pp.highly_variable_genes(adata, subset=True, n_top_genes=4000)
        # sc.pp.highly_variable_genes(adata, subset=True, n_top_genes=2000)
    sc.pp.pca(adata)
    if not isinstance(adata.X, np.ndarray):
        adata.X = adata.X.toarray()
    adata = doPCA1(adata)
    adata.write_h5ad('RawNormal.h5ad')
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


def doPCA(adata):  ### pca with denoising data
    adata.layers['normal_counts'] = adata.X.copy()
    adata.X = adata.layers['X_pert']
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata)
    return adata

def doPCA1(adata): ### pca with default expression 
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

def f_doMixScape(species, keepall= False):
    adata = sc.read_h5ad('raw.h5ad')
    rawTranscriptNums = adata.shape[1]
    adata = transID(adata, species)
    TranscriptNums = adata.shape[1]
    filterNoneNums, filterCells, filterMT, filterMinNums, adata1, allGeneAffilter = doMixScape1(adata, mtpercent=10, keepall=keepall)
    if len(adata1.obs['gene'].unique()) <= 1:
        RawCellNums, RawGeneNums, CellNums, GeneNums = fun1(adata, False)[0][0], fun1(adata, False)[1], 0, 0
        return RawCellNums, RawGeneNums, CellNums, GeneNums, 0, 0, 0, 0, filterNoneNums, filterCells, filterMT, filterMinNums, rawTranscriptNums, TranscriptNums
    strongP, weakP, adata_pert, adata_pert1 = classGene(adata1, threshold=.2)  ### 对 filter non-perturbed cells
    adata_filter = adata[adata_pert.obs_names, allGeneAffilter]   ### raw data after basic quality control
    adata1.write_h5ad('mixscape_hvg.h5ad')

    adata_pert = doPCA(adata_pert)
    adata_pert.write_h5ad('mixscape_hvg_filter.h5ad')  ### data after all quality control
    adata_filter.write_h5ad('filterRaw.h5ad')  ### raw data after filtering non-perturbed cells

    adata2 = doMixScape2(adata, mtpercent= 10, keepall=keepall)   ### keep all genes instead of hvg, time consuming
    adata2.write_h5ad('mixscape.h5ad')
    adata3 = adata2[adata_pert.obs_names, ]
    adata3.write_h5ad('mixscape_filter.h5ad')
    
    RawCellNums, RawGeneNums, CellNums, GeneNums = fun1(adata, False)[0][0], fun1(adata, False)[1], fun1(adata1, False)[0][0], fun1(adata1, False)[1]
    return RawCellNums, RawGeneNums, CellNums, GeneNums, len(strongP), len(weakP), adata_pert.shape[0], adata_pert1.shape[0], filterNoneNums, filterCells, filterMT, filterMinNums, rawTranscriptNums, TranscriptNums



def ff_doMixScape(species='Homo sapiens', keepall=False):  #Mus musculus
    fileout = 'dataInfo.tsv'
    with open(fileout, 'w') as fout:
        fout.write('data\tspecies\tRawCellNums\tRawGeneNums\tCellNums\tGeneNums\tstrongP\tweakP\tCellNums_strongP\tCellNums_weakP\tfilterNoneNums\tfilterCells\tfilterMT\tfilterMinNums\trawTranscriptNums\tTranscriptNums\n')
        result =  f_doMixScape(species, keepall)
        result = [str(i) for i in result]
        fout.write('{}\t{}\n'.format(species, '\t'.join(result)))


if __name__ == '__main__':
    print(sys.argv)
    path = sys.argv[1]
    os.chdir(path)
    species = sys.argv[2]
    if species == 'Hs':
        species = 'Homo sapiens'
    else:
        species = 'Mus musculus'
    ff_doMixScape(species)
