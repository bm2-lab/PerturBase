import subprocess, os, sys, re, glob
from collections import defaultdict
import numpy as np, pandas as pd
from multiprocessing import Pool
from sklearn.utils import shuffle
from tqdm import tqdm
from multiprocessing import Pool
from scipy.sparse import csr_matrix, issparse, spmatrix
from anndata import AnnData
from scipy import sparse
from scanpy.tools._utils import _choose_representation
import anndata as ad
try:
    import scanpy as sc
except:
    pass

import matplotlib.pyplot as plt

#### 把h5ad转换为 seurat object

'''
library("Seurat")
library("anndata")
library(reticulate)

path_to_python <- "/home/wzt/anaconda3/bin/python"
use_python(path_to_python)

data <- read_h5ad("raw.h5ad")
data <- CreateSeuratObject(counts = t(as.matrix(data$X)), meta.data = data$obs)

'''


def myPool(func, mylist, processes):
    with Pool(processes) as pool:
        results = list(tqdm(pool.imap(func, mylist), total=len(mylist)))
    return results


def preGene(x,  ctrl = 'CTRL'):   ####  ctrl指定对照的标签
    xs = x.split(',')
    xs1 = [i for i in xs if i !=  ctrl]
    if len(xs1) == 0:
        return 'CTRL'
    else:
        return ','.join(sorted(set(xs1)))


def dowmLoadData(x):
    pass
    


def perturbation_signature(
    adata: AnnData,
    pert_key: str,
    control: str,
    split_by = None,
    n_neighbors: int = 20,
    use_rep = None,
    n_pcs = None,
    batch_size = None,
    copy: bool = False,
    **kwargs,
):
    if copy:
        adata = adata.copy()

    if not isinstance(adata.X, np.ndarray):  ### 修改成array格式
        adata.X = adata.X.toarray()
    adata.layers["X_pert"] = adata.X.copy()   #### 修改成array的格式，防止后续计算出现内存不足等错误

    control_mask = adata.obs[pert_key] == control

    if split_by is None:
        split_masks = [np.full(adata.n_obs, True, dtype=bool)]
    else:
        split_obs = adata.obs[split_by]
        cats = split_obs.unique()
        split_masks = [split_obs == cat for cat in cats]

    R = _choose_representation(adata, use_rep=use_rep, n_pcs=n_pcs)

    for split_mask in split_masks:
        control_mask_split = control_mask & split_mask

        R_split = R[split_mask]
        R_control = R[control_mask_split]

        from pynndescent import NNDescent  # saves a lot of import time

        eps = kwargs.pop("epsilon", 0.1)
        nn_index = NNDescent(R_control, **kwargs)
        indices, _ = nn_index.query(R_split, k=n_neighbors, epsilon=eps)

        X_control = np.expm1(adata.X[control_mask_split])

        n_split = split_mask.sum()
        n_control = X_control.shape[0]

        if batch_size is None:
            col_indices = np.ravel(indices)
            row_indices = np.repeat(np.arange(n_split), n_neighbors)

            neigh_matrix = csr_matrix(
                (np.ones_like(col_indices, dtype=np.float64), (row_indices, col_indices)),
                shape=(n_split, n_control),
            )
            neigh_matrix /= n_neighbors
            adata.layers["X_pert"][split_mask] -= np.log1p(neigh_matrix @ X_control)
        else:
            is_sparse = issparse(X_control)
            split_indices = np.where(split_mask)[0]
            for i in range(0, n_split, batch_size):
                size = min(i + batch_size, n_split)
                select = slice(i, size)

                batch = np.ravel(indices[select])
                split_batch = split_indices[select]
                size = size - i
                # sparse is very slow
                means_batch = X_control[batch]
                means_batch = means_batch.toarray() if is_sparse else means_batch
                means_batch = means_batch.reshape(size, n_neighbors, -1).mean(1)
                adata.layers["X_pert"][split_batch] -= np.log1p(means_batch)
    if copy:
        return adata



# def runMulti(adata, processes=10):
#     ctrlandNone = adata[adata.obs['gene'].isin(['None', 'CTRL'])].copy()
#     genelist = adata.obs['gene'].unique()
#     genelist = [i for i in genelist if i not in ['None', 'CTRL']]
#     adata_list = [adata[adata.obs['gene'].isin([i, 'CTRL'])].copy() for i in genelist]
#     result = myPool(low_multiple_filter_perturb_and_cell, adata_list, processes=processes)
#     adata1 = ad.concat(result)
#     adata2 = adata1[adata1.obs['gene'] != 'CTRL']
#     adata3 = ad.concat([adata2, ctrlandNone])
#     return adata3


### 批量跑mixscape的脚本
def preData(adata, filterNone=True, minNums = 30, shuffle=True, filterCom=False,  seed = 42, mtpercent = 10,  min_genes = 200):  #### 为了测试聚类，最好不要进行排序
    adata.var_names_make_unique()
    adata = adata[~adata.obs.index.duplicated()] ### 去除重复的细胞
    if filterCom:
        tmp = adata.obs['gene'].apply(lambda x: True if ',' not in x else False);  adata = adata[tmp] ### 过滤组合扰动
    if filterNone:
        adata = adata[adata.obs["gene"] != "None"]
    filterNoneNums = adata.shape[0]
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=3)
    filterCells = adata.shape[0]

    if np.any([True if i.startswith('mt-') else False for i in adata.var_names]):
        adata.var['mt'] = adata.var_names.str.startswith('mt-')
    else:
        adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    if sum(adata.obs['pct_counts_mt'] < mtpercent) / adata.shape[0] <= 0.5: mtpercent += 5
    adata = adata[adata.obs.pct_counts_mt < mtpercent, :]
    filterMT = adata.shape[0]
    tmp = adata.obs['gene'].value_counts()   ###去除细胞数量太少的扰动
    genes = list(tmp[tmp >= minNums].index)  ###
    if 'CTRL' not in genes: genes += ['CTRL']
    adata = adata[adata.obs['gene'].isin(genes), :]
    filterMinNums = adata.shape[0]
    
    sc.pp.normalize_total(adata, target_sum=1e4)   ## 预处理, normalize
    sc.pp.log1p(adata)

    adata = adata[adata.obs.sort_values(by='gene').index,:]  ### 排序
    if shuffle:
        tmp = list(adata.obs.index)
        np.random.seed(seed); np.random.shuffle(tmp); adata = adata[tmp]
    return filterNoneNums, filterCells, filterMT, filterMinNums, adata

### 转换ID
def transID(adata, species):
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    if species == 'Mus musculus':
        dat = pd.read_csv('/home/wzt/database/ENSEMBL2SYMBOL_Mm.tsv', sep='\t')
    else:
        dat = pd.read_csv('/home/wzt/database/ENSEMBL2SYMBOL_Hs.tsv', sep='\t')
    if adata.var_names[0].startswith('ENS'):
        dat.set_index('ENSEMBL', inplace=True)
    else:
        dat.set_index('SYMBOL', inplace=True)
    dat = dat[~dat.index.duplicated()]
    #df, adata_var = dat.align(adata.var, join="inner", axis=0) ### 会去除基因
    #adata = adata[:, adata_var.index]
    #adata.var = adata.var.merge(df, left_index=True, right_index=True)
    adata.var = pd.merge(adata.var, dat, left_index=True, right_index=True,how='left') ### 不去除基因
    if adata.var_names[0].startswith('ENS'):
        adata.var['ENSEMBL'] = adata.var.index
        adata.var.set_index('SYMBOL', inplace=True)
        adata = adata[:, ~adata.var_names.isna()]
    adata.var = adata.var[['ENSEMBL', 'ENTREZID']]
    return adata


def gen_mpl_labels( adata, groupby, exclude=(), ax=None, adjust_kwargs=None, text_kwargs=None):
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

    adjust_text(texts, **adjust_kwargs) # type: ignore



def getStrongPerturb(Nums = 25):
    perturbationClass = pd.read_csv('perturbationClass.tsv', sep='\t')
    strongP = list(perturbationClass['perturbation'][:Nums])
    if 'CTRL' not in strongP: strongP.append('CTRL')
    return  strongP