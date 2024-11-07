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
import scanpy as sc


def myPool(func, mylist, processes):
    with Pool(processes) as pool:
        results = list(tqdm(pool.imap(func, mylist), total=len(mylist)))
    return results


def preGene(x,  ctrl = 'CTRL'):
    xs = x.split(',')
    xs1 = [i for i in xs if i !=  ctrl]
    if len(xs1) == 0:
        return 'CTRL'
    else:
        return ','.join(sorted(set(xs1)))


#### modified function of Mixscape to speed up
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

    if not isinstance(adata.X, np.ndarray):  ### convert to array
        adata.X = adata.X.toarray()
    adata.layers["X_pert"] = adata.X.copy()   #### convert to array to aovid out of memory

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

### basic quality control with scanpy standard pipeline
def preData(adata, filterNone=True, minNums = 30, shuffle=True, filterCom=False,  seed = 42, mtpercent = 10,  min_genes = 200):
    adata.var_names_make_unique()
    adata = adata[~adata.obs.index.duplicated()] ### drop cells with duplicated index
    if filterCom:
        tmp = adata.obs['gene'].apply(lambda x: True if ',' not in x else False);  adata = adata[tmp]
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
    if sum(adata.obs['pct_counts_mt'] < mtpercent) / adata.shape[0] <= 0.5: mtpercent += 5  ### set percentage to 15 if a lot of cell is filtered
    adata = adata[adata.obs.pct_counts_mt < mtpercent, :]
    filterMT = adata.shape[0]
    tmp = adata.obs['gene'].value_counts()   ### min counts 30
    genes = list(tmp[tmp >= minNums].index)
    if 'CTRL' not in genes: genes += ['CTRL']
    adata = adata[adata.obs['gene'].isin(genes), :]
    filterMinNums = adata.shape[0]
    
    sc.pp.normalize_total(adata, target_sum=1e4)   ##normalize
    sc.pp.log1p(adata)

    adata = adata[adata.obs.sort_values(by='gene').index,:]  ### sort 
    if shuffle:
        tmp = list(adata.obs.index)
        np.random.seed(seed); np.random.shuffle(tmp); adata = adata[tmp]
    return filterNoneNums, filterCells, filterMT, filterMinNums, adata

### preprocess gene id
def transID(adata, species):
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    current_dir = os.path.dirname(os.path.abspath(__file__))
    if species == 'Mus musculus':
        data_file_path = os.path.join(current_dir, 'id_mapping', 'ENSEMBL2SYMBOL_Mm.tsv')
        
    else:
        data_file_path = os.path.join(current_dir, 'id_mapping', 'ENSEMBL2SYMBOL_Hs.tsv')
    dat = pd.read_csv(data_file_path, sep='\t')
    if adata.var_names[0].startswith('ENS'):
        dat.set_index('ENSEMBL', inplace=True)
    else:
        dat.set_index('SYMBOL', inplace=True)
    dat = dat[~dat.index.duplicated()]
    adata.var = pd.merge(adata.var, dat, left_index=True, right_index=True,how='left')
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