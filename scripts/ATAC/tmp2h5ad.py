import pandas as pd
import numpy as np
from glob import glob
import scanpy as sc
import subprocess
import os
import operator
import sys
path = sys.argv[1]
os.chdir(path)
output_folder = '.'
feature_file = os.path.join(output_folder,'h5ad_metafeature_tem.csv')
metadata_file = os.path.join(output_folder,'h5ad_metadata_tem.csv')
matrix_tmp = os.path.join(output_folder,'h5ad_matrix_tem.mtx')



obs_dict = pd.read_csv(metadata_file)
feature_dict = pd.read_csv(feature_file)
feature_dict.drop(feature_dict.columns,axis=1,inplace=True)
data = sc.read_mtx(matrix_tmp)
data = data.T
data.var = feature_dict
data.obs = obs_dict
data.write_h5ad('filtered.h5ad')