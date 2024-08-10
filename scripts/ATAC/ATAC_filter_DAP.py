# path = '/home/sirm/graduation_design/snakemake_folder/Demo_data/ATAC/fragment_file'
# species = 'Hs'
import sys

path = sys.argv[1]
species = sys.argv[2]


import os
import subprocess 
from glob import glob
import pandas as pd

def make_enrichment(target_folder,avg_log2FC=1,p_val_adj=0.05,species = 'Hs'):
    DEG_result_file = target_folder.rstrip('/')+'/DEG_result.csv'
    DEG_target_file = target_folder.rstrip('/')+'/DEG.csv'
    tmp = target_folder.rstrip('/')+'/DEG_filter.csv'
    if os.path.isfile(tmp):
        cmd = 'mv {} {}'.format(tmp,DEG_target_file)
    else:
        print('process')
        data = pd.read_csv(DEG_result_file)
        data = data[data.apply(lambda x : True if float(x['p_val_adj'])<=p_val_adj and abs(x['avg_log2FC'])>=avg_log2FC  else False,axis=1)]
        data = data[['gene_name','perturb']]
        data.columns = ['DEG','Perturb']
        data.to_csv(DEG_target_file,index=None) 

           

method_list = ['Wilcoxon','ttest','LR']
folder = path
for method in method_list:
    folder_ss = folder.rstrip('/')+'/'+method
    print(folder_ss)
    if os.path.isdir(folder):
        make_enrichment(folder_ss,avg_log2FC=1,p_val_adj=0.05)