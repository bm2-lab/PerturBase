import pickle
import pandas as pd
import os
import subprocess
# import scanpy as sc
from glob import glob
import sys



def run_GO_single(GO_csv,GO_png1,GO_png2,R_env,Script_base_path):
    # cmd = '/home/wzt/.conda/envs/Rversion4.2/bin/Rscript /home/sirm/graduation_design/data_pre/single_data/enrichment_single_GO_barplot.R -f {} -o {}'.format(GO_csv,GO_png1)
    # subprocess.run(cmd,shell=True)
    cmd = '{}/bin/Rscript {}/scripts/enrichment_single_GO_dotplot.R -f {} -o {}'.format(R_env,Script_base_path,GO_csv,GO_png2)
    subprocess.run(cmd,shell=True)
def run_KEGG_single(GO_csv,GO_png1,GO_png2,R_env,Script_base_path):
    # cmd = '/home/wzt/.conda/envs/Rversion4.2/bin/Rscript /home/sirm/graduation_design/data_pre/single_data/enrichment_single_KEGG_barplot.R -f {} -o {}'.format(GO_csv,GO_png1)
    # subprocess.run(cmd,shell=True)
    cmd = '{}/bin/Rscript {}/scripts/enrichment_single_KEGG_dotplot.R -f {} -o {}'.format(R_env,Script_base_path,GO_csv,GO_png2)
    subprocess.run(cmd,shell=True)
def run_KEGG_all(GO_csv,GO_png,R_env,Script_base_path):
    cmd = '{}/bin/Rscript {}/scripts/enrichment_all_KEGG_dotplot.R -f {} -o {}'.format(R_env,Script_base_path,GO_csv,GO_png)
    subprocess.run(cmd,shell=True)
    
def run_GO_all(GO_csv,GO_png,R_env,Script_base_path):
    cmd = '{}/bin/Rscript {}/scripts/enrichment_all_GO_dotplot.R -f {} -o {}'.format(R_env,Script_base_path,GO_csv,GO_png)
    subprocess.run(cmd,shell=True)
def single_plot(x,R_env,Script_base_path):
    path = os.path.join(x,'single_gene')
    GO_csv_list = glob('{}/*_enrichment_GO.csv'.format(path))
    for GO_csv in GO_csv_list:
        print(GO_csv )
        GO_png1 = GO_csv.replace('_enrichment_GO.csv','_GO_barplot.png')
        GO_png2 = GO_csv.replace('_enrichment_GO.csv','_GO_dotplot.png')
        run_GO_single(GO_csv,GO_png1,GO_png2,R_env,Script_base_path)
    KEGG_csv_list = glob('{}/*_enrichment_KEGG.csv'.format(path))
    for KEGG_csv in KEGG_csv_list:
        print(KEGG_csv )
        GO_png1 = KEGG_csv.replace('_enrichment_KEGG.csv','_KEGG_barplot.png')
        GO_png2 = KEGG_csv.replace('_enrichment_KEGG.csv','_KEGG_dotplot.png')
        run_KEGG_single(KEGG_csv,GO_png1,GO_png2,R_env,Script_base_path)
def all_plot(x,R_env,Script_base_path):
    print('all_GO')
    GO_csv = os.path.join(x,'compareCluster_enrichment_GO.csv')
    GO_png = os.path.join(x,'GO_dotplot.png')
    run_GO_all(GO_csv,GO_png,R_env,Script_base_path)
    print('all_KEGG')
    KEGG_csv = os.path.join(x,'compareCluster_enrichment_KEGG.csv')
    KEGG_png = os.path.join(x,'KEGG_dotplot.png')
    run_KEGG_all(KEGG_csv,KEGG_png,R_env,Script_base_path)
def plot_method(x,R_env,Script_base_path):
    print(x)
    all_plot(x,R_env,Script_base_path)
    single_plot(x,R_env,Script_base_path)    



def plot_RNA_method(R_env,Script_base_path):
    method_list = ['Wilcoxon','sceptre','ttest','GSFA','scMageCK']
    for method in method_list:
        print(method)
        path = os.path.join('enrichment',method)
        if os.path.isdir(path):
            plot_method(path,R_env,Script_base_path)

def plot_ATAC_method(R_env,Script_base_path):
    method_list = ['LR','ttest','Wilcoxon']
    for method in method_list:
        print(method)
        path = os.path.join('enrichment',method)
        plot_method(path,R_env,Script_base_path)  
if __name__ == '__main__':
    path = sys.argv[1]
    species = sys.argv[2]
    R_env = sys.argv[3]
    Script_base_path = sys.argv[4]
    # path = '/home/sirm/graduation_design/snakemake_folder/Demo_data/ATAC/fragment_file'
    # species = 'Hs'
    os.chdir(path)
    plot_ATAC_method(R_env,Script_base_path)
    