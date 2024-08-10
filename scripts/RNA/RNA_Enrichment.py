import os
import subprocess
import pickle
import pandas as pd
import os
import subprocess
import scanpy as sc
from glob import glob
from multiprocessing import Pool
from tqdm import tqdm
import sys

def run_enrichment(csv_file,outputfile,species,perturbationClass,R_env,Script_base_path):
    print(csv_file)
    cmd = '{}/bin/Rscript {}/scripts/enrichment_full_png.R -f {} -s {} -t {} -o {} '.format(R_env,Script_base_path,csv_file,species,perturbationClass,outputfile)
    subprocess.run(cmd,shell=True)
    return(csv_file)

def move_KEGG_single(compare_KEGG,target_folder):
    data = pd.read_csv(compare_KEGG)
    if 'Perturb' not in data.columns:
        return 0
    perturb_list = set(data['Cluster'])
    for perturb in perturb_list:
        perturb = perturb.replace('/','')
        perturb_tmp = data[data['Cluster']==perturb]
        perturb_tmp = perturb_tmp.loc[:,['ID', 'Description', 'GeneRatio', 'BgRatio',
       'pvalue', 'p.adjust', 'qvalue', 'geneID', 'Count']]
        target_file = os.path.join(target_folder,'{}_enrichment_KEGG.csv'.format(perturb))
        perturb_tmp.to_csv(target_file,index=None)

def trans_pkl_csv(pkl_file,csv_file):
    with open(pkl_file,'rb') as f:
        data = pickle.load(f)
    data_list = []
    for i,j in data.items():
        if len(j)>0:
            x = pd.DataFrame(j,columns=['DEG']) 
            x['Perturb'] = i
            data_list.append(x) 
    if len(data_list) >0:      
        data = pd.concat(data_list,axis = 0)
        data.to_csv(csv_file,index=None)
    else:
        data = pd.DataFrame(data=None,columns=['DEG', 'Perturb'])
        data.to_csv(csv_file,index=None)

def run_RNA_KEGG():
    for method in ['Wilcoxon','sceptre','ttest','GSFA','scMageCK']:
        compare_KEGG = os.path.join('enrichment',method,'compareCluster_enrichment_KEGG.csv')
        target_folder =  os.path.join('enrichment',method,'single_gene')
        move_KEGG_single(compare_KEGG,target_folder)
        
def run_RNA_compare_cluster_GOKEGG(R_env,Script_base_path):
    perturbationClass = 'perturbationClass.tsv'
    cmd = 'rm -rf {}'.format('enrichment')
    subprocess.run(cmd,shell=True)
    cmd = 'mkdir {}'.format('enrichment')
    subprocess.run(cmd,shell=True)
    method_list = ['Wilcoxon','sceptre','ttest','GSFA','scMageCK']
    for method in method_list:
        print(method)
        pkl_file = method+'/deg.pkl'
        csv_file = method+'/DEG.csv'
        outputfile = 'enrichment/'+method
        if not os.path.isfile(pkl_file):
            print('{} do not have deg.pkl'.format(method))
            continue
        trans_pkl_csv(pkl_file,csv_file)
        cmd = 'mkdir {}'.format(outputfile)
        subprocess.run(cmd,shell=True)
        run_enrichment(csv_file,outputfile,species,perturbationClass,R_env,Script_base_path)

def run_enrichment_single_GO(csv_file,outputfile,R_env,Script_base_path):
    print(csv_file)
    cmd = '{}/bin/Rscript {}/scripts/enrichment_single_all.R -f {}   -o {}'.format(R_env,Script_base_path,csv_file,outputfile)
    subprocess.run(cmd,shell=True)
    return(csv_file)


def run_RNA_GO(R_env,Script_base_path):
    method_list = ['Wilcoxon','sceptre','ttest','GSFA','scMageCK']
    for method in method_list:
        print(method)
        csv_file = method+'/DEG.csv'
        outputfile = 'enrichment/'+method+'/single_gene'
        if not os.path.isfile(csv_file):
            print('{} do not have DEG.csv'.format(method))
            return
        cmd = 'mkdir {}'.format(outputfile)
        subprocess.run(cmd,shell=True)
        run_enrichment_single_GO(csv_file,outputfile,R_env,Script_base_path)
if __name__ == '__main__':
    path = sys.argv[1]
    species = sys.argv[2]
    R_env = sys.argv[3]
    Script_base_path = sys.argv[4]
    # path = '/home/sirm/graduation_design/pipeline_folder/Demo_data/RNA'
    # species = 'Hs'
    # R_env = '/home/wzt/.conda/envs/Rversion4.2'
    # Script_base_path = '/home/sirm/graduation_design/pipeline_folder'
    os.chdir(path)
    # compareCluster
    run_RNA_compare_cluster_GOKEGG(R_env,Script_base_path)
    # GO enrichment of single perturbation
    run_RNA_GO(R_env,Script_base_path)
    # KEGG enrichment of single perturbation
    run_RNA_KEGG()