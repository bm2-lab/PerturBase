# bash RNA_all.sh /home/pipeline_folder/Demo_data/RNA Hs 1>RNA_all.stdout 2>RNA_all.stderr 
set -e
path=$1
species=$2

#### Change the Python and R environment to your own path
#### Change the Script to your own path
#### Demo:
# Python_env='/home/anaconda3/envs/pertpyV5'
# R_env='/home/.conda/envs/Rversion4.2'
# Script_base_path='/home/pipeline_folder'
# msigdb_Signature_path='/home/Github_source_code/PerturBase/MsigDB_signature'
# source  /home/software/miniconda3/etc/profile.d/conda.sh
# chmod 777 /home/Github_source_code/PerturBase/RRA # make it executable
Python_env='/path/to/your/Python_env'
R_env='/path/to/your/R_env'
Script_base_path='/path/to/your/pipeline_folder'
msigdb_Signature_path='/path/to/your/MsigDB_signature'
chmod 777 /path/to/your/pipeline_folder/RRA 
source /path/to/your/miniconda3/etc/profile.d/conda.sh


conda activate $Python_env
# denoising
python -u ${Script_base_path}/scripts/RNA/denoising.py $path $species
# NP SP ratio and mixscape
python -u ${Script_base_path}/scripts/RNA/Mixscape.py $path $species
# analysis 
python -u ${Script_base_path}/scripts/RNA/analysis.py $path $species $R_env $Script_base_path $Python_env $msigdb_Signature_path
# Enrichment and plot
python -u ${Script_base_path}/scripts/RNA/RNA_Enrichment.py $path $species $R_env $Script_base_path
python -u ${Script_base_path}/scripts/RNA/RNA_Enrichment_plot.py $path $species $R_env $Script_base_path

