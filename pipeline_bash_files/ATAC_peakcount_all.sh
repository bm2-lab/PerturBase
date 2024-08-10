
set -e
path=$1
species=$2

#### Change the Python and R environment to your own path
#### Change the Script to your own path

Python_env='/home/anaconda3/envs/pertpyV5'
R_env='/home/.conda/envs/Rversion4.2'
Script_base_path='/home/pipeline_folder'
source  /home/software/miniconda3/etc/profile.d/conda.sh 

# Signac preprocessing
conda activate $R_env
Rscript  ${Script_base_path}/scripts/ATAC/peakcount2tmp.R $path $species
conda activate $Python_env
python -u ${Script_base_path}/scripts/ATAC/tmp2h5ad.py $path

# Differentially accessible peaks
conda activate $R_env
Rscript  ${Script_base_path}/scripts/ATAC/ATAC_DAP.R $path $species
conda activate $Python_env
python -u ${Script_base_path}/scripts/ATAC/ATAC_filter_DAP.py $path $species

# Cluster Distribution and Edistance and correlation
python -u ${Script_base_path}/scripts/ATAC/ATAC_Edist_ClusterDistr.py $path $species

# Enrichment and plot
python -u ${Script_base_path}/scripts/ATAC/ATAC_Enrichment.py $path $species $R_env $Script_base_path
python -u ${Script_base_path}/scripts/ATAC/ATAC_Enrichment_plot.py $path $species $R_env $Script_base_path
