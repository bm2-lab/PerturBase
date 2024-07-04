suppressMessages(library("Seurat"))
suppressMessages(library("anndata"))
suppressMessages(library("SCREE"))
suppressMessages(library("reticulate"))
suppressMessages(library("R.utils"))
suppressMessages(library(gridExtra))

path_to_python <- "/home/wzt/anaconda3/bin/python"
use_python(path_to_python)

source('/NFS_home/NFS_home_2/wzt/project/HC-CrisprScreen/myUtil.r')

data <- read_h5ad("rawHVG.h5ad")
data <- CreateSeuratObject(counts = t(as.matrix(data$to_df())), meta.data = data$obs)

sg_lib = 'barcode.txt'
mtx = data
sgRNA_quality_plot(mtx, sg_lib)

