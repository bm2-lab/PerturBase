suppressMessages(library("Seurat"))
suppressMessages(library("anndata"))
suppressMessages(library("reticulate"))
suppressMessages(library("R.utils"))
path_to_python <- "/home/wzt/anaconda3/bin/python"
use_python(path_to_python)

### rra模式，计算signature
### library("scMAGeCK") 不加载，直接用拷贝复制的函数，把permutation改成1000, 缩短计算的时间
source('/NFS_home/NFS_home_2/wzt/project/HC-CrisprScreen/myUtil.r')



Args = commandArgs(T)

dirName = Args[1]
species = Args[2]
setwd(dirName)

data <- read_h5ad("mixscape_hvg_filter_subset.h5ad")
data <- CreateSeuratObject(counts = t(as.matrix(data$layers['X_pert'])), meta.data = data$obs)  ### 用X_pert数据即矫正后的数据
adata <- ScaleData(data, do.center = F, do.scale = F)

bc_dox = read.table("barcode.txt", header = TRUE, as.is = TRUE, sep='\t')

if (species == 'hsa'){
    SIGNATURE = '/home/wzt/database/msigdb/hsa/h.all.v2023.1.Hs.symbols.gmt'
    scmageck_rra(BARCODE='barcode.txt', RDS=adata, SIGNATURE = SIGNATURE,  NEGCTRL="CTRL", KEEPTMP = FALSE, SAVEPATH = './scMageCK/')
} else{
    SIGNATURE = '/home/wzt/database/msigdb/mmu/mh.all.v2023.1.Mm.symbols.gmt'
    scmageck_rra(BARCODE='barcode.txt', RDS=adata, SIGNATURE = SIGNATURE,  NEGCTRL="CTRL", KEEPTMP = FALSE, SAVEPATH = './scMageCK/')
}



# ### 开始跑scMageck
# datList = read.table('/home/wzt/project/HC-CrisprScreen/results/dataInfo.tsv', sep='\t', header=T)
# for (index in 1:dim(datList)[1]){
#     dirName = datList$data[index]
#     species = datList$speics[index]
#     setwd(dirName)
#     print (dirName)
#     adata = readRDS('forSeurat.rds')
#     scmageck_rra(BARCODE='barcode.txt', RDS=adata, SIGNATURE = '/home/wzt/database/msigdb/hsa/test.gmt',  NEGCTRL="CTRL", KEEPTMP = FALSE)
#}


### 提供单个基因
###  rra_result1 <- scmageck_rra(BARCODE='barcode.txt', RDS=adata, LABEL ='scMageCK', GENE="MKI67", NEGCTRL="CTRL", KEEPTMP=FALSE, PATHWAY=FALSE, SAVEPATH = './scMageCK/')