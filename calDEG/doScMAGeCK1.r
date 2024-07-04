library("Seurat")
library("anndata")
library("reticulate")
library("scMAGeCK")

path_to_python <- "/home/wzt/anaconda3/bin/python"
use_python(path_to_python)


### 准备数据  把h5ad数据转换为rds数据， seurat格式的数据
### 用矫正过后的数据, 后续由于rds数据量比较大，写入太满，先暂时在doScMAGeCK2脚本中用的时候再转换
datList = read.table('/home/wzt/project/HC-CrisprScreen/results/dataInfo.tsv', sep='\t', header=T)
for (dirName in datList$data){
    setwd(dirName)
    print (dirName)
    data <- read_h5ad("mixscape.h5ad")
    data <- CreateSeuratObject(counts = t(as.matrix(data$layers['X_pert'])), meta.data = data$obs)
    adata <- ScaleData(data, do.center = F, do.scale = F)
    saveRDS(adata, 'forScMageck.rds')
}



### 原始数据
# datList = read.table('/home/wzt/project/HC-CrisprScreen/results/dataInfo.tsv', sep='\t', header=T)
# for (dirName in datList$data){
#     setwd(dirName)
#     print (dirName)
#     data <- read_h5ad("forScMageck.h5ad")
#     data <- CreateSeuratObject(counts = t(as.matrix(data$to_df())), meta.data = data$obs)
#     adata <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
#     adata <- ScaleData(adata, do.center = F, do.scale = F)
#     saveRDS(adata, 'forScMageck.rds')
# }

#rra_result1 <- scmageck_rra(BARCODE='barcode.txt', RDS=adata, LABEL ='scMageCK', GENE="MKI67", NEGCTRL="CTRL", KEEPTMP=FALSE, PATHWAY=FALSE)





