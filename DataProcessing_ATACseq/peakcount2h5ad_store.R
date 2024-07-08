library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
set.seed(1234)
options(warn=1)
library(Matrix)


args <- commandArgs(trailingOnly = TRUE)
fragment_file = args[1]
cbc_perturbation_dict = args[2]
species = 'Hs'

metadata = read.table(cbc_perturbation_dict, header=T, sep="\t",row.names = 1)

if (species =='Hs'){
    genome = "hg38"
    blackregions = blacklist_hg38
    ensdb = EnsDb.Hsapiens.v86
}else{
    genome = "mm10"
    blackregions = blacklist_mm10
    ensdb = EnsDb.Mmusculus.v79
    
}





print('tmp file create')
metadata = ATAC_object@meta.data
metadata = metadata['gene']
write.table (metadata, file = '/home/sirm/graduation_design/data_process/Signac/test/data/h5ad_metadata_tem.csv', sep =",",  col.names =TRUE, quote =TRUE,row.names =TRUE)

# matrix_data = as.matrix(ATAC_object@assays$peaks@data)
# write.table (matrix_data, file = '/home/sirm/graduation_design/data_process/Signac/test/data/matrix_tem.csv', sep =",",  col.names =TRUE, quote =TRUE,row.names =TRUE)

featuredata = ATAC_object@assays$peaks@meta.features
write.table (featuredata, file = '/home/sirm/graduation_design/data_process/Signac/test/data/h5ad_metafeature_tem.csv', sep =",",  col.names =TRUE, quote =TRUE,row.names =TRUE)
sparse.gbm <- Matrix(ATAC_object@assays$peaks@data, sparse = T )
writeMM(obj = sparse.gbm, file="/home/sirm/graduation_design/data_process/Signac/test/data/h5ad_matrix_tem.mtx")

## following we will use scanpy to read 3 tmp files to crate h5ad file.