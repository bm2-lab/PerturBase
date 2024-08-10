# path = '/home/sirm/graduation_design/snakemake_folder/Demo_data/ATAC/peakcount_file'
# species = 'Hs'
args <- commandArgs(trailingOnly = TRUE)
path = args[1]
species = args[2]


setwd(path)
h5ad_file = paste0('raw.h5ad')
output_folder = '.'
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
set.seed(1234)
options(warn=1)
library(SeuratDisk)
library(stringr)
library(Matrix)




# da_test_use_list = c('wilcox','t','LR')





if (species =='Hs'){
    genome = "hg38"
    blackregions = blacklist_hg38
    ensdb = EnsDb.Hsapiens.v86
}else{
    genome = "mm10"
    blackregions = blacklist_mm10
    ensdb = EnsDb.Mmusculus.v79
    
}

process_h5ad <- function(h5ad_file,genome){
    h5seurat_file = paste(str_sub(h5ad_file,end = -6),'h5seurat',sep = '.')
    Convert(h5ad_file, dest="h5seurat",
            assay = "peaks",
            overwrite=T)
    seurat_object <- LoadH5Seurat(h5seurat_file)
    seurat_object = UpdateSeuratObject(seurat_object)
    ##
    counts = seurat_object@assays$peaks@counts
    metadata = seurat_object@meta.data
    chrom_assay_fragment <- CreateChromatinAssay(
      counts = counts,
      sep = c("-", "-"),
      genome = genome,
      min.cells = 10,
      min.features = 200
    )
    ATAC_object <- CreateSeuratObject(
      counts = chrom_assay_fragment,
      assay = 'peaks', 
      project = 'ATAC',
      min.cells = 1,
      meta.data = metadata
    )
    ATAC_object = SetIdent(ATAC_object, value = "gene")
    return(ATAC_object)
    
    }
print('creating ATAC_data')
ATAC_object = process_h5ad(h5ad_file,genome)

ATAC_object = SetIdent(ATAC_object, value = "gene")
ATAC_object <- subset(
  x = ATAC_object,
  subset = nCount_peaks > 3000 &
    nCount_peaks < 30000 )
if ('None' %in% rownames(table(ATAC_object$gene))){
if (table(ATAC_object$gene)['None'] != 0){
    ATAC_object = subset(ATAC_object,idents=c('None'),invert=TRUE)
}
    }

if ('CTRL' %in% rownames(table(ATAC_object$gene))){
if (table(ATAC_object$gene)['CTRL'] != 0){    
    print('CTRL match')
}else{
    print('No CTRL in perturb gene!!! miss in the last')
}
}else{
    print('No CTRL in perturb gene!!! miss in the first')
}


ATAC_object<- RunTFIDF(ATAC_object)



print('tmp create')
metadata = ATAC_object@meta.data
metadata = metadata['gene']
write.table (metadata, file = paste0(output_folder,'/h5ad_metadata_tem.csv'), sep =",",  col.names =TRUE, quote =TRUE,row.names =TRUE)

# matrix_data = as.matrix(ATAC_object@assays$peaks@data)
# write.table (matrix_data, file = '/home/sirm/graduation_design/data_process/Signac/test/data/matrix_tem.csv', sep =",",  col.names =TRUE, quote =TRUE,row.names =TRUE)

featuredata = ATAC_object@assays$peaks@meta.features
write.table (featuredata, file = paste0(output_folder,'/h5ad_metafeature_tem.csv'), sep =",",  col.names =TRUE, quote =TRUE,row.names =TRUE)
sparse.gbm <- Matrix(ATAC_object@assays$peaks@data, sparse = T )
writeMM(obj = sparse.gbm, file=paste0(output_folder,'/h5ad_matrix_tem.mtx'))