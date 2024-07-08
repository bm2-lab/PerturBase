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

args <- commandArgs(trailingOnly = TRUE)
h5ad_file = args[1]
output_folder = args[2]

species = 'Hs'

# gene_features = c('GAPDH','B2M','TBP')
da_test_use_list = c('wilcox','t','LR')
da_latent_vars = 'nCount_peaks'




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