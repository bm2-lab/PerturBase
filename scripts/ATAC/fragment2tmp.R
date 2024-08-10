# path = '/home/sirm/graduation_design/snakemake_folder/Demo_data/ATAC/fragment_file'
# species = 'Hs'
args <- commandArgs(trailingOnly = TRUE)
path = args[1]
species = args[2]

setwd(path)
fragment_file = paste0('fragment.tsv.gz')
cbc_perturbation_dict = paste0('cbc_perturbation_dict.tsv')
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
library(Matrix)



# da_test_use_list = c('wilcox','t','LR')





if (species == "Hs") {
    genome <- "hg38"
    blackregions <- blacklist_hg38
    ensdb <- EnsDb.Hsapiens.v86
} else if (species == "Mm") {
    genome <- "mm10"
    blackregions <- blacklist_mm10
    ensdb <- EnsDb.Mmusculus.v79
} else {
    stop("species variable contains an unrecognized value.")
}

metadata = read.table(cbc_perturbation_dict, header=T, sep="\t",row.names = 1)


print('creating ATAC_data')
process_fragment <-function(fragment_file){
    # call peaks
    peaks =  CallPeaks(
      object = fragment_file,
      macs2.path = '/home/sirm/.conda/envs/Rversion4.2/bin/macs2'
    )
    # remove peaks on nonstandard chromosomes and in genomic blacklist regions
    peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
    # create fragments
    fragments <- CreateFragmentObject(fragment_file) # need indexed fragment file .tbi
    peaks_matrix = FeatureMatrix(
    fragments = fragments,
    features = peaks,
    cells = NULL,
    process_n = 2000,
    sep = c("-", "-"),
    verbose = TRUE
    )
    return(peaks_matrix)
}


counts = process_fragment(fragment_file)
chrom_assay_fragment <- CreateChromatinAssay(
  counts = counts,
  sep = c("-", "-"),
  genome = genome,
  fragments = fragment_file,
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
# change the level to gene
ATAC_object = SetIdent(ATAC_object, value = "gene")
print('QC')

## FRIP
total_fragments <- CountFragments(fragment_file)
rownames(total_fragments) <- total_fragments$CB
ATAC_object$fragments <- total_fragments[colnames(ATAC_object), "frequency_count"]
ATAC_object <- FRiP(
  object = ATAC_object,
  assay = 'peaks',
  total.fragments = 'fragments'
)

## blacklist_fraction
ATAC_object$blacklist_fraction <- FractionCountsInRegion(
  object = ATAC_object, 
  assay = 'peaks',
  regions = blackregions
)

## annotation
annotations <- GetGRangesFromEnsDb(ensdb = ensdb)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- genome
Annotation(ATAC_object) <- annotations

# compute nucleosome signal score per cell
ATAC_object = NucleosomeSignal(object = ATAC_object,assay = "peaks")

# compute TSS enrichment score per cell
ATAC_object <- TSSEnrichment(object = ATAC_object, fast = FALSE)

print('filtering')
ATAC_object <- subset(
  x = ATAC_object,
  subset = nCount_peaks > 3000 &
    nCount_peaks < 30000 &
    FRiP > 0.15 &
    blacklist_fraction < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 3
)


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