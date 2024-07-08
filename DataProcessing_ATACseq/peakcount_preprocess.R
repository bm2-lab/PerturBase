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
ATAC_object = SetIdent(ATAC_object, value = "gene")

# removing None subset
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

ATAC_object = SetIdent(ATAC_object, value = "gene")

# Normalization and linear dimensional reduction
ATAC_object<- RunTFIDF(ATAC_object)
ATAC_object <- FindTopFeatures(ATAC_object, min.cutoff = 'q0')
ATAC_object<- RunSVD(ATAC_object)

# run scale and PCA
ATAC_object <- ScaleData(ATAC_object)
ATAC_object<- RunPCA(ATAC_object)

ATAC_object <- RunUMAP(object = ATAC_object, reduction = 'pca', dims = 2:30,  reduction.name = "umappca",reduction.key = "UMAPpca_")
ATAC_object <- RunUMAP(object = ATAC_object, reduction = 'lsi', dims = 2:30,  reduction.name = "umaplsi",reduction.key = "UMAPlsi_")



# gene activity matrix
gene.activities <- GeneActivity(ATAC_object)
# add the gene activity matrix to the Seurat object as a new assay and normalize it
ATAC_object[['RNA']] <- CreateAssayObject(counts = gene.activities)
ATAC_object <- NormalizeData(
  object = ATAC_object,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(ATAC_object$nCount_RNA)
)
write.table (gene.activities, file =paste(output_folder,"gene_activities_raw.csv",sep = "/"), sep =",",  col.names =TRUE, quote =TRUE,row.names = TRUE)