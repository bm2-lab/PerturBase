library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
set.seed(1234)
options(warn=1)

#1.  covert fragment to peakcount object: ATAC_object

args <- commandArgs(trailingOnly = TRUE)
fragment_file = args[1]
cbc_perturbation_dict = args[2]
output_folder = args[3]


# fragment_file = '/home/sirm/graduation_design/HC-CrisprScreen/poolSC_data/07Perturb_ATAC/PRJNA714243/database/GM_LargeScreen_Rep4_fragment.tsv.gz'
# output_folder = './test_fragment'
# cbc_perturbation_dict = '/home/sirm/graduation_design/HC-CrisprScreen/poolSC_data/07Perturb_ATAC/PRJNA714243/database/GM_LargeScreen_Rep4_fragment_cbc_perturbation_dict.tsv'

species = 'Hs'

gene_features = c('GAPDH','B2M','TBP')
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

metadata = read.table(cbc_perturbation_dict, header=T, sep="\t",row.names = 1)


print('creating ATAC_data')
process_fragment <-function(fragment_file){
    # call peaks
    peaks =  CallPeaks(
      object = fragment_file,
    )
    # remove peaks on nonstandard chromosomes and in genomic blacklist regions
    peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
    # create fragments
    fragments <- CreateFragmentObject(fragment_file)
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



