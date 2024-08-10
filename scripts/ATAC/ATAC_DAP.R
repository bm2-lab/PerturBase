# path = '/home/sirm/graduation_design/snakemake_folder/Demo_data/ATAC/fragment_file'
# species = 'Hs'
args <- commandArgs(trailingOnly = TRUE)
path = args[1]
species = args[2]


setwd(path)
h5ad_file = paste0('filtered.h5ad')
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

metadata_file <- paste0('h5ad_metadata_tem.csv')
metafeature_file <- paste0('h5ad_metafeature_tem.csv')
matrix_file <- paste0( 'h5ad_matrix_tem.mtx')


sparse_matrix <- readMM(file = matrix_file)


cell_metadata <- read.csv(metadata_file, row.names = 1)


feature_metadata <- read.csv(metafeature_file, row.names = 1)


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

colnames(sparse_matrix) <- rownames(cell_metadata)
rownames(sparse_matrix) <- rownames(feature_metadata)


seurat_object <- CreateSeuratObject(counts = sparse_matrix, meta.data = cell_metadata)


seurat_object[["peaks"]] <- CreateAssayObject(counts = sparse_matrix, meta.features = feature_metadata)
counts = seurat_object@assays$peaks@counts
metadata = seurat_object@meta.data
chrom_assay_fragment <- CreateChromatinAssay(
  counts = counts,
  sep = c("-", "-"),
  genome = genome,
  min.cells = 0,
  min.features = 0
)
ATAC_object <- CreateSeuratObject(
  counts = chrom_assay_fragment,
  assay = 'peaks', 
  project = 'ATAC',
  min.cells = 0,
  meta.data = metadata
)
ATAC_object = SetIdent(ATAC_object, value = "gene")
annotations <- GetGRangesFromEnsDb(ensdb = ensdb)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- genome
Annotation(ATAC_object) <- annotations
# da for all perturb
ATAC_object = SetIdent(ATAC_object, value = "gene")
DefaultAssay(ATAC_object) <- 'peaks'
number_table = table(Idents(ATAC_object))
perturb_list = rownames(table(Idents(ATAC_object)))

# all the perturb with less than 30 cells must be filtered
tmp_perturb = c()
for (perturb in perturb_list){
    if (number_table[perturb]>30){
        tmp_perturb = append(tmp_perturb,perturb)
    }
}

perturb_list = tmp_perturb
da_test_use_list = c('Wilcoxon','ttest','LR')
da_test_mapping <- c(Wilcoxon = 'wilcox', ttest = 't', LR = 'LR')
da_latent_vars = 'nCount_peaks'
da_result_check=1

if ('CTRL' %in% perturb_list && length(perturb_list)>1){
    print('DA')
    for (da_test_use in da_test_use_list){
      print(da_test_use)  
if (!file.exists(paste(output_folder,da_test_use,sep = "/"))) {
  dir.create(paste(output_folder,da_test_use,sep = "/"), recursive = TRUE)
}
    #create data frame with 0 rows and 3 columns
    da_result <- data.frame(matrix(ncol = 7, nrow = 0))
    #provide column names
    colnames(da_result) <- c('p_val', 'avg_log2FC', 'pct.1','pct.2','p_val_adj','perturb','DA')

    for (perturb in perturb_list){
        if (perturb != 'CTRL'){
            print(perturb)
            da_peaks <- FindMarkers(
            object = ATAC_object,
            ident.1 = perturb,
            ident.2 = "CTRL",
            test.use = da_test_mapping[da_test_use],
            latent.vars = da_latent_vars
            )
            da_peaks$perturb = perturb
            da_peaks$DA = rownames(da_peaks)
            da_result = rbind(da_result,da_peaks)

        }    
    }
    

    # da_peaks[da_peaks$avg_log2FC > log2FC_thred&& da_result$p_val_adj < p_thred, ]
    # differient Acessibility
    diff_da = da_result

    write.table (diff_da, file =paste(output_folder,da_test_use,"Diff_accessibility.csv",sep = "/"), sep =",",  col.names =TRUE, quote =TRUE,row.names = FALSE)



    # map peaks to gene 
    #create data frame with 0 rows and 3 columns
    DEG_result <- data.frame(matrix(ncol = 8, nrow = 0))
    #provide column names
    colnames(DEG_result) <- c('tx_id', 'gene_name', 'gene_id','gene_biotype','type','closest_region','query_region','distance')

    print('DEG')

    for (perturb in rownames(table(diff_da$perturb))){
        print(perturb)
        sub_diff_da = diff_da[diff_da$perturb==perturb,]
    #     # print(length(rownames(sub_diff_da)))
    #     #rownames(sub_diff_da) = sub_diff_da$DA
        DEG = ClosestFeature(ATAC_object, regions = sub_diff_da$DA)
        DEG_result = rbind(DEG_result,DEG)

    }

    DEG_file = cbind(DEG_result,diff_da)
    write.table (DEG_file, file =paste(output_folder,da_test_use,"DEG_result.csv",sep = "/"), sep =",",  col.names =TRUE, quote =TRUE,row.names =FALSE)
    }
}else{
    da_result_check=0
    print('no CTRL or no perturb gene pass the check')
}