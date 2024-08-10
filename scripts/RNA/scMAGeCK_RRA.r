suppressMessages(library("Seurat"))
suppressMessages(library("anndata"))
suppressMessages(library("reticulate"))
suppressMessages(library("R.utils"))
Args = commandArgs(T)
species = Args[1]

path_to_python = Args[2] 
use_python(path_to_python)
myUtil_r = Args[3] 
msigdb_Signature_path = Args[4] 
Script_base_path = Args[5] 
### rra mode

### calculate signature
source(myUtil_r)
print(species)




data <- read_h5ad("mixscape_hvg_filter.h5ad")
data <- CreateSeuratObject(counts = t(as.matrix(data$layers['X_pert'])), meta.data = data$obs)
adata <- ScaleData(data, do.center = F, do.scale = F)

bc_dox = read.table("barcode.txt", header = TRUE, as.is = TRUE, sep='\t')

RRAPATH = paste0(Script_base_path,'/executable/RRA')
if (species == 'hsa'){
    SIGNATURE = paste0(msigdb_Signature_path,'/hsa/h.all.v2023.2.Hs.symbols.gmt')
    scmageck_rra(BARCODE='barcode.txt', RDS=adata, RRAPATH = RRAPATH ,SIGNATURE = SIGNATURE,  NEGCTRL="CTRL", KEEPTMP = FALSE, SAVEPATH = './scMageCK/')
} else{
    SIGNATURE = paste0(msigdb_Signature_path,'/mmu/mh.all.v2023.2.Mm.symbols.gmt')
    scmageck_rra(BARCODE='barcode.txt', RDS=adata, RRAPATH = RRAPATH , SIGNATURE = SIGNATURE,  NEGCTRL="CTRL", KEEPTMP = FALSE, SAVEPATH = './scMageCK/')
}