suppressMessages(library("Seurat"))
suppressMessages(library("SCREE"))
suppressMessages(library("anndata"))
suppressMessages(library("reticulate"))

Args = commandArgs(T)
species = Args[1]


path_to_python = Args[2] 
use_python(path_to_python)
myUtil_r = Args[3] 
# source('/NFS_home/NFS_home_2/wzt/project/HC-CrisprScreen/myUtil.r')
source(myUtil_r)
### lr mode



score_cut = 0.2
pval_cut = 0.01

data <- read_h5ad("mixscape_hvg_filter.h5ad")
data <- CreateSeuratObject(counts = t(as.matrix(data$layers['X_pert'])), meta.data = data$obs) ##
adata <- ScaleData(data, do.center = F, do.scale = F)

bc_dox = read.table("barcode.txt", header = TRUE, as.is = TRUE, sep='\t', check.names = FALSE)

results <- improved_scmageck_lr_edit(BARCODE = bc_dox, 
  RDS = adata, 
  NEGCTRL = "CTRL", 
  SELECT_GENE = NULL, 
  PERMUTATION = 100, 
  SAVEPATH = "scMageCK/", 
  LAMBDA = 0.01,
  LABEL = 'scMageckLR',
  NTC_baseline = TRUE)

score <- results[[1]][, -1]
pval <- results[[2]][, -1]


if (TRUE){
DE_gene_plot(score = score, 
  pval = pval, 
  top = 25,
  project = "", 
  prefix = "scMageCK", 
  label = "", 
  pval_cut = pval_cut, 
  score_cut = score_cut, 
  sort_by = "number", 
  #y_break = c(80, 150), 
  width = 10, 
  height = 7)
}

#### DEGs

if (species == 'hsa'){
  database = 'org.Hs.eg.db'
}else{
  database = 'org.Mm.eg.db'
}

de_genes <- data.frame(non = rep(0, ncol(score)),
                        up = rep(0,ncol(score)), 
                        down = rep(0, ncol(score)))
rownames(de_genes) <- colnames(score)
for (i in colnames(score)) {
    scmageck <- data.frame(score = score[, i], pval = pval[, i])
    scmageck$Difference <- "non"
    scmageck$Difference[scmageck$score > score_cut & scmageck$pval < pval_cut] <- "up"
    scmageck$Difference[scmageck$score < -score_cut & scmageck$pval < pval_cut] <- "down"
    a <- plyr::count(scmageck$Difference)
    rownames(a) <- a$x
    for (j in a$x) {
        de_genes[i, j] <- a[j, 2]
    }
}

write.table(de_genes, 'scMageCK/deg.tsv', sep='\t', quote = F, col.names = NA, row.names = T)
mylist = c()
for (i in colnames(score)){
  geneNumes = de_genes[i, 'up'] + de_genes[i, 'down']
  if (geneNumes >= 20){
    mylist = c(mylist, i)
  }
}
########

if (TRUE){
heatmap(score = score, 
  pval = pval, 
  remove_neg = TRUE, 
  prefix = "scMageCK", 
  width = "auto", 
  height = 7)
}

