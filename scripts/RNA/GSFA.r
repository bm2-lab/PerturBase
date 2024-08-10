suppressMessages(library("Seurat"))
suppressMessages(library("SCREE"))
suppressMessages(library("anndata"))
suppressMessages(library("reticulate"))
suppressMessages(library("sceptre"))
suppressMessages(library("Matrix"))
suppressMessages(library("data.table"))
suppressMessages(library("GSFA"))

Args = commandArgs(T)

species = Args[1]
path_to_python = Args[2] 
use_python(path_to_python)
myUtil_r = Args[3] 
# source('/NFS_home/NFS_home_2/wzt/project/HC-CrisprScreen/myUtil.r')
source(myUtil_r)




## Deviance residual transformation
data <- read_h5ad("mixscape_hvg_filter.h5ad")
scaled.gene_exp <- as.matrix(data$layers['X_pert']) ##


sample_names <- rownames(data)
rownames(scaled.gene_exp) <- sample_names
gene_names = colnames(data)
colnames(scaled.gene_exp) <- gene_names


G_mat = read.table('barcode_GSFA.txt', sep='\t', header = T, row.names = 1, check.names = F)
G_mat = t(as.matrix(G_mat))

fit <- fit_gsfa_multivar(Y = scaled.gene_exp, G = G_mat, 
                         K = 20,
                         prior_type = "mixture_normal", 
                         init.method = "svd",
                         verbose = T, return_samples = T)

saveRDS(fit, 'GSFA/fit.rds')
gibbs_PM <- fit$posterior_means
write.table(t(gibbs_PM$beta_pm), 'GSFA/beta_pm.tsv', sep = '\t', col.names = NA, quote = F, row.names = T)

### plot
if (TRUE){
  gibbs_PM <- fit$posterior_means   ##  gibbs_PM$beta_pm  effect size,     gibbs_PM$Gamma_pm   PIP  significant effects
  lfsr_mat <- fit$lfsr[, -ncol(fit$lfsr)]
  KO_names <- colnames(lfsr_mat)
  a = dotplot_beta_PIP(fit, target_names = KO_names, reorder_target = c(KO_names[KO_names!="CTRL"], "CTRL"))
  ggsave('GSFA/PerturbationEffects.pdf', plot = a, width = 7, height = 5)

}

fit = readRDS('GSFA/fit.rds')
gibbs_PM <- fit$posterior_means
effect_mat = gibbs_PM$W_pm %*% t(gibbs_PM$beta_pm[-nrow(gibbs_PM$beta_pm), ])
write.table(effect_mat, 'GSFA/effect.tsv', sep = '\t', col.names = NA, quote = F, row.names = T)


fit = readRDS('GSFA/fit.rds')
write.table(fit$lfsr, 'GSFA/lfsr.tsv', sep = '\t', col.names = NA, quote = F, row.names = T)

