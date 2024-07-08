library(ggplot2)
library(optparse)
option_list <- list(
  make_option(c("-f", "--file_csv"), type = "character", default = FALSE,
              action = "store", help = "your enrichment csv file"
  ),
  make_option(c("-o", "--output_png"), type = "character", default = FALSE,
              action = "store", help = "the output png file"
  ))
opt = parse_args(OptionParser(option_list = option_list, usage = "compareEnrichment for user"))
file_csv = opt$file_csv
output_png = opt$output_png

# file_csv = '/home/wzt/project/HC-CrisprScreen/poolSC_data/8ECCITE-seq/PRJNA641353/ECCITE/enrichment/GSFA/compareCluster_enrichment_GO.csv'


enrichment <- read.csv(file_csv,header = FALSE)
if (length(rownames(enrichment))<=1){
  q(save = 'no')
}
print('pass the cricter')
enrichment <- read.csv(file_csv,header = T)
enrichment$Perturb=enrichment$Cluster
tmp_perturb_max = 10

if(length(table(enrichment$Perturb))>10){
  tmp_perturb_max = 10
}else{
  tmp_perturb_max = length(table(enrichment$Perturb))
}

if (tmp_perturb_max >1){
  perturb_list = rownames(sort(table(enrichment$Perturb),decreasing = T)[1:tmp_perturb_max])
}else{
  perturb_list = c(rownames(table(enrichment$Perturb))[1])
}

enrichment = subset(enrichment,Perturb%in% perturb_list)



tmp_pathway_max = 30
if(length(table(enrichment$Description))>30){
  tmp_perturb_max = 30
}else{
  tmp_perturb_max = length(table(enrichment$Description))
}

if (tmp_pathway_max >1){
  description_list = rownames(sort(table(enrichment$Description),decreasing = T)[1:tmp_pathway_max])
}else{
  description_list = c(rownames(table(enrichment$Description))[1])
}

enrichment = subset(enrichment,Description%in% description_list)








shorten_names <- function(x, n_word=4, n_char=40){
  if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > 40))
  {
    if (nchar(x) > 40) x <- substr(x, 1, 40)
    x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)],
                     collapse=" "), "...", sep="")
    return(x)
  } 
  else
  {
    return(x)
  }
}

labels=(sapply(
  as.character(enrichment$Description),
  shorten_names))

enrichment$label = data.frame(labels)$labels


enrichment$log10_p.adj = log(enrichment$p.adjust,base = 10)
enrichment$"-log10_p.adj" = -log(enrichment$p.adjust,base = 10)

enrichment$GeneRatio <- as.character(enrichment$GeneRatio)


enrichment$GeneRatio <- sapply(enrichment$GeneRatio, function(x) eval(parse(text = x)))




label_list = rownames(sort(table(enrichment$label),decreasing = F))
enrichment$label = factor(enrichment$label,levels = label_list )


library("ggsci")
p = ggplot() +
  geom_point(data=enrichment, aes(y=Perturb, x=label,size = GeneRatio, color=`-log10_p.adj`))+
  theme_test() +
  xlab("KEGG term") + 
  theme(axis.text=element_text(face = "bold", color="gray50")) +
  labs(title = "The Most Enriched KEGG Terms")+ 
  theme(axis.text.x=element_text( angle = 90,hjust = 1))+
  coord_flip()+scale_fill_gsea(reverse=TRUE)+
  scale_color_gradient(low = "white", high = "red") +
  theme(axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold",hjust = 0.5))


png(output_png,res=300,height = 2000,width = 1800)
p
dev.off()







