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






enrichment <- read.csv(file_csv,header = FALSE)
if (length(rownames(enrichment))<=1){
  q(save = 'no')
}

print('pass the cricter')

enrichment <- read.csv(file_csv,header = T)




tmp_max = min(30,length(rownames(enrichment)))
enrichment = enrichment[1:tmp_max,]
enrichment$number = factor(1: length(rownames(enrichment)))


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

names(labels)<- 1:nrow(enrichment)






enrichment$log10_p.adj = log(enrichment$p.adjust,base = 10)

library("ggsci")
p = ggplot(data=enrichment, aes(x=number, y = Count, fill=log10_p.adj)) +
  geom_bar(stat="identity", width=0.8)  + 
  theme_test() + 
  scale_x_discrete(labels=labels) +
  xlab("KEGG term") + 
  theme(axis.text=element_text(face = "bold", color="gray50")) +
  labs(title = "The Most Enriched KEGG Terms")+ 
  theme(axis.text.x=element_text( hjust = 1))+
  coord_flip()+scale_y_continuous(expand = c(0,0))+ 
  scale_fill_gsea(reverse=TRUE)+
  theme(axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold",hjust = 0.5))

png(output_png,res=300,height = 1200,width = 2200)
p
dev.off()