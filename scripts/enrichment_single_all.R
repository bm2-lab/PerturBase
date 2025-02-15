library(optparse)
library(org.Hs.eg.db) 
library(org.Mm.eg.db)


library(clusterProfiler)

library(enrichplot)
library(ggplot2)
library(openxlsx)
library(ggplot2)
library(stringr)
library(enrichplot)



option_list <- list(
  make_option(c("-f", "--DEG_list_file"), type = "character", default = FALSE,
              action = "store", help = "your DEG_list_file"
  ),
  make_option(c("-s", "--species"), type = "character", default = 'Hs',
              action = "store", help = "species"
  ),
  make_option(c("-o", "--output_folder"), type = "character", default = FALSE,
              action = "store", help = "output_folder"
  ),
  make_option(c("-m", "--method"), type = "character", default = "GO,KEGG",
              action = "store", help = "enrichment method"
  ),  
  make_option(c("-A", "--pAdjustMethod"), type = "character", default = "BH",
              action = "store", help = "pAdjustMethod"
  ),  
  make_option(c("-p", "--pvalueCutoff"), type = "double", default = 0.01,
              action = "store", help = "pvalueCutoff"
  ),  
  make_option(c("-q", "--qvalueCutoff"), type = "double", default = 0.05,
              action = "store", help = "qvalueCutoff"
  ),  
  make_option(c("-g", "--query_gene"), type = "character", default = 'all',
              action = "store", help = "query_gene"
  )
)

opt = parse_args(OptionParser(option_list = option_list, usage = "compareEnrichment for user"))

DEG_list_file = opt$DEG_list_file
species = opt$species
output_folder = opt$output_folder
method = unlist(strsplit(opt$method, ","))
pAdjustMethod = opt$pAdjustMethod
pvalueCutoff = opt$pvalueCutoff
qvalueCutoff = opt$qvalueCutoff
DEG_list = read.table(DEG_list_file,sep=',',head=TRUE)
query_gene = opt$query_gene
colnames(DEG_list) <- c("SYMBOL", "Perturb")
if (length(rownames(DEG_list))==0){
    print('the gene is 0')
    q("no")
}

if (query_gene=='all'){
    query_gene = rownames(table(DEG_list$Perturb))
}else{
    query_gene = unlist(strsplit(opt$query_gene, ","))
}


DEG_list = subset(DEG_list,Perturb %in% query_gene)



if (species == 'Hs'){
        OrgDb = org.Hs.eg.db
        organism = "hsa"
    }else {
        OrgDb = org.Mm.eg.db
        organism = "mmu"
    }

    #print(DEG_list)
SYMBOL_to_ENTREZID =  bitr(DEG_list[['SYMBOL']], 
        fromType="SYMBOL", 
        toType="ENTREZID",  
        OrgDb=OrgDb) 
    #print(SYMBOL_to_ENTREZID)
missing_rate = (length(unique(DEG_list[[1]])) - length(SYMBOL_to_ENTREZID[[1]]))/length(unique(DEG_list[[1]]))
if ( missing_rate> 0.30){
    print('too many genes missing ENTREZID,please check your gene symbol or species')
}

DEG_list = merge(DEG_list,SYMBOL_to_ENTREZID)


print(length(query_gene))


for (gene_per in query_gene){
    print(gene_per)
    DEG_list_tmp = subset(DEG_list,Perturb == gene_per)
    if (length(rownames(DEG_list_tmp))==0){
    print(paste(gene_per,'the gene is 0'))
    next
    }
    gene_per = gsub("/","",gene_per)
    if ('GO' %in% method){
    GO<-enrichGO( DEG_list_tmp$ENTREZID,
              OrgDb = OrgDb,
              keyType = "ENTREZID",
              ont = "ALL",
                pAdjustMethod=pAdjustMethod,
              pvalueCutoff = pvalueCutoff,
              qvalueCutoff = qvalueCutoff
              )
    write.table (GO, file =paste(output_folder,paste(gene_per,"enrichment_GO.csv",sep='_'),sep = "/"), sep =",",  col.names =TRUE, quote =TRUE,row.names = FALSE)    
    
    }
}

    
    
    
