rm(list=ls())


outdir <- "~/GDAN//fpkm-uq/TCGA-HNSC/output/"




subtype <- "~/GDAN/fpkm-uq/TCGA-HNSC/subtype_2.txt"
geneannot_file <- "~/GDAN/scripts/geneannot.rds"

active_file_rds <- paste(outdir,"active.rds",sep="")
legacy_file_rds <- paste(outdir,"legacy.rds",sep="")
merged_file_rds <- paste(outdir,"merged_active_legacy.rds",sep="")
functions_file <- "~/GDAN/scripts/rsem_fpkm_functions.R"

cancer_var <- "HNSC"


subtype_var <- "Expression_subtype"
data.full <- readRDS(file=merged_file_rds)
data.ge <- subset(data.full,! row.names(data.full) %in% c("dataset",subtype_var,"entrezgene")) 

dataset <- data.full["dataset",]
data_df1 <- rbind(data.ge,cl=dataset)

dim(data_df1)
source(functions_file)
analysis(data_df1,"HNSC\n",outdir)




#subtype analysis

pat <- "-\\d\\d\\w-\\d\\d\\w-\\w\\w\\w\\w-\\d\\d$"

actdf <- readRDS(file=active_file_rds)
actdf <- add_subtype_row(actdf,subtype,subtype_var,"Barcode",pat)
rownames(actdf) <- gsub("\\.\\d+","",rownames(actdf), perl=TRUE)
actdf <- rbind(actdf,dataset=rep(38,dim(actdf)[2]))


legacy <- readRDS(file=legacy_file_rds)
en2hg <- readRDS(file=geneannot_file)

legacy$hgnc <- gsub("\\|\\d+$","",rownames(legacy),perl=TRUE)
legacy$entrezgene <- gsub("^[\\w\\-]+\\|","",rownames(legacy),perl=TRUE)
legacy <- legacy[!duplicated(legacy$entrezgene),]
legacy <- merge(legacy,en2hg,by.x='entrezgene',by.y='entrezgene')
legacy <- legacy[!duplicated(legacy$ensembl_gene_id),]
rownames(legacy) <- legacy$ensembl_gene_id
legacy <- rbind(legacy,dataset=rep(19,dim(legacy)[2]))
legacy <- add_subtype_row(legacy,subtype,subtype_var,"Barcode",pat)


mergeddf <- merge(legacy,actdf,by="row.names")
rownames(mergeddf) <- mergeddf$Row.names

dim(mergeddf)


mergeddf$Row.names <- NULL
mergeddf$hgnc <- NULL
mergeddf$ensembl_gene_id <- NULL
mergeddf$entrezgene <- NULL

fileout=paste(outdir,"merged_with_subtype.rds",sep="")
saveRDS(mergeddf,file=fileout)

m <- readRDS(file=paste(outdir,"merged_with_subtype.rds",sep=""))

plotname <- paste(cancer_var,subtype_var)
class_analysis(m,"dataset",subtype_var,plotname,outdir,cancer_var)

