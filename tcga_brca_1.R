rm(list=ls())

legacy_dir <- "~/GDAN/fpkm-uq/TCGA-BRCA/legacy/data/" 
active_dir <- "~/GDAN/fpkm-uq/TCGA-BRCA/active/data/"
geneannot_file <- "~/GDAN/scripts/geneannot.rds"
outdir <- "~/GDAN//fpkm-uq/TCGA-BRCA/output/"
functions_file <- "~/GDAN/scripts/rsem_fpkm_functions.R"

source(functions_file)



file_pat <- "TCGA-*"

ftm <- proc.time()

ptm <- proc.time()
active <- create_df_active_v2 (active_dir,file_pat)
fileout <- paste(outdir,"active.rds",sep="")
saveRDS(active,file=fileout)
proc.time() - ptm

ptm <- proc.time()
legacy <- create_df_active_v2 (legacy_dir,file_pat)
fileout <- paste(outdir,"legacy.rds",sep="")
saveRDS(legacy,file=fileout)
proc.time() - ptm


dim(active)
dim(legacy)


#load(en2hg)
en2hg <- readRDS(file=geneannot_file)

legacy$hgnc <- gsub("\\|\\d+$","",rownames(legacy),perl=TRUE)
legacy$entrezgene <- gsub("^[\\w\\-]+\\|","",rownames(legacy),perl=TRUE)
legacy <- legacy[!duplicated(legacy$entrezgene),]
legacy <- merge(legacy,en2hg,by.x='entrezgene',by.y='entrezgene')
legacy <- legacy[!duplicated(legacy$ensembl_gene_id),]
rownames(legacy) <- legacy$ensembl_gene_id
legacy <- rbind(legacy,dataset=rep(19,dim(legacy)[2]))

rownames(active) <- gsub("\\.\\d+","",rownames(active), perl=TRUE)
active <- rbind(active,dataset=rep(38,dim(active)[2]))

mergeddf <- merge(legacy,active,by="row.names")
rownames(mergeddf) <- mergeddf$Row.names

dim(mergeddf)


mergeddf$Row.names <- NULL
mergeddf$hgnc <- NULL
mergeddf$ensembl_gene_id <- NULL
mergeddf$entrezgene <- NULL
mergeddf$hgnc_symbol <- NULL



fileout <- paste(outdir,"merged_active_legacy.rds",sep="")
#saveRDS(actdf,file=fileout)
saveRDS(mergeddf,file=fileout)
proc.time() - ftm
