library(parallel)
createmerged <- function (df1,df2) {
  
  rsem <- df1
  fpkm <- df2
  
  inrsem <- setdiff(colnames(rsem),colnames(fpkm))
  for (i in inrsem) { rsem[,i] <- NULL}
  
  #remove samples that are in fpkm but not in rsem
  infpkm <- setdiff(colnames(fpkm),colnames(rsem))
  for (i in infpkm) { fpkm[,i] <- NULL}
  mergedDF <- merge(rsem,fpkm,by=0)
  rownames(mergedDF) <- mergedDF$Row.names
  mergedDF$Row.names <- NULL
  
  return(mergedDF)
  
  
  
}



analysis <- function (M,plotname,outdir) {

  par(pch=".")
  
  
  mergedDF <- NULL
  mergedDF <- M



 #work with only samples that are in both data set
  mergedDF <-  mergedDF[grep("[xy]$",colnames(mergedDF),perl=TRUE)]

  full <- NULL
  full <- mergedDF
  
  # corelation
  samples <- NULL
  samples <- colnames(mergedDF)[grep("[x|y]$",colnames(mergedDF),perl=TRUE)]
  samples <- gsub(".[x,y]$","",samples,perl=T)
  samples <- unique(samples)
  corf <- c()
  length(samples)

  for (i in samples) { 
	corf <- append(corf,cor(mergedDF[,paste(i,".x",sep="")],mergedDF[,paste(i,".y",sep="")],method="spearman")) 
	}

  fileout <- paste(outdir,"corf.rds",sep="")
  saveRDS(corf,file=fileout)


  classv <- NULL
  classv <- unlist(as.list(full["cl",]))

  #MAD
  fullmad <- NULL
  fullmad <- apply(mergedDF[(rownames(full)!= "cl"),],1,function (x,y) tapply(x,y,mad),classv)
  
  #MEDIAN
  fullmedian <- NULL
  fullmedian <- apply(mergedDF[(rownames(full)!= "cl"),],1,function (x,y) tapply(x,y,median),classv)

  #gene counts
  fg4 <- function (x) { ifelse( x>4,1,0) }
  fg0 <- function (x) { ifelse( x>0,1,0) }
  
  genecount0 <- as.data.frame(apply(mergedDF[(rownames(mergedDF) != "cl"),], 2, fg0))
  genecount4 <- as.data.frame(apply(mergedDF[(rownames(mergedDF) != "cl"),], 2, fg4))
  
  rsemgc0 <- colSums(genecount0[,grepl(".x",colnames(genecount0))])
  rsemgc4 <- colSums(genecount4[,grepl(".x",colnames(genecount4))])
  
  fpkmgc0 <- colSums(genecount0[,grepl(".y",colnames(genecount0))])
  fpkmgc4 <- colSums(genecount4[,grepl(".y",colnames(genecount4))])
  
  genecountplot <- data.frame(rsemgc0,fpkmgc0)
  
  myc <- makeCluster(16)

  wilcoxp <- parApply(myc,full[!(row.names(full) == "cl"),],1,function(x,y) wilcox.test(x~y)$p.value,t(full["cl",]))
  fileout <- paste(outdir,"wilcoxp_raw.rds",sep="")
  saveRDS(wilcoxp, file=fileout)

  corrected_wilcoxp <- p.adjust(wilcoxp,method="fdr")
  fileout <- paste(outdir,"wilcoxp_fdr.rds",sep="")
  saveRDS(corrected_wilcoxp, file=fileout)
  
  
  classv <- unlist(as.list(full["cl",]))
  fullmean <- apply(full[(rownames(full)!= "cl"),],1,function (x,y) tapply(x,y,mean),classv)
  #fullmean <- parApply(myc,full[(rownames(full)!= "cl"),],1,function (x,y) tapply(x,y,mean),classv)
  fm <- as.data.frame(fullmean)
  log2fc <- log2(fm["38",]) - log2(fm["19",])
  fm["log2fc",] <- log2fc
  tfm <- t(fm)
  tfmdf <- as.data.frame(tfm)
  vol <- NULL
  vol <- data.frame(tfmdf,corrected_wilcoxp)
  vol <- data.frame(tfmdf[(rownames(tfmdf) != "cl"),],corrected_wilcoxp)
  #plot(vol$log2fc,-log10(vol$corrected_wilcoxp))

  fileout <- paste(outdir,"volcano_plotdata.rds",sep="")
  saveRDS(vol, file=fileout)
  
  
  stopCluster(myc)
  
  #PLOT1 correlation
  library(vioplot)
  vioplot(corf)
  mt <- paste(plotname,"Sample pairwise corelation.")
  mt <- paste(c(mt," N=",length(samples)),collapse="")
  title(ylab="spearman coefficient",main=mt)

  vol$AofMA <- (log2(vol$X38) + log2(vol$X19))/2
  #plot(vol$AofMA,vol$log2fc)
  
#basalluminalradio(mg)
#PLOT2
genecountplot <- data.frame(rsemgc0,fpkmgc0)
fileout <- paste(outdir,"genecount0_plotdata.rds",sep="")
saveRDS(genecountplot,file=fileout)
mt <- paste (plotname, "Gene Counts Per Sample with at least one read")
boxplot(rsemgc0,fpkmgc0,ylim=c(12000,20000),names=c("legacy","current"),outline=FALSE,main=mt)

#mt <- paste (plotname, "Gene Counts Per Sample with at least 5 reads")
#boxplot(rsemgc4,fpkmgc4,ylim=c(12000,20000),names=c("legacy","current"), main=mt)

#mt <- paste (plotname, "Legacy Vs Current Gene Expression Median Values")
#boxplot(log10(fullmedian["19",]),log10(fullmedian["38",]),names=c("legacy","current"),ylab="log10 gene median",main=mt)

#mt <- paste (plotname, "Legacy Vs Current Gene Expression mean Values")
#boxplot(log10(fullmean["19",]),log10(fullmean["38",]),names=c("legacy","current"),ylab="log10 gene mean",main=mt)


#mt <- paste (plotname, "Gene Mean Absolute Deviation(MAD)")
#smoothScatter(log2(fullmad["19",]),log2(fullmad["38",]),main=mt ,xlab="log2. MAD. Legacy",ylab="log2. MAD. Current",xlim=c(0,25),ylim=c(0,25))
#plotfitlines(log2(fullmad["19",]),log2(fullmad["38",]))


genemedianplot <- data.frame(x=log2(fullmedian["19",]),y=log2(fullmedian["38",]))
fileout <- paste(outdir,"genemedian_plotdata.rds",sep="")
saveRDS(genemedianplot, file=fileout)

mt <- paste (plotname, "Gene Median")
#smoothScatter(log2(fullmedian["19",]),log2(fullmedian["38",]),main=mt,xlab="log2. Median. Legacy",ylab="log2. Median. Current",xlim=c(0,25),ylim=c(0,25))
smoothScatter(log2(fullmedian["19",]),log2(fullmedian["38",]),main=mt,xlab="log2. Median. Legacy",ylab="log2. Median. Current")
plotfitlines(log2(fullmedian["19",]),log2(fullmedian["38",]))

#PLOT3 MAPLOT
mt <- paste (plotname, "MA plot")
plot(vol$AofMA,vol$log2fc, main=mt, xlab ="log2. Gene Average across the two datasets", ylab="log2. fold change")
abline(h=0)

mt <- paste (plotname, "Volcano plot")
plot(vol$log2fc,-log10(vol$corrected_wilcoxp), main=mt, xlab="log2. fold change",ylab="-log10. p-values",pch=".")
  
}
  
  
  
  

basalluminalratio <- function(classdf,datadf,plotname) {

basal <- classdf
mergedDF <- datadf



basal$ncall[basal$Call == "Basal"] <- 1
basal$ncall[basal$Call != "Basal"] <- 0
b <- data.frame(sample=basal$Barcode,call=basal$ncall)
rownames(b) <- basal$Barcode
b$sample<- NULL
rownames(b) <- gsub("-","\\.",rownames(b))

tb <- t(b)
mtb <- merge(tb,tb,by=0)
rownames(mtb) <- c("call")


mergedDF2 <- rbind(mergedDF,mtb[,names(mergedDF)],make.row.names="T")

ct <- unlist(as.list(mergedDF2["call",]))

legacy <- mergedDF2[ grep("x$",colnames(mergedDF2),perl=TRUE) ]
current <- mergedDF2[ grep("y$",colnames(mergedDF2),perl=TRUE) ]

ct <- unlist(as.list(legacy["call",]))
legacybasal <- apply(legacy,1, function(x,y) tapply(x,y,mean,na.rm=TRUE),ct)

ct <- unlist(as.list(current["call",]))
currentbasal <- apply(current,1, function(x,y) tapply(x,y,mean,na.rm=TRUE),ct)


tlb <- t(legacybasal)
colnames(tlb) <- c("luminal","basal")
xlb <- log2(tlb[,"basal"]) - log2(tlb[,"luminal"]) 

tcb <- t(currentbasal)
colnames(tcb) <- c("luminal","basal")
ylb <- log2(tcb[,"basal"]) - log2(tcb[,"luminal"])

#PLOT4 
mt <- paste(plotname,"Relative Change. BRCA Basal vs non-Basal")
plot(xlb,ylb,main=mt ,xlab="Legacy. Mean. Log2(basal/non-basal)",ylab="Current. Mean. Log2(basal/non-basal)",pch=".",col="green")
plotfitlines(xlb,ylb)
abline(1,1,lty=2)
abline(-1,1,lty=2)

basalplot <- data.frame(xlb,ylb)

fileout <- paste(outdir,"ratioplot.rds",sep="")
saveRDS(basalplot, file=fileout)

k <- data.frame(xlb,ylb)
k1 <- k[((k$xlb - k$ylb > 1) | (k$xlb - k$ylb < -1)),]
k1 <- k1[!is.na(k1$ylb),]
k1 <- k1[!is.na(k1$xlb),]
k1 <- k1[is.finite(k1$ylb),]
k1 <- k1[is.finite(k1$xlb),]
points(k1,pch=".",col="blue")
return(k1)


}

plotfitlines <- function (x,y) {
	plotfit <- data.frame(x=x,y=y)
	plotfitfinite <- plotfit[is.finite(rowSums(plotfit)),]
	lsmad <- lm(plotfitfinite$y~plotfitfinite$x)
	mtext <- paste ("adjusted r square", substr(summary(lsmad)$adj.r.squared,1,5))
	abline(lsmad,col="blue")
	abline(0,1)
	mtext(mtext,side=3)
        return(summary(lsmad)$adj.r.squared)
	}

scalematrix <- function (df,classrowname) {
  

  mergedDF <- df
  c <- classrowname
 
  
  full <- mergedDF[(row.names(mergedDF) != c),]
  full <- scale(full)
  full <- rbind(full,cl=mergedDF[c,])

  return(full)
}


add_subtype_row <- function(adf,filename,var_subtype,var_barcode,pat) {

st <- NULL
tst <- NULL


st <- read.table(file=filename,sep="\t",header=TRUE)
st <- data.frame(barcode=st[,var_barcode],subtype=st[,var_subtype])
st$barcode <- as.character(st$barcode)
#st$subtype <- as.numeric(st$subtype)
rownames(st) <- st$barcode
st$barcode <- NULL

tst <- t(st)
tst <- as.data.frame(tst,stringsAsFactors = FALSE)

for (i in names(adf)) { v <- gsub(pat,"",i)
        if ( pmatch(v, colnames(tst),nomatch=0) > 0)
        {
                adf[var_subtype,i] <- tst["subtype",v]
                }
        else { adf[,i] <- NULL  }


}

return(adf)
}



ratio_analysis <- function (m,dataset_var,subtype_var,plotname) {


        legacy <- NULL
        legacy.ge <- NULL
        legacy.class <- NULL
        legacy.ge.mean <- NULL
        legacy.ge.median <- NULL
        ratio_legacy <- NULL
        legacy_class_ratio <- NULL

        active <- NULL
        active.ge <- NULL
        active.class <-NULL
        active.ge.mean <- NULL
        active.ge.median <- NULL
        ratio_active <- NULL
        active_class_ratio <- NULL

        legacy <- m[grep("\\.x$",colnames(m),perl=TRUE)]
        legacy.ge <- legacy[which(!rownames(m) %in% c(dataset_var,subtype_var)),]
        legacy.class <- legacy[which(rownames(m) %in% subtype_var ),]
        legacy.ge.mean <- apply(legacy.ge,1,function(x,y) tapply(x,y,mean),unlist(legacy.class))
        legacy.ge.median <- apply(legacy.ge,1,function(x,y) tapply(x,y,median),unlist(legacy.class))


        legacy.ge.mean <- as.data.frame(legacy.ge.mean)
        legacy.ge.mean <- t(legacy.ge.mean)
        legacy.ge.mean <- as.data.frame(legacy.ge.mean)
        legacy.ge.mean$log2ratio  <- log2 (legacy.ge.mean$`1` +1) - log2 (legacy.ge.mean$`2` + 1)




        active <- m[grep("\\.y$",colnames(m),perl=TRUE)]
        active.ge <- active[which(!rownames(m) %in% c(dataset_var,subtype_var)),]
        active.class <- active[which(rownames(m) %in% subtype_var),]
        active.ge.mean <- apply(active.ge,1,function(x,y) tapply(x,y,mean),unlist(active.class))
        active.ge.median <- apply(active.ge,1,function(x,y) tapply(x,y,median),unlist(active.class))


        active.ge.mean <- as.data.frame(active.ge.mean)
        active.ge.mean <- t(active.ge.mean)
        active.ge.mean <- as.data.frame(active.ge.mean)
        active.ge.mean$log2ratio  <- log2 (active.ge.mean$`1` +1) - log2 (active.ge.mean$`2` + 1)


        #rp is data for ratioplot
        rp <- data.frame(l=legacy.ge.mean$log2ratio,a=active.ge.mean$log2ratio,k=rownames(legacy.ge.mean))
        rp$color <- "blue"
        rp$color[(rp$l - rp$a) > 1 | (rp$l - rp$a) < -1 ] = "red"

        mt <- paste (plotname, "Ratio")




        plot(rp$l,rp$a,col=rp$color,main=mt,ylab="gene count ratio by class. Active",xlab="gene count ratio by class. Legacy")
        abline(1,1,lty=2)
        abline(-1,1,lty=2)
        plotfitlines(rp$l,rp$a)


        #am <- t(active.ge.median)
        #plot(density(log2(am[,"1"]+1)))
        #plot(density(log2(am[,"2"]+1)))

        #am <- t(legacy.ge.median)
        #plot(density(log2(am[,"1"]+1)))
        #plot(density(log2(am[,"2"]+1)))



        }


create_df_active_v2 <- function(dir_name,pattern) {

        files  <- list.files(path=dir_name, pattern=pattern,full.names=FALSE)
        mydf    <- do.call(cbind,lapply(files,function(fn)read.table(paste(dir_name,fn,sep=""),header=TRUE, sep="\t")[,2]))
        genes <- read.table(paste(dir_name,files[1],sep=""), header=TRUE, sep="\t")[,1]
        rownames(mydf) <- genes

        colnames(mydf) <- files
       mydf <- as.data.frame(mydf)
        return(mydf)

}






class_analysis <- function (m,dataset_var,subtype_var,plotname,outdir,cancer_var) {


        legacy <- NULL
        legacy.ge <- NULL
        legacy.class <- NULL
        legacy.ge.mean <- NULL
        legacy.ge.median <- NULL
        ratio_legacy <- NULL
        legacy_class_ratio <- NULL

        active <- NULL
        active.ge <- NULL
        active.class <-NULL
        active.ge.mean <- NULL
        active.ge.median <- NULL
        ratio_active <- NULL
        active_class_ratio <- NULL

        legacy <- m[grep("\\.x$",colnames(m),perl=TRUE)]
        legacy.ge <- legacy[which(!rownames(m) %in% c(dataset_var,subtype_var)),]
        legacy.class <- legacy[which(rownames(m) %in% subtype_var ),]
        legacy.ge.mean <- apply(legacy.ge,1,function(x,y) tapply(x,y,mean),unlist(legacy.class))



        legacy.ge.mean <- as.data.frame(legacy.ge.mean)
        legacy.ge.mean <- t(legacy.ge.mean)
        legacy.ge.mean <- as.data.frame(legacy.ge.mean)

        
        active <- m[grep("\\.y$",colnames(m),perl=TRUE)]
        active.ge <- active[which(!rownames(m) %in% c(dataset_var,subtype_var)),]
        active.class <- active[which(rownames(m) %in% subtype_var),]
        active.ge.mean <- apply(active.ge,1,function(x,y) tapply(x,y,mean),unlist(active.class))


        active.ge.mean <- as.data.frame(active.ge.mean)
        active.ge.mean <- t(active.ge.mean)
        active.ge.mean <- as.data.frame(active.ge.mean)


       # can do upto 4 classes.

       # if dataset has only 2 classes
	if (! is.na(colnames(legacy.ge.mean)[2] )) {

               col1 <- colnames(legacy.ge.mean)[1]
               col2 <- colnames(legacy.ge.mean)[2]

               class12 <- paste(c(col1,"-",col2),sep="",collapse="")
               


        	legacy.ge.mean[,class12]  <- log2 (legacy.ge.mean[,col1] ) - log2 (legacy.ge.mean[,col2] )
        	active.ge.mean[,class12]  <- log2 (active.ge.mean[,col1] ) - log2 (active.ge.mean[,col2] )

		varname <- paste(cancer_var,"-class12abs",sep="")
        	#cmp <- data.frame (abs(legacy.ge.mean[,class12] - active.ge.mean[,class12]),k=rownames(legacy.ge.mean))
        	cmp <- data.frame (abs(legacy.ge.mean[,class12] - active.ge.mean[,class12]))
		colnames(cmp)[1] <- varname
                rownames(cmp) = rownames(legacy.ge.mean)

        	rp <- data.frame(l=legacy.ge.mean[,class12],a=active.ge.mean[,class12])
                mt <- paste(c(plotname,"class",class12),collapse=" ")
		ratio_plot(rp,mt,"Legacy","Active")

  	
		fileout <- paste(outdir,"relative_change_class12.rds",sep="")
  		saveRDS(rp,file=fileout)
		legacy.ge.mean[,class12] <- NULL
				}

       # if dataset has only 3 classes
	if (! is.na(colnames(legacy.ge.mean)[3] )) {

               col1 <- colnames(legacy.ge.mean)[1]
               col2 <- colnames(legacy.ge.mean)[2]
               col3 <- colnames(legacy.ge.mean)[3]

               class13 <- paste(c(col1,"-",col3),sep="",collapse="")
               class23 <- paste(c(col2,"-",col3),sep="",collapse="")


        	legacy.ge.mean[,class13]  <- log2 (legacy.ge.mean[,col1] ) - log2 (legacy.ge.mean[,col3] )
        	legacy.ge.mean[,class23]  <- log2 (legacy.ge.mean[,col2] ) - log2 (legacy.ge.mean[,col3] )


        	active.ge.mean[,class13]  <- log2 (active.ge.mean[,col1] ) - log2 (active.ge.mean[,col3] )
        	active.ge.mean[,class23]  <- log2 (active.ge.mean[,col2] ) - log2 (active.ge.mean[,col3] )

        	rp <- data.frame(l=legacy.ge.mean[,class13],a=active.ge.mean[,class13],k=rownames(legacy.ge.mean))
                mt <- paste(c(plotname,"class",class13),collapse=" ")
		ratio_plot(rp,mt,"Legacy","Active")

        	rp <- data.frame(l=legacy.ge.mean[,class23],a=active.ge.mean[,class23],k=rownames(legacy.ge.mean))
                mt <- paste(c(plotname,"class",class23),collapse=" ")
		ratio_plot(rp,mt,"Legacy","Active")

		varname <- paste(cancer_var,"-class13abs",sep="")
        	cmp <- cbind(cmp,abs(legacy.ge.mean[,class13] - active.ge.mean[,class13]))
		colnames(cmp)[2] <- varname

		varname <- paste(cancer_var,"-class23abs",sep="")
        	cmp <- cbind(cmp,abs(legacy.ge.mean[,class23] - active.ge.mean[,class23]))
		colnames(cmp)[3] <- varname

		legacy.ge.mean[,class13] <- NULL
		legacy.ge.mean[,class23] <- NULL
		}

       # if dataset has only 4 classes
	if (! is.na(colnames(legacy.ge.mean)[4] )) {
               col1 <- colnames(legacy.ge.mean)[1]
               col2 <- colnames(legacy.ge.mean)[2]
               col3 <- colnames(legacy.ge.mean)[3]
               col4 <- colnames(legacy.ge.mean)[4]

               class14 <- paste(c(col1,"-",col4),sep="",collapse="")
               class24 <- paste(c(col2,"-",col4),sep="",collapse="")
               class34 <- paste(c(col3,"-",col4),sep="",collapse="")

        	legacy.ge.mean[,class14]  <- log2 (legacy.ge.mean[,col1] ) - log2 (legacy.ge.mean[,col4] )
        	legacy.ge.mean[,class24]  <- log2 (legacy.ge.mean[,col2] ) - log2 (legacy.ge.mean[,col4] )
        	legacy.ge.mean[,class34]  <- log2 (legacy.ge.mean[,col3] ) - log2 (legacy.ge.mean[,col4] )

        	active.ge.mean[,class14]  <- log2 (active.ge.mean[,col1] ) - log2 (active.ge.mean[,col4] )
        	active.ge.mean[,class24]  <- log2 (active.ge.mean[,col2] ) - log2 (active.ge.mean[,col4] )
        	active.ge.mean[,class34]  <- log2 (active.ge.mean[,col3] ) - log2 (active.ge.mean[,col4] )

        	rp <- data.frame(l=legacy.ge.mean[,class14],a=active.ge.mean[,class14],k=rownames(legacy.ge.mean))
                mt <- paste(c(plotname,"class",class14),collapse=" ")
		ratio_plot(rp,mt,"Legacy","Active")

        	rp <- data.frame(l=legacy.ge.mean[,class24],a=active.ge.mean[,class24],k=rownames(legacy.ge.mean))
                mt <- paste(c(plotname,"class",class24),collapse=" ")
		ratio_plot(rp,mt,"Legacy","Active")

        	rp <- data.frame(l=legacy.ge.mean[,class34],a=active.ge.mean[,class34],k=rownames(legacy.ge.mean))
                mt <- paste(c(plotname,"class",class34),collapse=" ")
		ratio_plot(rp,mt,"Legacy","Active")

		varname <- paste(cancer_var,"-class14abs",sep="")
        	cmp <- cbind(cmp,abs(legacy.ge.mean[,class14] - active.ge.mean[,class14]))
		colnames(cmp)[4] <- varname

		varname <- paste(cancer_var,"-class24abs",sep="")
        	cmp <- cbind(cmp,abs(legacy.ge.mean[,class24] - active.ge.mean[,class24]))
		colnames(cmp)[5] <- varname

		varname <- paste(cancer_var,"-class34abs",sep="")
        	cmp <- cbind(cmp,abs(legacy.ge.mean[,class34] - active.ge.mean[,class34]))
		colnames(cmp)[6] <- varname

		legacy.ge.mean[,class14] <- NULL
		legacy.ge.mean[,class24] <- NULL
		legacy.ge.mean[,class34] <- NULL

				}
  	fileout <- paste(outdir,"DFofRelativeChangeSubtypes.rds",sep="")
  	saveRDS(cmp,file=fileout)


    

        }

ratio_plot <- function (df,mt,xtext,ytext) {

        rp <- NULL
	rp <- df
        rp$color <- "blue"
        rp$color[(rp$l - rp$a) > 1 | (rp$l - rp$a) < -1 ] = "red"
        #mt <- paste (plotname, "Ratio")
        plot(rp$l,rp$a,col=rp$color,main=mt,ylab=ytext,xlab=xtext)
        abline(1,1,lty=2)
        abline(-1,1,lty=2)
        plotfitlines(rp$l,rp$a)


	}

V2_add_subtype_row <- function(adf,filename,var_subtype,var_barcode,pat) {

st <- NULL
tst <- NULL


st <- read.table(file=filename,sep="\t",header=TRUE)
st <- data.frame(barcode=st[,var_barcode],subtype=st[,var_subtype])
st$barcode <- as.character(st$barcode)
#st$subtype <- as.numeric(st$subtype)
rownames(st) <- st$barcode
st$barcode <- NULL

tst <- t(st)
tst <- as.data.frame(tst,stringsAsFactors = FALSE)

for (i in names(adf)) { v <- gsub(pat,"",i)
        if ( pmatch(v, colnames(tst),nomatch=0) > 0)
        {
                adf[var_subtype,i] <- tst["subtype",v]
                }
        else { adf[,i] <- NULL  }


}

return(adf)
}


maplot <- function(adf,group_var)  {

#adf is a dataframe
#group_var - name of the row that contains the grouping information in adf

       adf.ge <- NULL
       adf.class <- NULL
       adf.ge.mean <- NULL


	 par(pch=".")

	adf.ge <- adf[ !rownames(adf) %in% c(group_var),]
	adf.class <- adf[ rownames(adf) %in% c(group_var),]


	adf.ge.mean <- apply(adf.ge,1, function(x,y) tapply(x,y,mean),unlist(adf.class))

	adf.ge.mean <- t(adf.ge.mean)
	adf.ge.mean <- as.data.frame(adf.ge.mean)

	col1 <- colnames(adf.ge.mean)[1]
	col2 <- colnames(adf.ge.mean)[2]

 

        

	adf.ge.mean[,"M"] <- log2(adf.ge.mean[,col2] ) - log2(adf.ge.mean[,col1] )
	adf.ge.mean[,"A"] <- ( log2(adf.ge.mean[,col1] ) + log2(adf.ge.mean[,col2] ))/2


	plot(adf.ge.mean$A, adf.ge.mean$M, xlab="Gene Average across both datasets", ylab="log ratio across data sets")


       	xt <- paste("Gene Average for Dataset ",col1)
	plot(log2(adf.ge.mean[,col1]), adf.ge.mean$M, xlab=xt, ylab="log ratio across data sets")

       	xt <- paste("Gene Average for Dataset ",col2)
	plot(log2(adf.ge.mean[,col2]), adf.ge.mean$M, xlab=xt, ylab="log ratio across data sets")

	return(adf.ge.mean)


	}


median_min_max <- function(adf,outdir,cancer_var) {

	full <- NULL
	full <- adf

	medvar <- paste(cancer_var,"-median-",sep="")
	maxvar <- paste(cancer_var,"-max-")
	minvar <- paste(cancer_var,"-min-")
	meanvar <- paste(cancer_var,"-mean-")

	classv <- unlist(as.list(full["cl",]))

	fmedian <- apply(full[(rownames(full)!= "cl"),],1,function (x,y) tapply(x,y,median),classv)
	fmedian <- as.data.frame(fmedian)
	fmedian <- t(fmedian)
	fmedian <- as.data.frame(fmedian)
	colnames(fmedian) <- paste(medvar,colnames(fmedian),sep="")
	fmedian$ensembl <- rownames(fmax)

  	fileout <- paste(c(outdir,cancer_var,"-DFMedianofcohorts.rds"),sep="",collapse="")
  	saveRDS(fmedian,file=fileout)

	
	#fmax <- apply(mergedDF[(rownames(full)!= "cl"),],1,function (x,y) tapply(x,y,max),classv)
	#fmax <- as.data.frame(fmax)
	#fmax <- t(fmax)
	#fmax <- as.data.frame(fmax)
	#colnames(fmax) <- paste(maxvar,colnames(fmax),sep="")
	#fmax$ensembl <- rownames(fmax)

	




	}
