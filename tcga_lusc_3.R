rm(list=ls())

outdir <- "~/GDAN//fpkm-uq/TCGA-LUSC/output/"


plot_cor_file <- paste(outdir,"corf.rds",sep="")
plot_genemedian_file <- paste(outdir,"genemedian_plotdata.rds",sep="")
plot_genecount0_file <- paste(outdir,"genecount0_plotdata.rds",sep="")
plot_relativechange12_file <- paste(outdir,"relative_change_class12.rds",sep="")
plot_outputfile <- paste(outdir,"analysis_plots.png",sep="")


corf <- readRDS(file=plot_cor_file)
genecountplot <- readRDS(file=plot_genecount0_file)
genemedianplot <- readRDS(file=plot_genemedian_file)
relativechangeplot <- readRDS(file=plot_relativechange12_file)




cancer_var <- "LUSC"


png(plot_outputfile,width=8, height=8, units = 'in', res=300)

par(mfrow=c(2,2))

#PLOT1
mt <- paste(cancer_var,"Pairwise Sample correlation", sep= " ")
plot(density(corf),xlab="spearman coefficient.",main=mt)


#PLOT2
plotfit <- genemedianplot
plotfitfinite <- plotfit[is.finite(rowSums(genemedianplot)),]
lsmad <- lm(plotfitfinite$y~plotfitfinite$x)
plotfitfinite$resi <- resid(lsmad)
plotfitfinite$color <- "blue"
plotfitfinite$color [plotfitfinite$resi > 2 | plotfitfinite$resi < -2 ] = "red"
mt <- paste(cancer_var,"Log2 Gene Median", sep= " ")
plot(plotfitfinite$x,plotfitfinite$y,col=plotfitfinite$color,pch=".",main=mt ,xlab="Legacy",ylab="Current")
abline(lsmad,col="blue")
abline(0,1)
mtext <- paste ("adj r-square", substr(summary(lsmad)$adj.r.squared,1,5))
legend("bottomright",legend=c(mtext),cex=0.75)

#PLOT3
mt <-  paste(cancer_var,"Number of Genes detected", sep= " ")
boxplot(genecountplot$rsemgc0,genecountplot$fpkmgc0,names=c("Legacy","Current"),outline=FALSE,main=mt,ylab="No. of Genes")

#PLOT4
plotfit <- relativechangeplot
plotfitfinite <- plotfit[is.finite(rowSums(relativechangeplot)),]
lsmad <- lm(plotfitfinite$a~plotfitfinite$l)
plotfitfinite$color <- "blue"
plotfitfinite$color [(plotfitfinite$l - plotfitfinite$a) > 1 | (plotfitfinite$l - plotfitfinite$a) < -1   ] = "red"

mt <- "Relative Change.\nLUSC Basal vs Classical"
plot(plotfitfinite$l,plotfitfinite$a,col=plotfitfinite$color,pch=".",main=mt,xlab="Legacy. Mean. Log2(basal/classical)",ylab="Current. Mean. Log2(basal/classical)")

abline(lsmad,col="blue")
abline(0,1)
abline(1,1,lty=2)
abline(-1,1,lty=2)
mtext <- paste ("adj r-square", substr(summary(lsmad)$adj.r.squared,1,5))
legend("bottomright",legend=c(mtext),cex=0.75)

dev.off()
