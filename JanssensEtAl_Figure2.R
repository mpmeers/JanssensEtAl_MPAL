#!usr/bin/Rscript

## Collect arguments
args <- commandArgs(TRUE)
 
## Default setting when no arguments passed
if(length(args) < 2) {
  args <- c("--help")
}
 
## Help section
if("--help" %in% args) {
  cat("
			DAlChrE: Generate PCA plots from count data
			
      Arguments:
			--input=someValue   - Input count file from DAlChrE
			--output=someValue   - Output prefix
")
 
  q(save="no")
}

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1
invis <- gc(verbose=FALSE) 

## Arg1 default
#if(is.null(args[1])){
if(is.null(argsL$input) | is.null(argsL$output)) {
  stop("Argument is missing!
			DAlChrE: Generate PCA plots from count data
			
      Arguments:
			--input=someValue   - Input count file from DAlChrE
			--output=someValue   - Output prefix
")
 
  q(save="no")
}

input<-read.table(argsL$input, header=TRUE)

library(MASS)
sums<-apply(input[,c(4:ncol(input))], 1, sum) ## Added 12/19/20 to add a lower bound filter for total signal in peak across all samples
fit<-fitdistr(input$end-input$start, densfun="lognormal")
fit2<-fitdistr(sums, densfun="lognormal")     ## Added 12/19/20 to add a lower bound filter for total signal in peak across all samples
summodel<-rlnorm(nrow(input), fit2$estimate["meanlog"], fit$estimate["sdlog"]) ## Added 12/19/20 to add a lower bound filter for total signal in peak across all samples
input<-input[sums > as.numeric(quantile(summodel, 0.01)),] ## Added 12/19/20 to add a lower bound filter for total signal in peak across all samples
frame<-data.frame(data=input$end-input$start, model=rlnorm(nrow(input), fit$estimate["meanlog"], fit$estimate["sdlog"]))
ecdf<-ecdf(log(frame$data))(seq(min(log(frame$data)),max(log(frame$data)),length=1000))
ecdf.model<-pnorm(seq(min(log(frame$data)),max(log(frame$data)),length=1000), fit$estimate["meanlog"], fit$estimate["sdlog"])
dist2d<-function(a,b,c){v1<- b - c; v2<- a - b; m<-cbind(v1,v2); d<-det(m)/sqrt(sum(v1*v1))}
distances<-data.frame(thresholds=seq(min(log(frame$data)),max(log(frame$data)),length=1000), ecdf=ecdf, model=ecdf.model)
distances$dist<-apply(distances, 1, function(x) dist2d(c(x[2],x[3]),c(0,0),c(1,1)))
if(max(distances$dist) > 0.05){
	input<-input[input$end-input$start > exp(distances$thresholds[distances$dist==max(distances$dist)]),]
}

counts<-input[,c(4:ncol(input))]
#tfidf<-log(apply(counts, 1, function(x) length(x)/length(x[x>0]))) ## New on 6/26/20
#counts<-counts*tfidf ## New on 6/26/20
#counts<-counts[apply(counts,1,sum) > 0,] ## New on 6/26/20
#counts_norm<-counts*(1000/(input$end-input$start))
counts_norm<-counts
#counts_norm_log<-log10(counts_norm+min(counts_norm[counts_norm != 0]))
counts_norm_log<-asinh(counts_norm) ## New on 6/28/20
counts_norm_log_scale<-as.data.frame(scale(counts_norm_log))
outbed<-cbind(input[,c(1:3)],counts_norm_log_scale)
write.table(outbed, file=paste(argsL$output, ".normalized.bed", sep=""), quote=FALSE, sep="\t", row.names=FALSE)
counts_norm_log_scale_pca<-prcomp(counts_norm_log_scale)
comps<-as.data.frame(counts_norm_log_scale_pca$rotation)
vars<-round((((counts_norm_log_scale_pca$sdev)^2)/sum((counts_norm_log_scale_pca$sdev)^2))*100, 2)
var1<-paste("PC1 (",vars[1],"% variance)", sep="")
var2<-paste("PC2 (",vars[2],"% variance)", sep="")
library(ggplot2)
p<-ggplot(comps, aes(x=PC1, y=PC2, label=rownames(comps))) + geom_point(size=2) + geom_text(hjust="inward") + theme_bw() + xlab(var1) + ylab(var2)
ggsave(paste(argsL$output, ".PCA.pdf", sep=""), p)
library(Rtsne)
#counts_norm_log_scale_pca_tsne<-Rtsne(counts_norm_log_scale_pca$x, initial_dims=length(vars[vars > 1]), pca=FALSE, check_duplicates=FALSE)
counts_norm_log_scale_pca_tsne<-Rtsne(counts_norm_log_scale_pca$x, initial_dims=length(vars[vars > 1]), perplexity=round(nrow(counts_norm_log_scale_pca$x)^(1/2), 0), pca=FALSE, check_duplicates=FALSE)
library(densityClust)
denstest<-densityClust(counts_norm_log_scale_pca_tsne$Y)
counts_norm_log_scale_pca_tsne_dclust<-findClusters(denstest, rho=quantile(denstest$rho, 0.95), delta=quantile(denstest$delta, 0.95))
tsneframe<-as.data.frame(counts_norm_log_scale_pca_tsne$Y)
tsneframe$clust<-counts_norm_log_scale_pca_tsne_dclust$clusters
tsneframe$clust<-factor(tsneframe$clust)
tsnetemp<-data.frame(peaks=counts_norm_log_scale_pca_tsne_dclust$peaks, clust=c(1:length(counts_norm_log_scale_pca_tsne_dclust$peaks)))
tsnetemp<-tsnetemp[order(tsnetemp$peaks),]
#tsnelabels<-tsneframe[rownames(tsneframe) %in% tsnetemp$peaks,c(1:2)]
#tsnelabels$clust<-tsnetemp$clust
#tsnelabels$clust<-factor(tsnelabels$clust)
tsnelabels<-as.data.frame(aggregate(tsneframe, by=list(tsneframe$clust), FUN=mean)[,c(1:3)])
colnames(tsnelabels)<-c("clust", "V1", "V2")
q<-ggplot(tsneframe, aes(x=V1, y=V2, color=clust)) + geom_point(size=0.2, pch=16) + theme_light() + xlab("t-SNE 1") + ylab("t-SNE 2") + geom_text(data=tsnelabels, aes(x=V1, y=V2, label=clust), size=5, color="black") + theme(legend.position="none")
ggsave(paste(argsL$output, ".tsne_clust_scatter.pdf", sep=""), q, width=8, height=7)
#tsneout<-rbind(input[,c(1:3)], tsneframe)
clust.info<-aggregate(counts_norm_log_scale, by=list(tsneframe$clust),FUN=mean)
write.table(clust.info, file=paste(argsL$output, ".tsne_clust_meanvalues.bed", sep=""), quote=FALSE, sep="\t", row.names=FALSE)
clust.info.mat<-as.matrix(clust.info[,2:ncol(clust.info)])
len<-nrow(clust.info.mat)/ncol(clust.info.mat)
library(gplots)
library(heatmap3)
library(RColorBrewer)
RowColors<-colorRampPalette(c("blue", "forestgreen", "yellow", "orange", "red", "purple"))(n=length(table(tsneframe$clust)))
my_palette<-colorRampPalette(c("blue","yellow"))(n=99)
col_breaks<-seq(quantile(clust.info.mat, 0.05),quantile(clust.info.mat,0.95),length=100)
pdf(paste(argsL$output, ".tsne_clust_heatmap_mean.pdf", sep=""), height=round(nrow(clust.info)*0.1,0)+4, width=round(ncol(clust.info.mat)*0.1,0)+4)
#out<-heatmap3(clust.info.mat, col=my_palette, RowSideColors=RowColors)
out2<-heatmap.2(clust.info.mat, scale="none", trace="none", col=my_palette, RowSideColors=RowColors)
dev.off()
#out2<-heatmap.2(clust.info.mat, scale="none", trace="none", col=my_palette)
out<-heatmap3(clust.info.mat, col=my_palette, RowSideColors=RowColors)
clustframe<-cbind(tsneframe$clust, counts_norm_log_scale)
scores<-apply(clustframe, 1, function(x) dist(rbind(x[2:length(x)],clust.info[x[1],c(2:ncol(clust.info))]), method="euclidean"))
tsneout<-data.frame(chr=input$chr, start=input$start, end=input$end, tsne1=tsneframe$V1, tsne2=tsneframe$V2, cluster=tsneframe$clust, score=scores)
write.table(tsneout, file=paste(argsL$output, ".tsne_clust_clusters.bed", sep=""), quote=FALSE, sep="\t", row.names=FALSE)
tsneout<-tsneout[order(match(tsneout$cluster, out$rowInd), tsneout$score),]
fullRowColors<-RowColors[tsneout$cluster]
counts_norm_log_scale<-counts_norm_log_scale[order(match(tsneout$cluster, out$rowInd), tsneout$score),]
pdf(paste(argsL$output, ".tsne_clust_heatmap_full.pdf", sep=""), height=round(nrow(clust.info)*0.1,0)+4, width=round(ncol(clust.info.mat)*0.1,0)+4)
heatmap3(as.matrix(counts_norm_log_scale), Rowv=NA, Colv=as.dendrogram(out2$colDendrogram), col=my_palette, RowSideColors=fullRowColors, useRaster=FALSE)
#heatmap.2(as.matrix(counts_norm_log_scale), trace="none", scale="none", Rowv=NA, Colv=as.dendrogram(out2$colDendrogram), col=my_palette, RowSideColors=fullRowColors, useRaster=FALSE)
dev.off()

