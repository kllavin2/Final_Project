#Get the data set for Reference series GSE6631
#NCBI GDS2520
#annotation file GPL8300
library(Biobase)
library(GEOquery)
gds <- getGEO("GDS2520")
# Convert to expression set
eset <- GDS2eSet(gds, do.log2 = TRUE)
#get sample names
smp <- sampleNames(eset)
#get gene annotation information
gpl <- getGEO("GPL8300")
#create a datatable 
MA <- GDS2MA(gds, GPL = gpl)
dat <-data.frame(MA$M)
#create column names
nm <- MA[["targets"]][["sample"]]
ds <-MA[['targets']][["disease.state"]]
ds<-revalue(ds, c("head and neck squamous carcinoma"="carcinoma",
                  "normal"="normal"))
colnames(dat)<- c(paste(ds,nm))
#create rownames
rn <- MA$genes$ID
row.names(dat)<- rn
nrm <- dat[c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43)]
carc <-dat[c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44)]
dat <-cbind(nrm,carc)

#find outliers
#correlation plot
library(gplots)
#calculate correlation
dat.cor <- cor(dat)
dat.cor <-as.matrix(dat.cor)
#set layout
dev.off()
layout(matrix(c(1,1,1,1,1,1,1,1,2,2),
              5, 
              2, 
              byrow = TRUE))
par(oma=c(5,7,1,1))
#correlation plot
cx <- rev(colorpanel(25,
                     "blue","white","red"))
leg <- seq(min(dat.cor,na.rm=T),
           max(dat.cor,na.rm=T),length=10)

image(dat.cor,
      main="Correlation plot of HNSC Normal Vs Carcinoma Samples",
      axes=F,col=cx)
axis(1,at=seq(0,1,length=ncol(dat.cor)),
     label=dimnames(dat.cor)[[2]],
     cex.axis=0.75,las=2)
axis(2,at=seq(0,1,length=ncol(dat.cor)),
     label=dimnames(dat.cor)[[2]],
     cex.axis=0.75,las=2)

image(as.matrix(leg),
      col=cx,axes=F)
tmp <- round(leg,2)
axis(1,at=seq(0,1,length=length(leg)),
     labels=tmp,cex.axis=1)
#hierarchical clustering plot
#transpose the data
dat.t <- t(dat)
#get pairwise distances
dat.dist <-dist(dat.t)
#calculate and plot hierarchial tree
plot(hclust(dat.dist),
     labels = dimnames(dat)[[2]],
     main = "Hierarchial Clustering Dendrogram of HNSC Carcinoma and Normal Samples",
     cex.main= 0.6)
dat.mean <- apply(dat,2,mean)
#calculate sample standard deviations
dat.sd <- apply(dat,2,sd)
#calculate cv of samples
dat.cv <- dat.sd/dat.mean
#create cv plot
plot(dat.mean,
     dat.cv,
     main="HNSC Carcinoma/Normal-Sample CV vs. Mean",
     xlab="Mean",
     ylab="CV",
     col='blue',
     cex=0.75,
     type="n")
text(dat.mean,
     dat.cv,
     label=dimnames(dat)[[2]],
     cex=0.5)

# average correlation plot
dat.avg <- apply(dat.cor,1,mean)
par(oma=c(6,0.1,0.1,0.1))
plot(c(1,length(dat.avg)),
     range(dat.avg),
     type="n",
     xlab="",
     ylab="Avg r",
     main="Avg correlation of HNSC Carcinoma and normal samples",
     axes=F)
points(dat.avg,
       bg="red",
       col=1,
       pch=21,
       cex=1.5)
axis(1,at=c(1:length(dat.avg)),labels=dimnames(dat)[[2]],las=2,cex=0.5)
axis(2)
abline(v=seq(0.5,21.5,1),col="grey")
#Outliers: GSM153813, GSM153834, and GSM153826
#Remove these outliers from the data set
dat.out <- dat[-c(1,33,29)]
#look at distribution of the data
#look at this again http://www.cureffi.org/2013/08/23/gene-expression-analysis-qc-pipeline-in-r/
ddat <-as.numeric(as.matrix((dat.out)))

hist(ddat, breaks = 100, 
     col="green", 
     main = "Histogram of expression levels",
     xlab = "Expression Level")
     
#Find the upper quartile
group <- as.numeric(ddat>quantile(as.numeric(as.matrix(ddat)),0.75))
expgrp <- data.frame(dat.out[,-1],group)
#rowmeans of data
rwavg <- rowMeans(dat)

