#Get the data set for Reference series GSE3292
#NCBI GDS1667
#annotation file GPL570
library(Biobase)
library(GEOquery)
gds <- getGEO("GDS1667")
# Convert to expression set
eset <- GDS2eSet(gds, do.log2 = TRUE)
#get sample names
smp <- sampleNames(eset)
#get gene annotation information
gpl <- getGEO(filename="GPL570.annot.gz")
#create a datatable 
MA <- GDS2MA(gds, GPL = gpl)
dat <-data.frame(MA$M)
#create column names
test <- MA[["targets"]][["sample"]]
test2 <-MA[['targets']][["infection"]]
colnames(dat)<- c(paste(test,test2))
#create rownames
rn <- MA$genes$ID
row.names(dat)<- rn
#find outliers
#correlation plot
library(gplots)
dat.cor <- cor(dat)
dat.cor <-as.matrix(dat.cor)
layout(matrix(c(1,1,1,1,1,1,1,1,2,2), 5, 2, byrow = TRUE))
par(oma=c(5,7,1,1))

cx <- rev(colorpanel(25,"blue","white","red"))
leg <- seq(min(dat.cor,na.rm=T),max(dat.cor,na.rm=T),length=10)

image(dat.cor,main="Correlation plot of HPV Pos and HPV Neg samples",axes=F,col=cx)
axis(1,at=seq(0,1,length=ncol(dat.cor)),label=dimnames(dat.cor)[[2]],cex.axis=0.9,las=2)
axis(2,at=seq(0,1,length=ncol(dat.cor)),label=dimnames(dat.cor)[[2]],cex.axis=0.9,las=2)

image(as.matrix(leg),col=cx,axes=F)
tmp <- round(leg,2)
axis(1,at=seq(0,1,length=length(leg)),labels=tmp,cex.axis=1)
#PCA


#Hierarchial Clustering
#transpose the data
dat.t <- t(dat)
#get pairwise distance
dat.dist <- dist(dat.t)
#calculate and plot
par(mar = rep(2, 4))
plot(hclust(dat.dist), 
     labels= dimnames(dat)[[2]], 
     main="Hierarchical Clustering Dendrogram of HPV Pos and NEG Tumors",
     cex.main= 0.75)
