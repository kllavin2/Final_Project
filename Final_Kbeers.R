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
#create rownames
rn <- MA$genes$ID
row.names(dat)<- rn
#find outliers
