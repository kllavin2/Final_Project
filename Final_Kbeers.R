#Get the data set for Reference series GSE3292
#NCBI GDS1667
library(GEOquery)
getGEOSuppFiles("GSE3292")
gsd <- getGEO("GSE3292", GSEMatrix = TRUE)
