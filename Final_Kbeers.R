#Get the data set for Reference series GSE3292
#NCBI GDS1667
#annotation file GPL570
library(Biobase)
library(GEOquery)
gds <- getGEO(filename = "GDS1667_full.soft.gz")
# get annotation file
gpl <- getGEO(filename = "GPL570.annot")
# Transform the object into an expression set object




