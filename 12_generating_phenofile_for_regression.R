# create file for PRS regression 


export PRS_R=/scratch/c.c23045409/dissertation/postGWAS/PRS/Ranalysis
export BDRN=/scratch/c.c23045409/dissertation/postGWAS/PRS/test/raw
export QC=/scratch/c.c23045409/dissertation/postGWAS/PRS/QC/targetdata

cd $PRS_R
module load R/4.4.0
R 

######### adding pheno cols to fam file ##########

library(tidyverse)
library(data.table)

pheno <- read.csv(file = paste0((Sys.getenv("BDRN")),"/BDRN_pheno.csv"), header=TRUE)
fam <- read.table(file = paste0((Sys.getenv("QC")),"/M_BDRN.qc.fam"), header=TRUE)
sex <- read.csv(file = paste0((Sys.getenv("BDRN")),"/BDRN_sex.csv"), header=TRUE) # for array information
eigenvec <- fread(file = paste0((Sys.getenv("QC")),"/M_BDRN.PCA.eigenvec"))
PRS <- fread(file = paste0((Sys.getenv("PRS_R")),"/sumPRScs.txt"), fill=TRUE)

colnames(fam) <- c("IID", "FID", "PID", "MID", "sex", "pheno") # assigning col names to fam file
names(pheno)[2] <- "IID2" # renaming cols in pheno file
names(pheno)[1] <- "Pheno_ID"

##### add PCs #####
mergePC <- inner_join(fam, eigenvec, "IID")
mergedArray <- inner_join(mergePC, sex, "IID")

##### add phenotype variables #####
mergedArray$IID2 <- sub("_.*", "", mergedArray$IID) # create new column that is first half of IID

mergeSel <- mergedArray %>%
  select(IID, ARRAY, PHENO, Pheno_ID, sex, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, IID2)

pheno <- pheno %>%
  select(1:12) # selecting required columns from phenotype file

mergeB <- inner_join(mergeSel, pheno, "IID2")
mismatched_phenoID <- mergeB[mergeB$Pheno_ID.x != mergeB$Pheno_ID.y, ] # check for rows that do not have a matching pheontype ID 
ID_to_remove <- mismatched_phenoID$IID # 2 removed
filteredB <- mergeB[!mergeB$IID %in% ID_to_remove, ] # remove these rows 

#select cols
selB <- filteredB %>%
  select(-Pheno_ID.y) %>%
  rename(Pheno_ID = Pheno_ID.x)

# add PRS colum
merged <- inner_join(PRS, selB, by = "IID") #3446 rows 

write.csv(merged, "PRS_complete.txt", col.names = T, row.names = F, quote = F) 















