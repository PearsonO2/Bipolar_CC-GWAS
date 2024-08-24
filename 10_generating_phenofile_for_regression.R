# create file for PRS regression 

###### adding sex info to BDRN.fam file #######
library(tidyverse)
library(tidyverse)
library(data.table)

setwd("~/Documents/Masters/Dissertation/week6/PRS")

sex <- read.csv("BDRN_sex.csv")
fam <- read.table("M_BDRN.fam")

sex$Sex <- ifelse(sex$Sex == "Male", 1, ifelse(sex$Sex == "Female", 2, NA))
names(fam)[1] <- "IID"

merged <- inner_join(fam, sex, "IID")

fam_sex <- merged %>%
  select(IID, V2, V3, V4, Sex, V6)

#remove duplicates
fam_unique <- fam_sex %>% 
  distinct(IID, .keep_all = TRUE) # no duplicates

write.table(fam_sex, "M_BDRN_sex.fam", col.names = F, row.names = F, quote = F)


qc <- fread("M_BDRN.fam", col.names)

######### adding pheno cols to fam file ##########

library(tidyverse)
library(data.table)

setwd("~/Documents/Masters/Dissertation/week6/pre-PRS")

pheno <- read.csv("BDRN_pheno.csv")
fam <- read.table("M_BDRN.fam")
sex <- read.csv("BDRN_sex.csv")
eigenvec <- fread("../../week8/M_BDRN.PCA.eigenvec")

names(fam)[1] <- "IID"
names(pheno)[2] <- "IID2"
names(pheno)[1] <- "Pheno_ID"
sex$Sex <- ifelse(sex$Sex == "Male", 1, ifelse(sex$Sex == "Female", 2, NA))

mergeA <- inner_join(fam, sex, "IID")

##### add PCs #####
mergePC <- inner_join(mergeA, eigenvec, "IID")

##### add phenotype variables #####
mergePC$IID2 <- sub("_.*", "", mergePC$IID)

mergePC <- mergePC %>%
  select(IID, ARRAY, PHENO, Pheno_ID, Sex, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, IID2)

pheno <- pheno %>%
  select(1:12)

mergeB <- inner_join(mergePC, pheno, "IID2")
mismatched_phenoID <- mergeB[mergeB$Pheno_ID.x != mergeB$Pheno_ID.y, ]
ID_to_remove <- mismatched_phenoID$IID
filteredB <- mergeB[!mergeB$IID %in% ID_to_remove, ]

#select cols
selB <- filteredB %>%
  select(-Pheno_ID.y) %>%
  rename(Pheno_ID = Pheno_ID.x)

write.table(selB, "phenotypes_ID.txt", col.names=T)

phenotypes <- read.table("phenotypes_ID.txt")
write.csv(phenotypes, "test.txt", col.names = T, row.names = F, quote = F)












