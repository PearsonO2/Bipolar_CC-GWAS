# BDRN QC
# as per https://choishingwan.github.io/PRS-Tutorial/target/

export PRS_DATA=/scratch/c.c23045409/dissertation/postGWAS/PRS/DATA
export PRS_dir=/scratch/c.c23045409/dissertation/postGWAS/PRS/dir
export BDRN=/scratch/c.c23045409/dissertation/postGWAS/PRS/test
export LD=/scratch/c.c23045409/dissertation/ccgwas_input/LDSR/ref
export PRS=/scratch/c.c23045409/dissertation/postGWAS/PRS
PRScs=/scratch/c.c23045409/dissertation/postGWAS/PRS/PRScs
export PRS_LDSR=/scratch/c.c23045409/dissertation/postGWAS/PRS/LDSR
export QC=/scratch/c.c23045409/dissertation/postGWAS/PRS/QC/targetdata


cd $BDRN

# Need to add sex information to BDRN.fam file 

    module load R/4.4.0
    R

    library(tidyverse)
    library(data.table)

    sex <- read.csv(file = paste0((Sys.getenv("BDRN")),"/BDRN_sex.csv"), header=TRUE)
    fam <- read.table(file = paste0((Sys.getenv("BDRN")),"/M_BDRN.fam"), header=TRUE)

    sex$Sex <- ifelse(sex$Sex == "Male", 1, ifelse(sex$Sex == "Female", 2, NA)) # recoding sex information
    colnames(fam) <- c("IID", "FID", "PID", "MID", "sex", "pheno") # assigning col names to fam file

    merged <- inner_join(fam, sex, "IID")

    fam_sex <- merged %>%
    select(IID, V2, V3, V4, Sex, V6)

    #remove duplicates
    fam_unique <- fam_sex %>% 
    distinct(IID, .keep_all = TRUE) # no duplicates

    write.table(fam_sex, "M_BDRN_sex.fam", col.names = F, row.names = F, quote = F)


    #need to updated bed and bim files 
        #need to update IDs
        cut -d' ' -f2 ${BDRN}/M_BDRN_sex.fam > ${BDRN}/updated_ids.txt # extracting samples with sex info

        module purge
        module load plink/2.0
        plink2 --bfile ${BDRN}/M_BDRN --keep-fam ${BDRN}/updated_ids.txt --make-bed --out ${BDRN}/M_BDRN_filtered # extracting samples with sex info from bed, bim and fam files
        #adding sex info 
        plink2 --bfile ${BDRN}/M_BDRN_filtered --fam ${BDRN}/M_BDRN_sex.fam --make-bed --out ${BDRN}/M_BDRN_updated # creating bed, bim and fam file for 

#QC of target data
    #standard QC'ing
    plink2 --bfile ${BDRN}/M_BDRN_updated --maf 0.1 --geno 0.01 --mind 0.01 --hwe 1e-6 --write-snplist --make-just-fam --out ${QC}/M_BDRN.qc
    
    #pruning 
    plink2 --bfile ${BDRN}/M_BDRN_updated --keep ${QC}/M_BDRN.qc.fam --extract ${QC}/M_BDRN.qc.snplist --indep-pairwise 200 50 0.25 --out ${QC}/M_BDRN.qc

    # interogate heterozygosity 
    module purge 
    module load plink/1.9
    plink --bfile ${BDRN}/M_BDRN_updated  --extract ${QC}/M_BDRN.qc.prune.in --keep ${QC}/M_BDRN.qc.fam --het --out ${QC}/M_BDRN.qc
    
    module purge
    module load R/4.4.0
    cd $QC
    R 
    library(data.table)
    dat <- read.table(file = paste0((Sys.getenv("QC")),"/M_BDRN.qc.het"), header =TRUE)
    #0
        # Calculate mean and standard deviation of the column F
        mean_F <- mean(dat$F, na.rm = TRUE)
        sd_F <- sd(dat$F, na.rm = TRUE)
        # Filter the data frame based on the conditions
        valid <- dat[dat$F <= mean_F + 3 * sd_F & dat$F >= mean_F - 3 * sd_F, ]

    # print FID and IID for valid samples
    fwrite(valid[,c("FID","IID")], "M_BDRN.valid.sample", sep="\t") 
    q() 

    #mismatching SNPs 
    R 
    library(data.table)
    library(magrittr)

    bim <- fread(file = paste0((Sys.getenv("BDRN")),"/M_BDRN_updated.bim")) %>%
    setnames(., colnames(.), c("CHR", "SNP", "CM", "BP", "B.A1", "B.A2")) %>%
    .[,c("B.A1","B.A2"):=list(toupper(B.A1), toupper(B.A2))]  # And immediately change the alleles to upper cases

    # Read in summary statistic data (require data.table v1.12.0+)
        CCGWAS <- fread(file = paste0((Sys.getenv("PRS_DATA")),"/PRS.out.results"), fill=TRUE) %>%
        .[,c("AE","NEA"):=list(toupper(EA), toupper(NEA))] # And immediately change the alleles to upper cases
    
    # Read in QCed SNPs
    qc <- fread(file = paste0((Sys.getenv("QC")),"/M_BDRN.qc.snplist"), fill=FALSE, header=F) 

    # identify SNPs that need flipping 
        # Merge summary statistic with target
        info <- merge(bim, CCGWAS, by=c("SNP", "CHR", "BP")) %>%
        .[SNP %in% qc[,V1]] # And filter out QCed SNPs

    # Function for calculating the complementary allele
    complement <- function(x){
    switch (x,
        "A" = "T",
        "C" = "G",
        "T" = "A",
        "G" = "C",
        return(NA))} 

    # Get SNPs that have the same alleles across base and target
    info.match <- info[B.A1 == B.A1 & B.A2 == B.A2, SNP]
    # Identify SNPs that are complementary between base and target
    com.snps <- info[sapply(B.A1, complement) == B.A1 &
                    sapply(B.A2, complement) == B.A2, SNP]
    # Now update the bim file
    bim[SNP %in% com.snps, c("B.A1", "B.A2") :=
        list(sapply(B.A1, complement),
            sapply(B.A2, complement))]

    #Identify SNPs that require recoding in the target
    recode.snps <- info[B.A1==B.A2 & B.A2==B.A1, SNP]
    # Update the bim file
    bim[SNP %in% recode.snps, c("B.A1", "B.A2") :=
        list(B.A2, B.A1)]

    # identify SNPs that need recoding & complement
    com.recode <- info[sapply(B.A1, complement) == B.A2 &
                    sapply(B.A2, complement) == B.A1, SNP]
    # Now update the bim file
    bim[SNP %in% com.recode, c("B.A1", "B.A2") :=
        list(sapply(B.A2, complement),
            sapply(B.A1, complement))]
    # Write the updated bim file
    fwrite(bim[,c("SNP", "B.A1")], "M_BDRN.a1", col.names=F, sep="\t")

    #Identify SNPs that have different allele in base and target
    mismatch <- bim[!(SNP %in% info.match |
                    SNP %in% com.snps |
                    SNP %in% recode.snps |
                    SNP %in% com.recode), SNP]
    write.table(mismatch, "M_BDRN.mismatch", quote=F, row.names=F, col.names=F)
    q() 

# sex check - data does not include sex chromosome, 7 ambiguius though 
    module load plink/1.9
    plink \
    --bfile ${BDRN}/M_BDRN_updated \
    --extract ${QC}/M_BDRN.qc.prune.in \
    --keep ${QC}/M_BDRN.valid.sample \
    --check-sex \
    --out ${QC}/M_BDRN.qc

    #interrogate in R
    R 
    library(data.table)
    valid <- fread(file = paste0((Sys.getenv("QC")),"/M_BDRN.valid.sample")
    dat <- fread(file = paste0((Sys.getenv("QC")),"/M_BDRN.qc.sexcheck")[FID%in%valid$FID]
    fwrite(dat[STATUS=="OK",c("FID","IID")], "M_BDRN.qc.valid", sep="\t") 
    fwrite("M_BDRN.qc.valid", sep="\t")
    q() 

# relatedness 
    plink    \
    --bfile ${BDRN}/M_BDRN_updated  \
    --extract ${QC}/M_BDRN.qc.prune.in \
    --keep ${QC}/M_BDRN.valid.sample \
    --rel-cutoff 0.125 \
    --out ${QC}/M_BDRN.qc

# generate final file 
plink2 \
    --bfile ${BDRN}/M_BDRN_updated  \
    --make-bed \
    --keep ${QC}/M_BDRN.qc.rel.id \
    --out ${QC}/M_BDRN.qc \
    --extract ${QC}/M_BDRN.qc.snplist \
    --exclude ${QC}/M_BDRN.mismatch 
