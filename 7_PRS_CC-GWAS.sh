# PRS CC-GWAS

export PRS_DATA=/scratch/c.c23045409/dissertation/postGWAS/PRS/DATA
export PRS_dir=/scratch/c.c23045409/dissertation/postGWAS/PRS/dir
export BDRN=/scratch/c.c23045409/dissertation/postGWAS/PRS/test
export CCGWAS=/home/c.c23045409/dissertation/CCGWAS-master
export LD=/scratch/c.c23045409/dissertation/ccgwas_input/LDSR/ref

#CCGWAS prep 
awk '$2 !~ /^rs/' daner_pgc3_BDI_noBDRN > BDI_sort.txt 
awk '{$2 = gensub(/*:([0-9]+)_.*/, "rs//1", "g", $2); print;}' daner_pgc3_BDI_noBDRN > BDI_sorted.txt #sorting out SNP names

    #making base data input files for CC-GWAS
    cd $PRS_DATA
    module load R/4.4.0
    R 

    library(data.table)
    library(R.utils)
    library(tidyverse) 

    BPI <- fread(file = paste0((Sys.getenv("PRS_DATA")),"/daner_pgc3_BDI_noBDRN"), fill=TRUE) # load BPII file
    nrow(BPI) # 7436527
    head(BPI)

    BPi_wrang <- BPI %>%
    rename(EA=A1, NEA=A2, FRQ=FRQ_U_155439) %>%
    mutate(Neff = Neff_half * 2) %>%
    select(SNP, CHR, BP, EA, NEA, FRQ, OR, SE, P, Neff) 

    head(BPi_wrang) # check file 

    fwrite(BPi_wrang, "PGC3_BDI_noBDRN.neff.gz", compress = "gzip") # save file

    BPII <- fread(file = paste0((Sys.getenv("PRS_DATA")),"/daner_pgc3_BDII_noBDRN"), fill=TRUE) # load BPII file 
    nrow(BPII) #6751102
    head(BPII)

    BPii_wrang <- BPII %>%
    rename(EA=A1, NEA=A2, FRQ=FRQ_U_73009) %>%
    mutate(Neff = Neff_half * 2) %>%
    select(SNP, CHR, BP, EA, NEA, FRQ, OR, SE, P, Neff)

    head(BPii_wrang) # check file 
    fwrite(BPii_wrang, "PGC3_BDII_noBDRN.neff.gz", compress = "gzip") #save file 

#perform CCGWAS (still in R environment)

    library(devtools)
    #install_github("wouterpeyrot/CCGWAS")
    library(CCGWAS)
    source(file = paste0((Sys.getenv("CCGWAS")), "/R/CCGWAS.R"))

    CCGWAS( outcome_file = "PRS.out" , A_name = "BPI" , B_name = "BPII" , 
            sumstats_fileA1A0 = "PGC3_BDI_noBDRN.neff.gz" ,
            sumstats_fileB1B0 = "PGC3_BDII_noBDRN.neff.gz" ,
            K_A1A0 = 0.006 , K_A1A0_high = 0.01 , K_A1A0_low = 0.003 ,  
            K_B1B0 = 0.004 , K_B1B0_high = 0.01 , K_B1B0_low = 0.002 , 
            h2l_A1A0 = 0.1687 , h2l_B1B0 = 0.0693 , rg_A1A0_B1B0 = 0.9571 , intercept_A1A0_B1B0 = 0.171 , m = 8788,  
            N_A1 = 19578 , N_B1 = 4508 , N_A0 = 155439 , N_B0 = 73009 , N_overlap_A0B0 = 45169)


