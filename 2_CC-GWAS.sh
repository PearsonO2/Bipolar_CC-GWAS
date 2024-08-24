#final CC-GWAS

export DATA=/scratch/c.c23045409/dissertation/ccgwas_input/data 
export CCGWAS_out=/scratch/c.c23045409/dissertation/ccgwas_output
export CCGWAS=/home/c.c23045409/dissertation/CCGWAS-master

#prepare data files
    #Column names should be: SNP, CHR, BP, EA, NEA, FRQ, OR, SE, P, Neff. 
    ##EA denotes the effective allele, 
    ##NEA denotes the non effective allele, 
    ##FRQ denotes the frequency of the EA in controls, 
    ##OR denotes the OR per EA, 
    ##SE denotes the standard error of log(OR), 
    ##Neff labels the effective sample size.

    cd $DATA
    module load R/4.4.0
    R 

    library(data.table)
    library(tidyverse)
    BPI <- fread(file = paste0((Sys.getenv("DATA")),"/BPi.txt"), fill=TRUE, header = T) #load data 
    nrow(BPI)
    head(BPI)

    BPi_wrang <- BPI %>%
        select(SNP, CHR, BP, A1, A2, FRQ_U_449978, OR, SE, P, Neff) %>%
        rename(EA=A1, NEA=A2, FRQ=FRQ_U_449978)

    head(BPi_wrang) # check cols 

    fwrite(BPi_wrang, "PGC3_BDI_ccgwas.neff.gz", compress = "gzip") # save file


    BPII <- fread(file = paste0((Sys.getenv("DATA")),"/BPii.txt"), fill=TRUE) # load data 
    nrow(BPII) 
    head(BPII)

    BPii_wrang <- BPII %>%
        select(SNP, CHR, BP, A1, A2, FRQ_U_364075, OR, SE, P, Neff) %>%
        rename(EA=A1, NEA=A2, FRQ=FRQ_U_364075)

    head(BPii_wrang) # check cols 

    fwrite(BPii_wrang, "PGC3_BDII_ccgwas.neff.gz", compress = "gzip") # save file 

#perform CCGWAS 
    library(data.table)
    library(R.utils)
    library(devtools)
    #install_github("wouterpeyrot/CCGWAS")
    library(CCGWAS)
    source(file = paste0((Sys.getenv("CCGWAS")), "/R/CCGWAS.R"))

    CCGWAS( outcome_file = "test_final.out" , A_name = "BPI" , B_name = "BPII" , 
            sumstats_fileA1A0 = "PGC3_BDI_ccgwas.neff.gz" ,
            sumstats_fileB1B0 = "PGC3_BDII_ccgwas.neff.gz" ,
            K_A1A0 = 0.006 , K_A1A0_high = 0.01 , K_A1A0_low = 0.003 ,  
            K_B1B0 = 0.004 , K_B1B0_high = 0.01 , K_B1B0_low = 0.002 , 
            h2l_A1A0 = 0.1842 , h2l_B1B0 = 0.0908 , rg_A1A0_B1B0 = 0.8456 , intercept_A1A0_B1B0 = 0.1817 , m = 8788 ,  
            N_A1 = 25060 , N_B1 = 16781 , N_A0 = 449978 , N_B0 = 364075 , N_overlap_A0B0 = 53251 )


    quit()

    mv $DATA/test_final.* $CCGWAS_out