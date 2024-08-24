# post CC-GWAS analysis

export DATA=/scratch/c.c23045409/dissertation/ccgwas_input/data 
export CCGWAS_out=/scratch/c.c23045409/dissertation/ccgwas_output
export annot=/scratch/c.c23045409/dissertation/postGWAS/annot
export SMR=/scratch/c.c23045409/dissertation/postGWAS/SMR
export SMR_REF=/scratch/c.c23045409/MET588/Assessment/SMR
export LDclump=/scratch/c.c23045409/dissertation/postGWAS/LD_clumping
export LD_REF=/scratch/c.c23045409/dissertation/postGWAS/postGWAS_FINEMAP/ref_hrc
export FINEMAP=/scratch/c.c23045409/dissertation/postGWAS/postGWAS_FINEMAP


# annotate the significant SNPs 
    cd $annot
    ## use awk to modify gene list to give a 35kb/10kb window around each gene
    awk '($5 == "+"){print $2,($3-35000),($4+10000),$6}' $annot/NCBI37.3.gene.loc > $annot/GeneLocWindow # for positive strand 
    awk '($5 == "-"){print $2,($3-10000),($4+35000),$6}' $annot/NCBI37.3.gene.loc >> $annot/GeneLocWindow # for negative strand
    sed 's/-3991/0/' $annot/GeneLocWindow > $annot/GeneLocWindowB

    # select significant SNPs to annotate only 
    gunzip $CCGWAS_out/test_final.out.results.gz
    awk -F'\t' 'NR==1 || $12 == 1 {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11}' OFS='\t' $CCGWAS_out/test_final.out.results  > $CCGWAS_out/ccgwas.sig.csv
    gzip $CCGWAS_out/test_final.out.results

    # annotate the GWAS results with genes in PLINK
    alias plink="/scratch/c.c23045409/MET588/MET588_finemapping_AFP/software/plink"
    plink --annotate $CCGWAS_out/ccgwas.sig.csv  ranges=$annot/GeneLocWindowB --out ${annot}/ccgwas

    #annotate with attributes 
    plink --annotate $annot/ccgwas.annot attrib=$annot/snp129.attrib.gz --out $annot/ccgwas.attrib # 14 missense annotations 

    #annotate with eQTLs
    plink --annotate $annot/ccgwas.attrib.annot attrib=$annot/Brain_Cortex.eQTL.V7.SNP.Gene.txt --out $annot/ccgwas.eQTL # zero annotations 

# SMR
    #need file with cols: SNP, A1, A2, FRQ_A_67390, b, SE, P, Neff

    awk -F'\t' '{print $1, $2, $3, $4, $5, $6, $7, $8}' OFS='\t' $CCGWAS_out/ccgwas.sig.csv > ${SMR}/ccgwas.SMR.sig.txt

    cd $SMR
    module load R/4.4.0
    R

    library(tidyverse)
    library(data.table)

    BPii <- fread(file = paste0((Sys.getenv("DATA")),"/BPii.txt"), header = T)
    BPi <- fread(file = paste0((Sys.getenv("DATA")),"/BPi.txt"), header = T)
    SMR <- fread(file = paste0((Sys.getenv("SMR")),"/ccgwas.SMR.sig.txt"), header = T)

    BPIb <- select(BPi, "SNP", "FRQ_A_25060", "Neff")
    BPIIb <- select(BPii, "SNP", "FRQ_A_6781", "Neff")

    SMR <- inner_join(SMR, BPIb, by="SNP")
    SMR <- inner_join(SMR, BPIIb, by="SNP")

    SMR_frq_neff <- SMR %>%
        mutate(FRQ_avg = (FRQ_A_25060+FRQ_A_6781)/2) %>%
        mutate(Neff = (Neff.x+Neff.y)/2) %>%
        mutate(FRQ = 1-FRQ_avg)

    SMR_cols <- SMR_frq_neff %>%
        select("SNP", "EA", "NEA", "Neff", "FRQ", "OLS_pval", "OLS_se", "OLS_beta") %>%
        rename(SE = OLS_se) %>%
        rename(P = OLS_pval) %>%
        rename(b = OLS_beta) %>%
        rename(A2 = EA) %>%
        rename(A1 = NEA)

    write.table(SMR_cols, "SMR_input.txt", col.names = T, row.names = F, quote = F)
    quit()


    # perform SMR
    ./smr_Linux --bfile $SMR_REF/g1000_eur --gwas-summary $SMR/SMR_input.txt --beqtl-summary $SMR_REF/westra_eqtl_hg19 --out $SMR/SMR --thread-num 10

    # 1 SNP (rs4719336) was annotated 
    # probeID: ILMN_2358069 
    # ProbeChr: 7       
    # Gene: MAD1L1 
    # Probe_bp: 1937909          
    # topSNP: rs4719336
    # topSNP_chr: 7 
    # topSNP_bp: 1916397       
    # A1: c     
    # A2: T     
    # Freq: 0.423459    
    # b_GWAS: 0.00806 
    # se_GWAS: 0.00127 
    # p_GWAS: 2.080000e-10  
    # b_eQTL: -0.249295  
    # se_eQTL: 0.0192953 
    # p_eQTL: 3.462299e-38     
    # b_SMR: -0.0323312  
    # se_SMR: 0.00567578  
    # p_SMR: 1.224040e-08   
    # p_HEIDI NA
    # nsnp_HEIDI NA

## LD clumping using PLINK
    # create SNP list 
    gunzip $CCGWAS_out/test_final.out.results.gz
    awk 'NR>1 {print $1}' $CCGWAS_out/test_final.out.results > $LDclump/snp_list.txt
    gzip $CCGWAS_out/test_final.out.results.gz

        #extract SNP list from each chromosome 

        for CHR in {1..22}; do
            plink --bfile $LD_REF/HRC.r1-1.EUR.ref_chr${CHR} \
                --extract $LDclump/snp_list.txt \
                --make-bed \
                --out $LDclump/filtered_chr${CHR}
        done

            #create a merge_list file
        cd $LDclump
        for CHR in {1..22}; do
            echo "filtered_chr${CHR}" >> merge_list.txt
        done

        # make bed for all SNPs in the merged list 
        plink --merge-list merge_list.txt \
        --make-bed \
        --out $LDclump/filtered_all

            rm filtered_chr*

    # Clumping 
         # only need SNP and p-value 
        gunzip $CCGWAS_out/test_final.out.results.gz
        awk 'NR>0 {print $1, $8}' $CCGWAS_out/test_final.out.results > $LDclump/ccgwas_clump
        sed '1s/OLS_pval/P/' $LDclump/ccgwas_clump > $LDclump/ccgwas_B
        gzip $CCGWAS_out/test_final.out.results

        plink --bfile $LDclump/filtered_all \
        --clump $LDclump/ccgwas_B \
        --clump-r2 0.1 --clump-p1 0.001 --clump-p2 0.05 --clump-kb 3000 \
        --out $LDclump/CCGWAS_clump_0.1

        # parse clump file to obtain chromosome ranges for each clump 
        cd $LDclump
            module load R/4.4.0
            R 
            library(dplyr)
            library(data.table)

            output <- read.table(file = paste0((Sys.getenv("LDclump")),"/CCGWAS_clump_0.1.clumped"), header=TRUE)

            # convert list on output$SP2 into a separate list 
            # top 5 loci 
            for (i in 1:5) {
                clump <- output$SP2[i]
                clump_cleaned <- gsub("\\(1\\)", "", clump)
                clump_listed <- unlist(strsplit(clump_cleaned, ","))
        
                # Convert clump_listed to a data frame
                clump_df <- data.frame(clump_listed)
        
                # Construct the file name based on the iteration index
                filename <- paste0("clump2_", i, ".csv")
        
                # Write the clump_df to a CSV file
                write.csv(clump_df, file = filename, row.names = FALSE, col.names = FALSE)}

            quit()

# FINE MAPPING 
 alias finemap="/scratch/c.c23045409/MET588/MET588_finemapping_AFP/software/finemap_v1.4.2_x86_64"

# preare LD clumps for FINEMAP
# Finemap needs .z file with these headers: "rsid","chromosome","position","allele1","allele2","maf","beta","se"
    cd $FINEMAP

    module load R/4.4.0
    R 
    library(dplyr)
    library(data.table) 

    BPii <- fread(file = paste0((Sys.getenv("DATA")),"/BPii.txt"), header = T)
    BPi <- fread(file = paste0((Sys.getenv("DATA")),"/BPi.txt"), header = T)
    CCGWAS <- fread(file = paste0((Sys.getenv("CCGWAS_out")),"/test_final.out.results.gz"), header = T)
            # nrow(CCGWAS) 5563208 CHECK THIS

    # Finemap needs .z file with these headers: 
        ## "rsid","chromosome","position","allele1","allele2","maf","beta","se"

    # selecting for rows 
    CCGWAS_restricted <- select(CCGWAS, 1:8) # need CHR, SNP, BP, EA, NEA, OLS_beta, OLS_se, OLS_pval
    BPIb <- select(BPi, "SNP", "FRQ_A_25060")
    BPIIb <- select(BPii, "SNP","Nca","FRQ_A_6781")

    # join datasets 
    CCGWAS_fm <- inner_join(CCGWAS_restricted, BPIb, by="SNP")
    CCGWAS_fm <- inner_join(CCGWAS_fm, BPIIb, by="SNP")

    # only extract SNPs from clump 
        #clump1
    clump1 <- read.csv(file = paste0((Sys.getenv("LDclump")),"/clump2_1.csv"), header = T)
    clump1 <- rename(clump1, SNP = clump_listed)
    CCGWAS_clump <- inner_join(CCGWAS_fm, clump1, by="SNP")

    # create FRQ column that is an average of both 
        CCGWAS_z <- CCGWAS_clump %>%
        mutate(FRQ = (FRQ_A_25060+FRQ_A_6781)/2)

    # calculate minor allele frequency # 
        CCGWAS_z$FRQ[CCGWAS_z$FRQ > 0.5] <- 1 -
        CCGWAS_z$FRQ[CCGWAS_z$FRQ > 0.5]

    # select for and rename rows 
    CCGWAS_z <- select(CCGWAS_z, c("SNP","CHR","BP","EA","NEA","FRQ","OLS_beta","OLS_se"))
    colnames(CCGWAS_z) <- c("rsid","chromosome","position","allele1","allele2","maf","beta","se")

    write.table(CCGWAS_z,"CCGWAS1.z",col.names = T, row.names = F, quote = F)

    #clump2
    clump2 <- read.csv(file = paste0((Sys.getenv("LDclump")),"/clump2_2.csv"), header = T)
    clump2 <- rename(clump2, SNP = clump_listed)
    CCGWAS_clump <- inner_join(CCGWAS_fm, clump2, by="SNP")

    # create FRQ column that is an average of both 
        CCGWAS_z <- CCGWAS_clump %>%
        mutate(FRQ = (FRQ_A_25060+FRQ_A_6781)/2)

    # calculate minor allele frequency # 
        CCGWAS_z$FRQ[CCGWAS_z$FRQ > 0.5] <- 1 -
        CCGWAS_z$FRQ[CCGWAS_z$FRQ > 0.5]

    # select for and rename rows 
    CCGWAS_z <- select(CCGWAS_z, c("SNP","CHR","BP","EA","NEA","FRQ","OLS_beta","OLS_se"))
    colnames(CCGWAS_z) <- c("rsid","chromosome","position","allele1","allele2","maf","beta","se")

    write.table(CCGWAS_z,"CCGWAS2.z",col.names = T, row.names = F, quote = F)

    #clump3
    clump3 <- read.csv(file = paste0((Sys.getenv("LDclump")),"/clump2_3.csv"), header = T)
    clump3 <- rename(clump3, SNP = clump_listed)
    CCGWAS_clump <- inner_join(CCGWAS_fm, clump3, by="SNP")

    # create FRQ column that is an average of both 
        CCGWAS_z <- CCGWAS_clump %>%
        mutate(FRQ = (FRQ_A_25060+FRQ_A_6781)/2)

    # calculate minor allele frequency # 
        CCGWAS_z$FRQ[CCGWAS_z$FRQ > 0.5] <- 1 -
        CCGWAS_z$FRQ[CCGWAS_z$FRQ > 0.5]

    # select for and rename rows 
    CCGWAS_z <- select(CCGWAS_z, c("SNP","CHR","BP","EA","NEA","FRQ","OLS_beta","OLS_se"))
    colnames(CCGWAS_z) <- c("rsid","chromosome","position","allele1","allele2","maf","beta","se")

    write.table(CCGWAS_z,"CCGWAS3.z",col.names = T, row.names = F, quote = F)

    #clump4
    clump4 <- read.csv(file = paste0((Sys.getenv("LDclump")),"/clump2_4.csv"), header = T)
    clump4 <- rename(clump4, SNP = clump_listed)
    CCGWAS_clump <- inner_join(CCGWAS_fm, clump4, by="SNP")

    # create FRQ column that is an average of both 
        CCGWAS_z <- CCGWAS_clump %>%
        mutate(FRQ = (FRQ_A_25060+FRQ_A_6781)/2)

    # calculate minor allele frequency # 
        CCGWAS_z$FRQ[CCGWAS_z$FRQ > 0.5] <- 1 -
        CCGWAS_z$FRQ[CCGWAS_z$FRQ > 0.5]

    # select for and rename rows 
    CCGWAS_z <- select(CCGWAS_z, c("SNP","CHR","BP","EA","NEA","FRQ","OLS_beta","OLS_se"))
    colnames(CCGWAS_z) <- c("rsid","chromosome","position","allele1","allele2","maf","beta","se")

    write.table(CCGWAS_z,"CCGWAS4.z",col.names = T, row.names = F, quote = F)

    #clump5
    clump5 <- read.csv(file = paste0((Sys.getenv("LDclump")),"/clump2_5.csv"), header = T)
    clump5 <- rename(clump5, SNP = clump_listed)
    CCGWAS_clump <- inner_join(CCGWAS_fm, clump5, by="SNP")

    # create FRQ column that is an average of both 
        CCGWAS_z <- CCGWAS_clump %>%
        mutate(FRQ = (FRQ_A_25060+FRQ_A_6781)/2)

    # calculate minor allele frequency # 
        CCGWAS_z$FRQ[CCGWAS_z$FRQ > 0.5] <- 1 -
        CCGWAS_z$FRQ[CCGWAS_z$FRQ > 0.5]

    # select for and rename rows 
    CCGWAS_z <- select(CCGWAS_z, c("SNP","CHR","BP","EA","NEA","FRQ","OLS_beta","OLS_se"))
    colnames(CCGWAS_z) <- c("rsid","chromosome","position","allele1","allele2","maf","beta","se")

    write.table(CCGWAS_z,"CCGWAS5.z",col.names = T, row.names = F, quote = F)
    quit()

    ##perform finemapping
    #1. calculate LD of clumps using reference data
    #2. prepare FINEMAP infile 
    #3. run FINEMAP
    # 792643 = sample size

    #clump1
    plink --bfile ${LDclump}/filtered_all --a1-allele $FINEMAP/CCGWAS1.z 4 5 --extract $FINEMAP/CCGWAS1.z --r square spaces --out $FINEMAP/CCGWAS1 
    echo "z;ld;snp;config;cred;log;n_samples" > $FINEMAP/CCGWAS1.infile
    echo "CCGWAS1.z;CCGWAS1.ld;CCGWAS1.snp;CCGWAS1.config;CCGWAS1.cred;CCGWAS1.log;792643" >> $FINEMAP/CCGWAS1.infile
    finemap --log --sss --in-files $FINEMAP/CCGWAS1.infile --n-causal-snps 1
        # credible set of 31 in CCGWAS1.cred1

    #2
    plink --bfile ${LDclump}/filtered_all --a1-allele $FINEMAP/CCGWAS2.z 4 5 --extract $FINEMAP/CCGWAS2.z --r square spaces --out $FINEMAP/CCGWAS2
    echo "z;ld;snp;config;cred;log;n_samples" > $FINEMAP/CCGWAS2.infile
    echo "CCGWAS2.z;CCGWAS2.ld;CCGWAS2.snp;CCGWAS2.config;CCGWAS2.cred;CCGWAS2.log;792643" >> $FINEMAP/CCGWAS2.infile
    finemap --log --sss --in-files $FINEMAP/CCGWAS2.infile --n-causal-snps 1
        #credible set of 84 in CCGWAS2.cred1

    #3
    plink --bfile ${LDclump}/filtered_all --a1-allele $FINEMAP/CCGWAS3.z 4 5 --extract $FINEMAP/CCGWAS3.z --r square spaces --out $FINEMAP/CCGWAS3
    echo "z;ld;snp;config;cred;log;n_samples" > $FINEMAP/CCGWAS3.infile
    echo "CCGWAS3.z;CCGWAS3.ld;CCGWAS3.snp;CCGWAS3.config;CCGWAS3.cred;CCGWAS3.log;792643" >> $FINEMAP/CCGWAS3.infile
    finemap --log --sss --in-files $FINEMAP/CCGWAS3.infile --n-causal-snps 1
        #credible set of 11 in CCGWAS3.cred1

    #4
    plink --bfile ${LDclump}/filtered_all --a1-allele $FINEMAP/CCGWAS4.z 4 5 --extract $FINEMAP/CCGWAS4.z --r square spaces --out $FINEMAP/CCGWAS4
    echo "z;ld;snp;config;cred;log;n_samples" > $FINEMAP/CCGWAS4.infile
    echo "CCGWAS4.z;CCGWAS4.ld;CCGWAS4.snp;CCGWAS4.config;CCGWAS4.cred;CCGWAS4.log;792643" >> $FINEMAP/CCGWAS4.infile
    finemap --log --sss --in-files $FINEMAP/CCGWAS4.infile --n-causal-snps 1
        #credible set of 19 in CCGWAS4.cred1
    
    #5
    plink --bfile ${LDclump}/filtered_all --a1-allele $FINEMAP/CCGWAS5.z 4 5 --extract $FINEMAP/CCGWAS5.z --r square spaces --out $FINEMAP/CCGWAS5
    echo "z;ld;snp;config;cred;log;n_samples" > $FINEMAP/CCGWAS5.infile
    echo "CCGWAS5.z;CCGWAS5.ld;CCGWAS5.snp;CCGWAS5.config;CCGWAS5.cred;CCGWAS5.log;792643" >> $FINEMAP/CCGWAS5.infile
    finemap --log --sss --in-files $FINEMAP/CCGWAS5.infile --n-causal-snps 1
        #credible set of 19 in CCGWAS5.cred1

