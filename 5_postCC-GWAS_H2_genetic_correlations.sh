# post CC-GWAS heritability estimates and genetic correlations 

export DATA=/scratch/c.c23045409/dissertation/ccgwas_input/data 
export CCGWAS_out=/scratch/c.c23045409/dissertation/ccgwas_output
export annot=/scratch/c.c23045409/dissertation/postGWAS/annot
export LDSR_post=/scratch/c.c23045409/dissertation/postGWAS/LDSR
export LDAK_post=/scratch/c.c23045409/dissertation/postGWAS/LDAK
export sex_strat=/scratch/c.c23045409/dissertation/postGWAS/sex_strat
export LD=/scratch/c.c23045409/dissertation/ccgwas_input/LDSR/ref

## heritability estimates
    #LDSR

        cd $LDSR_post
        module load R/4.4.0
        R 
        library(dplyr)
        library(data.table) 

        BPii <- fread(file = paste0((Sys.getenv("DATA")),"/BPii.txt"), header = T)
        BPi <- fread(file = paste0((Sys.getenv("DATA")),"/BPi.txt"), header = T)
        CCGWAS <- fread(file = paste0((Sys.getenv("CCGWAS_out")),"/test_final.out.results.gz"), header = T)

        CCGWAS_ldsr <- CCGWAS %>%
            select("SNP", "CHR", "BP", "EA", "NEA", "Exact_beta", "Exact_se", "Exact_pval")

        BPIb <- select(BPi, "SNP", "INFO", "FRQ_A_25060", "Neff")
        BPIIb <- select(BPii, "SNP", "INFO", "FRQ_A_6781", "Neff")

        CCGWAS <- inner_join(CCGWAS_ldsr, BPIb, by="SNP")
        CCGWAS_LDSR <- inner_join(CCGWAS, BPIIb, by="SNP")

        #create INFO, FRQ and Neff columns that are averages of both
        CCGWAS_LDSR <- CCGWAS_LDSR %>%
            mutate(INFO = (INFO.x+INFO.y)/2) %>% 
            mutate(Neff = (Neff.x+Neff.y)/2) %>%
            mutate(FRQ = (FRQ_A_25060+FRQ_A_6781)/2) %>%
            mutate(OR = exp(Exact_beta))

        # select cols for LDSR
        ccgwas_LDSR <- CCGWAS_LDSR %>%
            select("SNP", "CHR", "BP", "EA", "NEA", "OR", "Exact_se", "Exact_pval", "INFO", "FRQ", "Neff") %>%
            rename(SE = Exact_se) %>%
            rename(P = Exact_pval)

        write.table(ccgwas_LDSR, "postGWAS_ldsr_input.txt", col.names = T, row.names = F, quote = F)
        quit()

        module load ldsc/1.0.1

        munge_sumstats.py --sumstats $LDSR_post/postGWAS_ldsr_input.txt --N-col Neff --snp SNP --a1 EA --a2 NEA --p P --signed-sumstats OR,1 --out $LDSR_post/LDSR_munge

        # liability scale 
        ldsc.py \
        --h2 $LDSR_post/LDSR_munge.sumstats.gz  \
        --ref-ld-chr ${LD}/eur_w_ld_chr/ \
        --w-ld-chr ${LD}/eur_w_ld_chr/ \
        --samp-prev 0.5 \
        --pop-prev 0.005 \
        --out $LDSR_post/CCGWAS_liab


    #LDAK
    #summary stats file needs to have the following column headers:
    #predictor (format Chr:BP)
    #A1
    #A2
    #n (total sample size)
    #Z (effect size/se)

    # restrict your summary statistics MAF > 0.05
        cd $LDAK_post 
        module load R/4.4.0
        R 
        library(dplyr)
        library(data.table)

        BPii <- fread(file = paste0((Sys.getenv("DATA")),"/BPii.txt"), header = T)
        BPi <- fread(file = paste0((Sys.getenv("DATA")),"/BPi.txt"), header = T)

        BPii_maf <- BPii %>%
            mutate(maf = (1-(FRQ_A_6781+FRQ_U_364075)/2))  
        BPii_maf_filt <- BPii_maf %>%
            filter(maf > 0.05)                            

        BPi_maf <- BPi %>%
            mutate(maf = (1-(FRQ_A_25060+FRQ_U_449978)/2)) 
        BPi_maf_filt <- BPi_maf %>%
            filter(maf > 0.05)                             

        # restrict CCGWAS output to SNPs with MAF > 0.5
        ccgwas <- fread(file = paste0((Sys.getenv("CCGWAS_out")),"/test_final.out.results.gz"), header = T)

        ccgwas_inner <- inner_join(ccgwas,BPi_maf_filt, by = "SNP") 
        ccgwas_inner_cols <- ccgwas_inner %>%
            select(SNP, CHR.x, BP.x, EA, NEA, Exact_beta, Exact_se)

        ccgwas_inner2 <- inner_join(ccgwas_inner_cols,BPii_maf_filt, by = "SNP") 
        ccgwas_maf <- ccgwas_inner2 %>%
            select(SNP, CHR.x, BP.x, EA, NEA, Exact_beta, Exact_se) 

        #calculate Z
        ccgwas_maf$Z <- ccgwas_maf$Exact_beta/ccgwas_maf$Exact_se
        median(ccgwas_maf$Z) # 0.04328165
        hist(dt$Z) # Output seems reasonable 

        #compute 'Predictor'
        ccgwas_maf <- ccgwas_maf %>%
            mutate(Predictor = paste(CHR.x, BP.x, sep = ":")) %>%
            arrange(CHR.x, BP.x)

        ##### remove MHC region #####  chr 6, 25-35Mb
        xMHC <- ccgwas_maf %>%
            filter(CHR.x != 6 & (BP.x < 25000000 | BP.x > 35000000))

        #select columns
        ldak <- xMHC %>%
            select(Predictor, EA, NEA, Z) %>%
            rename(A1 = EA, A2 = NEA)

        ##### check for duplicate positions #####
        duplicate_positions <- ldak %>%
            group_by(Predictor) %>%
            filter(n() > 1) %>%
            distinct(Predictor)  
            
            nrow(duplicate_positions) #660 rows 

        duplicate_positions_vector <- duplicate_positions$Predictor

        #filter out rows 
        filtered_data <- ldak %>%
            filter(!(ldak[[1]] %in% duplicate_positions_vector))

        # Save the filtered dataset back to a file
        write.table(filtered_data, "ccgwas_maf_gt_0.05.raw", sep = "\t", row.names = FALSE, quote = FALSE)
        quit()

    #run LDAK on home server
        /Users/oliviap/Documents/Masters/Dissertation/week5/LDAK
        chmod a+x ldak5.2.mac
        ./ldak5.2.mac

            # LDAK-thin tagging file
        ./ldak5.2.mac --sum-hers snpher2 --summary ccgwas_maf_gt_0.05.raw --tagfile ldak.thin.hapmap.gbr.tagging --fixed-n 792643 --check-sums NO
            # BLD-LDAK tagging file
        ./ldak5.2.mac --sum-hers snpher4 --summary ccgwas_maf_gt_0.05.raw --tagfile bld.ldak.hapmap.gbr.tagging --fixed-n 792643 --check-sums NO

## Estimate genetic correlations with other disorders using LDSR
    #MDD
    export MDD=/home/c.c23045409/dissertation/postGWAS/postGWAS_xtLDSR/MDD
    export MDD_scratch=/scratch/c.c23045409/dissertation/postGWAS/postGWAS_xtLDSR/MDD

        munge_sumstats.py --sumstats $MDD_scratch/MDD_sumstats.txt --snp SNP --a1 A1 --a2 A2 --p P --signed-sumstats OR,1 --out $MDD/MDD

    ldsc.py \
    --rg $LDSR_post/LDSR_munge.sumstats.gz,$MDD/MDD.sumstats.gz \
    --ref-ld-chr $LD/eur_w_ld_chr/ \
    --w-ld-chr $LD/eur_w_ld_chr/ \
    --samp-prev 0.5,0.318 \
    --pop-prev 0.005,0.15 \
    --out $MDD/BPvMDD 

    #SCZ
    export SCZ=/home/c.c23045409/dissertation/postGWAS/postGWAS_xtLDSR/SCZ
    export SCZ_scratch=/scratch/c.c23045409/dissertation/postGWAS/postGWAS_xtLDSR/SCZ

    munge_sumstats.py --sumstats $SCZ_scratch/SCZ_sumstats.txt --snp SNP --a1 A1 --a2 A2 --p PVAL --signed-sumstats OR,1 --out $SCZ/SCZ

    ldsc.py \
    --rg $LDSR_post/LDSR_munge.sumstats.gz,${SCZ}/SCZ.sumstats.gz \
    --ref-ld-chr ${LD}/eur_w_ld_chr/ \
    --w-ld-chr ${LD}/eur_w_ld_chr/ \
    --samp-prev 0.5,0.4086372 \
    --pop-prev 0.005,0.012 \
    --out $SCZ/BPvSCZ 

    #CDG
    export CDG=/home/c.c23045409/dissertation/postGWAS/postGWAS_xtLDSR/CDG
    export CDG_scratch=/scratch/c.c23045409/dissertation/postGWAS/postGWAS_xtLDSR/CDG

    munge_sumstats.py --sumstats $CDG_scratch/CDG_sumstats.txt --snp SNP --a1 A1 --a2 A2 --p PVAL --signed-sumstats OR,1 --out $CDG/CDG

    ldsc.py \
    --rg $LDSR_post/LDSR_munge.sumstats.gz,$CDG/CDG.sumstats.gz  \
    --ref-ld-chr ${LD}/eur_w_ld_chr/ \
    --w-ld-chr ${LD}/eur_w_ld_chr/ \
    --samp-prev 0.5,0.3203901 \
    --pop-prev 0.005,0.273 \
    --out $CDG/BPvCDG 

## sex-stratified genetic correlations
    # sumstats file format:SNP, CHR, POS, REF, ALT, FRQ, OR, SE, P, N 

    cd $sex_strat
    module load R/4.4.0
    R 

    #set up envitonment 
    library(dplyr)
    library(data.table) 

    #load files
    BPii_m <- fread(file = paste0((Sys.getenv("sex_strat")),"/daner_pgc4_males_BDII_only_newPC"), header = T)
    BPii_f <- fread(file = paste0((Sys.getenv("sex_strat")),"/daner_pgc4_females_BDII_only_newPC"), header = T)
    BPi_f <- fread(file = paste0((Sys.getenv("sex_strat")),"/daner_pgc4_females_BDI_only_newPC"), header = T)
    BPi_m <- fread(file = paste0((Sys.getenv("sex_strat")),"/daner_pgc4_males_BDI_only_newPC"), header = T)

    #select columns 
    BPii_mb <- select(BPii_m, "SNP", "CHR", "BP", "A1", "A2", "OR", "SE", "P", "Nca", "Nco")
    BPii_fb <- select(BPii_f, "SNP", "CHR", "BP", "A1", "A2", "OR", "SE", "P", "Nca", "Nco")
    BPi_mb <- select(BPi_m, "SNP", "CHR", "BP", "A1", "A2", "OR", "SE", "P", "Nca", "Nco")
    BPi_fb <- select(BPi_f, "SNP", "CHR", "BP", "A1", "A2", "OR", "SE", "P", "Nca", "Nco")

    #save new files 
    write.table(BPii_mb, "BPii_m_sumstats.txt", col.names = T, row.names = F, quote = F)
    write.table(BPii_fb, "BPii_f_sumstats.txt", col.names = T, row.names = F, quote = F)
    write.table(BPi_mb, "BPi_m_sumstats.txt", col.names = T, row.names = F, quote = F)
    write.table(BPi_fb, "BPi_f_sumstats.txt", col.names = T, row.names = F, quote = F)
    quit()

    module load ldsc/1.0.1

    #BPii_male
    munge_sumstats.py --sumstats $sex_strat/BPii_m_sumstats.txt --N-cas-col Nca --N-con-col Nco --snp SNP --a1 A1 --a2 A2 --p P --signed-sumstats OR,1 --out $sex_strat/BPii_m

    ldsc.py \
    --rg $LDSR_post/LDSR_munge.sumstats.gz,$sex_strat/BPii_m.sumstats.gz  \
    --ref-ld-chr ${LD}/eur_w_ld_chr/ \
    --w-ld-chr ${LD}/eur_w_ld_chr/ \
    --samp-prev 0.5,0.1024183 \
    --pop-prev 0.005,0.004 \
    --out $sex_strat/BPvBPii_m

    #BPii_female
    munge_sumstats.py --sumstats $sex_strat/BPii_f_sumstats.txt --N-cas-col Nca --N-con-col Nco --snp SNP --a1 A1 --a2 A2 --p P --signed-sumstats OR,1 --out $sex_strat/BPii_f

    ldsc.py \
    --rg $LDSR_post/LDSR_munge.sumstats.gz,$sex_strat/BPii_f.sumstats.gz  \
    --ref-ld-chr ${LD}/eur_w_ld_chr/ \
    --w-ld-chr ${LD}/eur_w_ld_chr/ \
    --samp-prev 0.5,0.1577917 \
    --pop-prev 0.005,0.004 \
    --out $sex_strat/BPvBPii_f

        #BPi_female
    munge_sumstats.py --sumstats $sex_strat/BPi_f_sumstats.txt --N-cas-col Nca --N-con-col Nco --snp SNP --a1 A1 --a2 A2 --p P --signed-sumstats OR,1 --out $sex_strat/BPi_f

        ldsc.py \
    --rg $LDSR_post/LDSR_munge.sumstats.gz,$sex_strat/BPi_f.sumstats.gz  \
    --ref-ld-chr ${LD}/eur_w_ld_chr/ \
    --w-ld-chr ${LD}/eur_w_ld_chr/ \
    --samp-prev 0.5,0.3096903 \
    --pop-prev 0.005,0.006 \
    --out $sex_strat/BPvBPi_f

    #BPi_male
    munge_sumstats.py --sumstats $sex_strat/BPi_m_sumstats.txt --N-cas-col Nca --N-con-col Nco --snp SNP --a1 A1 --a2 A2 --p P --signed-sumstats OR,1 --out $sex_strat/BPi_m

    ldsc.py \
    --rg $LDSR_post/LDSR_munge.sumstats.gz,$sex_strat/BPi_m.sumstats.gz  \
    --ref-ld-chr ${LD}/eur_w_ld_chr/ \
    --w-ld-chr ${LD}/eur_w_ld_chr/ \
    --samp-prev 0.5,0.2618788 \
    --pop-prev 0.005,0.006 \
    --out $sex_strat/BPvBPi_m