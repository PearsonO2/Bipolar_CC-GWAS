# post PRS CC-GWAS QC

#as per https://choishingwan.github.io/PRS-Tutorial/
    # from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7612115/

export PRS_DATA=/scratch/c.c23045409/dissertation/postGWAS/PRS/DATA
export PRS_dir=/scratch/c.c23045409/dissertation/postGWAS/PRS/dir
export BDRN=/scratch/c.c23045409/dissertation/postGWAS/PRS/test
LD=/scratch/c.mpmlh/MET588_h2_rg/MET588_LSH/
export PRS=/scratch/c.c23045409/dissertation/postGWAS/PRS
PRScs=/scratch/c.c23045409/dissertation/postGWAS/PRS/PRScs
export PRS_LDSR=/scratch/c.c23045409/dissertation/postGWAS/PRS/LDSR
export BD_QC=/scratch/c.c23045409/dissertation/postGWAS/PRS/QC/basedata

# all on GWAS dataset 
#Heritablity check 
    module load ldsc/1.0.1
    module load R/4.4.0

    # getting sumsets into correct format 
    cd $PRS_LDSR
    R 
    library(dplyr)
    library(data.table)

    #need a file with SNP, CHR, POS, REF, ALT, FRQ, INFO, BETA, SE, P, LOG10P, N 

    BPI <- fread(file = paste0((Sys.getenv("PRS_DATA")),"/daner_pgc3_BDI_noBDRN"), fill=TRUE)
    BPII <- fread(file = paste0((Sys.getenv("PRS_DATA")),"/daner_pgc3_BDII_noBDRN"), fill=TRUE)
    CCGWAS <- fread(file = paste0((Sys.getenv("PRS_DATA")),"/PRS.out.results.gz"), fill=TRUE)

    CCGWAS_ldsr <- CCGWAS %>%
        select("SNP", "CHR", "BP", "EA", "NEA", "Exact_beta", "Exact_se", "Exact_pval")


    BPIb <- select(BPI, "SNP", "INFO", "FRQ_A_19578", "Neff_half")
    BPIIb <- select(BPII, "SNP", "INFO", "FRQ_A_4508", "Neff_half")
    
    CCGWAS <- inner_join(CCGWAS_ldsr, BPIb, by="SNP")
    CCGWAS_LDSR <- inner_join(CCGWAS, BPIIb, by="SNP")

    #create an INFO column for average of both
    CCGWAS_LDSR <- CCGWAS_LDSR %>%
        mutate(INFO = (INFO.x+INFO.y)/2) %>%
        mutate(FRQ = (FRQ_A_19578+FRQ_A_4508)/2) %>%
        mutate(Neff = Neff_half.x+Neff_half.y) %>%
        mutate(OR = exp(Exact_beta))

    # select cols for LDSR
    ccgwas_LDSR <- CCGWAS_LDSR %>%
        select("SNP", "CHR", "BP", "EA", "NEA", "OR", "Exact_se", "Exact_pval", "INFO", "FRQ", "Neff") %>%
        rename(SE = Exact_se) %>%
        rename(P = Exact_pval)

    write.table(ccgwas_LDSR, "PRS_ldsr_input.txt", col.names = T, row.names = F, quote = F)
    quit()

    munge_sumstats.py --sumstats $PRS_LDSR/ldsr_input.txt --N-col Neff --snp SNP --a1 EA --a2 NEA --p P --signed-sumstats OR,1 --out $PRS_LDSR/PRS
    
    # heritability estimates on the liable scale 
    ldsc.py \
    --h2 $PRS_LDSR/PRS.sumstats.gz  \
    --ref-ld-chr $LD/eur_w_ld_chr/ \
    --w-ld-chr $LD/eur_w_ld_chr/ \
    --samp-prev 0.5 \
    --pop-prev 0.005 \
    --out $PRS_LDSR/PRSccgwas_LDSR

        #SNP-h2 == 0.0015 (0.0011)
        #should be > 0.05 

# select for MAF>0.01 and INFO>0.8
    cd $PRS_LDSR
    #create a maf column for BPI and BPII sumstats without BDRN
    awk 'NR==1 {print $2, "MAF"} NR>1 {n_a = 4508; n_u = 73009; overall_freq = ($6 * n_a + $7 * n_u) / (n_a + n_u); maf = (overall_freq < 0.5) ? overall_freq : 1 - overall_freq; print $2, maf;}' $PRS_DATA/daner_pgc3_BDII_noBDRN > $BD_QC/BDII_maf.txt
    awk 'NR==1 {print $2, "MAF"} NR>1 {n_a = 19578; n_u = 155439; overall_freq = ($6 * n_a + $7 * n_u) / (n_a + n_u); maf = (overall_freq < 0.5) ? overall_freq : 1 - overall_freq; print $2, maf;}' $PRS_DATA/daner_pgc3_BDI_noBDRN > $BD_QC/BDI_maf.txt
    #add MAF cols to ldsr_input.txt
    awk 'NR==FNR {maf[$1]=$2; next} NR==1 {print $0, "MAF_BPi"} NR>1 {print $0, maf[$1]}' $BD_QC/BDI_maf.txt $PRS_LDSR/ldsr_input.txt > $BD_QC/prs_ccgwas_mafi.txt
    awk 'NR==FNR {maf[$1]=$2; next} NR==1 {print $0, "MAF_BPii"} NR>1 {print $0, maf[$1]}' $BD_QC/BDII_maf.txt $BD_QC/prs_ccgwas_mafi.txt > $BD_QC/prs_ccgwas_mafi_ii.txt

    #calucalte MAF_avg
    awk 'NR==1 {print $0, "MAF_avg"} NR>1 {maf_avg=($12+$13)/2; print $0, maf_avg}' $BD_QC/prs_ccgwas_mafi_ii.txt > $BD_QC/prs_ccgwas_maf.txt
    
    #select SNPs with MAF_avg > 0.01 and INFO > 0.8
    awk 'NR==1 || ($14 > 0.01) && ($9 > 0.8) {print}' $BD_QC/prs_ccgwas_maf.txt > $BD_QC/prs_ccgwas_maf_info.txt
    # Removed 74134 SNPs (4911919-4837785)

# duplicate SNPs
    awk '{seen[$3]++; if(seen[$3]==1){ print}}' $BD_QC/prs_ccgwas_maf_info.txt > $BD_QC/prs_gwas_maf_info_xdup.txt
    #removed 62004 (4837785-4775781)

#ambiguous SNPs
    awk '!( ($4=="A" && $5=="T") || ($4=="T" && $5=="A") || ($4=="G" && $5=="C") || ($4=="C" && $5=="G")) {print}' $BD_QC/prs_gwas_maf_info_xdup.txt > $BD_QC/prs_gwas_maf_info_xdup_amb.txt
    #no ambiguous SNPs

# remove SNPs in extended MHC locus of CHR6 25MB-35MB
    awk '!(($2 == 6) && ($3 >= 25000000) && ($3 <= 35000000))' $BD_QC/prs_gwas_maf_info_xdup_amb.txt > $BD_QC/prs_gwas_maf_info_xdup_amb_xMHC.txt
    # 26033 SNPs removed 
