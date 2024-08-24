# pre CC-GWAS for PRS 

export PRS_DATA=/scratch/c.c23045409/dissertation/postGWAS/PRS/DATA
export PRS_dir=/scratch/c.c23045409/dissertation/postGWAS/PRS/dir
export BDRN=/scratch/c.c23045409/dissertation/postGWAS/PRS/test
LD=/scratch/c.mpmlh/MET588_h2_rg/MET588_LSH/

#mixer 
     # getting sumsets into correct format
    munge_sumstats.py --sumstats $PRS_DATA/daner_pgc3_BDI_noBDRN --N-cas-col Nca --N-con-col Nco --snp SNP --a1 A1 --a2 A2 --p P --signed-sumstats OR,1 --out $PRS_dir/BPI
    munge_sumstats.py --sumstats $$PRS_DATA/daner_pgc3_BDII_noBDRN --N-cas-col Nca --N-con-col Nco --snp SNP --a1 A1 --a2 A2 --p P --signed-sumstats OR,1 --out $PRS_dir/BPII

    cd /home/$USER/containers/mixer
    module load singularity 
    sbatch mixerPRS.job

    #making the images
    cd /home/c.c23045409/containers/mixer/c.c23045409/containers/mixer
    module load singularity 
    singularity shell --bind /home/c.c23045409/containers/mixer:/mnt /home/c.c23045409/containers/mixer/singularity/mixer.sif

    # creating the figures
    python /tools/mixer/precimed/mixer_figures.py combine --json BPiPRS.fit.rep@.json  --out BPiPRS.fit
    python /tools/mixer/precimed/mixer_figures.py combine --json BPiiPRS.fit.rep@.json  --out BPiiPRS.fit 
    python /tools/mixer/precimed/mixer_figures.py one --json BPiPRS.fit.json BPiiPRS.fit.json --out BPi_and_BPiiPRS.fit --trait1 BPi BPii --statistic mean std --ext png #univariate
    exit 

    #ultimatelty original number will be used for CCGWAS as sample withouht BRDN is too poorly powered for an accurate assessement by mixer 

# heritability estimates on the liability scale and genetic correlation
    module load ldsc/1.0.1

    #add Neff column (Neff_half*2) using awk 
    awk 'NR==1 {print $0, "Neff"; next} {print $0, $NF * 2}' $PRS_DATA/daner_pgc3_BDI_noBDRN > $PRS_DATA/daner_pgc3_BDI_noBDRN_Neff
    awk 'NR==1 {print $0, "Neff"; next} {print $0, $NF * 2}' $PRS_DATA/daner_pgc3_BDII_noBDRN > $PRS_DATA/daner_pgc3_BDII_noBDRN_Neff

    # #add Neff column (Neff_half*2) using R 
    R
    library(data.table)
    library(tidyverse)

    BPI <- fread(file = paste0((Sys.getenv("PRS_DATA")),"/daner_pgc3_BDI_noBDRN"), fill=TRUE)
    BPII <- fread(file = paste0((Sys.getenv("PRS_DATA")),"/daner_pgc3_BDII_noBDRN"), fill=TRUE)

    head(BPI)
    head(BPII)

    BPI_neff <- BPI %>%
    mutate(Neff = Neff_half * 2) 

    BPII_neff <- BPII %>%
    mutate(Neff = Neff_half * 2)

    write.table(BPI_neff, "daner_pgc3_BDI_noBDRN_Neff_R", col.names = T, row.names = F, quote = F)
    write.table(BPII_neff, "daner_pgc3_BDII_noBDRN_Neff_R", col.names = T, row.names = F, quote = F)


    # getting sumsets into correct format
        munge_sumstats.py --sumstats $PRS_DATA/daner_pgc3_BDI_noBDRN_Neff --N-col Neff --snp SNP --a1 A1 --a2 A2 --p P --signed-sumstats OR,1 --out $PRS_dir/BPI_prs
        munge_sumstats.py --sumstats $PRS_DATA/daner_pgc3_BDII_noBDRN_Neff --N-col Neff --snp SNP --a1 A1 --a2 A2 --p P --signed-sumstats OR,1 --out $PRS_dir/BPII_prs

        #checking R files 
        munge_sumstats.py --sumstats $PRS_DATA/daner_pgc3_BDI_noBDRN_Neff_R --N-col Neff --snp SNP --a1 A1 --a2 A2 --p P --signed-sumstats OR,1 --out $PRS_dir/BPI_prs_R
        munge_sumstats.py --sumstats $PRS_DATA/daner_pgc3_BDII_noBDRN_Neff_R --N-col Neff --snp SNP --a1 A1 --a2 A2 --p P --signed-sumstats OR,1 --out $PRS_dir/BPII_prs_R

        #sample prevalance estimate
        ##BPI = 19578/155439 = 0.1259529
        ##BPII = 4508/73009 =  0.06174581

        #population prevalence estimate
        ##BPI = 0.6%
        ##BPII = 0.4% 

    ldsc.py \
        --rg $PRS_dir/BPI_prs_R.sumstats.gz,$PRS_dir/BPII_prs_R.sumstats.gz \
        --ref-ld-chr $LD/eur_w_ld_chr/ \
        --w-ld-chr $LD/eur_w_ld_chr/ \
        --samp-prev 0.5,0.5 \
        --pop-prev 0.006,0.004 \
        --out $PRS_dir/BPIvII_prs_R

#Heritability of phenotype 1
#---------------------------
#Total Liability scale h2: 0.1687 (0.0091)
#Lambda GC: 1.3169
#Mean Chi^2: 1.3817
#Intercept: 1.0337 (0.011)
#Ratio: 0.0882 (0.0289)

#Heritability of phenotype 2/2
#-----------------------------
#Total Liability scale h2: 0.0693 (0.0168)
#Lambda GC: 1.0926
#Mean Chi^2: 1.0889
#Intercept: 1.0429 (0.0083)
#Ratio: 0.4825 (0.0936)

#Genetic Covariance
#------------------
#Total Liability scale gencov: 0.1035 (0.0092)
#Mean z1*z2: 0.2917
#Intercept: 0.171 (0.0068)

#Genetic Correlation
#-------------------
#Genetic Correlation: 0.9571 (0.1148)
#Z-score: 8.3353
#P: 7.7329e-17