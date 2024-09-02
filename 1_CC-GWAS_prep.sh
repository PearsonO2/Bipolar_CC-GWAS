# Determining CC-GWAS input parameters

export DATA=/scratch/c.c23045409/dissertation/ccgwas_input/data 
export LDSR=/scratch/c.c23045409/dissertation/ccgwas_input/LDSR 
export LD=/scratch/c.c23045409/dissertation/ccgwas_input/LDSR/ref
export LDAK=/scratch/c.c23045409/dissertation/postGWAS/LDAK
export MIX=/home/c.c23045409/containers/mixer
export MIX_out=/home/c.c23045409/containers/mixer/c.c23045409/containers/mixer
export MIX_dir=/scratch/c.c23045409/dissertation/ccgwas_input/mixer 

## SNP-based heritability on the libability scale 
    ###LDSR
    module load ldsc/1.0.1

    # munge sumstats
    munge_sumstats.py --sumstats $DATA/BPi.txt --N-col Neff --snp SNP --a1 A1 --a2 A2 --p P --signed-sumstats OR,1 --out $LDSR/BPI_neff
    munge_sumstats.py --sumstats $DATA/BPii.txt --N-col Neff --snp SNP --a1 A1 --a2 A2 --p P --signed-sumstats OR,1 --out $LDSR/BPII_neff

        #BDI
        ldsc.py \
        --h2 $LDSR/BPI_neff.sumstats.gz \
        --ref-ld-chr $LD/eur_w_ld_chr/ \
        --w-ld-chr $LD/eur_w_ld_chr/ \
        --samp-prev 0.5 \
        --pop-prev 0.006 \
        --out $LDSR/BPI_Neff

        #BDII
        ldsc.py \
        --h2 $LDSR/BPII_neff.sumstats.gz \
        --ref-ld-chr $LD/eur_w_ld_chr/ \
        --w-ld-chr $LD/eur_w_ld_chr/ \
        --samp-prev 0.5 \
        --pop-prev 0.004 \
        --out $LDSR/BPII_Neff

    # cross trait LDSR

    ldsc.py \
        --rg $LDSR/BPI_neff.sumstats.gz,$LDSR/BPII_neff.sumstats.gz \
        --ref-ld-chr $LD/eur_w_ld_chr/ \
        --w-ld-chr $LD/eur_w_ld_chr/ \
        --samp-prev 0.5,0.5 \
        --pop-prev 0.006,0.004 \
        --out $LDSR/BPIvII_neff


        #calculating p-vales 
        module load R/4.4.0
        R 
        # Heritability estimate
        h2_est <- 0.116 
        # Standard error
        se <- 0.01
        null <- 0.00 #null = 1 for intercept

        # Calculate Z score
        z_score <- (h2_est - null) / se
        z_score

        # Log p-value for better precision
        log_p_value <- pnorm(z_score, lower.tail = FALSE, log.p = TRUE)
        exact_p_value <- exp(log_p_value)
        exact_p_value

    ### LDAK 
    #summary stats file needs to have the following column headers:
    #predictor (format Chr:BP)
    #A1
    #A2
    #n (total sample size)
    #Z (effect size/se)

    cd $LDAK

    # pre_LDAK prep in R
        module load R/4.4.0
        R 
        library(dplyr)
        library(data.table)
        library(R.utils)

        BPii <- fread(file = paste0((Sys.getenv("DATA")),"/BPii.txt"), header = T)
        BPi <- fread(file = paste0((Sys.getenv("DATA")),"/BPi.txt"), header = T)

    # restrict your summary statistics MAF > 0.05
        BPii_maf <- BPii %>%
            mutate(maf = (1-(FRQ_A_6781+FRQ_U_364075)/2)) #7188236 rows 

        BPii_maf_filt <- BPii_maf %>%
            filter(maf > 0.05)                            # 5418771 rows 

        BPi_maf <- BPi %>%
            mutate(maf = (1-(FRQ_A_25060+FRQ_U_449978)/2)) #7391594 rows 

        BPi_maf_filt <- BPi_maf %>%
            filter(maf > 0.05)                             #5431057 rows 

    #calculate Z
        ##BPi
        BPi_maf_filt$Z <- log(BPi_maf_filt$OR)/BPi_maf_filt$SE
        median(BPi_maf_filt$Z) # -0.0409147

        ##BPii
        BPii_maf_filt$Z <- log(BPii_maf_filt$OR)/BPii_maf_filt$SE
        median(BPii_maf_filt$Z) # -0.07114342

    #compute 'Predictor' column
        ##BPi
        BPi_pred <- BPi_maf_filt %>%
            mutate(Predictor = paste(CHR, BP, sep = ":")) %>%
            arrange(CHR, BP)

        #BPii
        BPii_pred <- BPii_maf_filt %>%
            mutate(Predictor = paste(CHR, BP, sep = ":")) %>%
            arrange(CHR, BP)

    ##### remove MHC region #####  chr 6, 25-35Mb
        ##BPI
        BPi_xMHC <- BPi_pred %>%
            filter(CHR != 6 & (BP < 25000000 | BP > 35000000))

        ##BPii
        BPii_xMHC <- BPii_pred %>%
            filter(CHR != 6 & (BP < 25000000 | BP > 35000000))

    #select columns
        ##BPi
        BPi_ldak <- BPi_xMHC %>%
            select(Predictor, A1, A2, Z)

        ##BPii
        BPii_ldak <- BPii_xMHC %>%
            select(Predictor, A1, A2, Z)

    ##### check for duplicate positions #####
        #BPi
        duplicate_positions_BPi <- BPi_ldak %>%
            group_by(Predictor) %>%
            filter(n() > 1) %>%
            distinct(Predictor)  #1832 rows 

        duplicate_positions_vector_BPi <- duplicate_positions_BPi$Predictor

        #filter out rows 
        BPi_filtered <- BPi_ldak %>%
            filter(!(BPi_ldak[[1]] %in% duplicate_positions_vector_BPi))

        ##BPii
        duplicate_positions_BPii <- BPii_ldak %>%
            group_by(Predictor) %>%
            filter(n() > 1) %>%
            distinct(Predictor) #1787 rows 
        
        duplicate_positions_vector_BPii <- duplicate_positions_BPii$Predictor

        #filter out rows 
        BPii_filtered <- BPii_ldak %>%
            filter(!(BPii_ldak[[1]] %in% duplicate_positions_vector_BPii))

    # Save the filtered dataset back to a file
        write.table(BPi_filtered, "BPi_maf_gt_0.05.raw", sep = "\t", row.names = FALSE, quote = FALSE)
        write.table(BPii_filtered, "BPii_maf_gt_0.05.raw", sep = "\t", row.names = FALSE, quote = FALSE)
    quit()

    ## running LDAK
        # cannot be run on VScode (run on terminal)
        cd $LDAK
        chmod a+x ldak5.2.mac
        ./ldak5.2.mac

        #LDAK-thin 
        ./ldak5.2.mac --sum-hers snpher-thin_BPI --summary pre-CCGWAS_LDAK/BPi_maf_gt_0.05.raw --tagfile ldak.thin.hapmap.gbr.tagging --fixed-n 792643 --check-sums NO
        ./ldak5.2.mac --sum-hers snpher-thin_BPII --summary pre-CCGWAS_LDAK/BPii_maf_gt_0.05.raw --tagfile ldak.thin.hapmap.gbr.tagging --fixed-n 792643 --check-sums NO

        #BDL-LDAK
        # using BLD-LDAK tagging file (It is a slightly more complex model of heritability which also takes into account functional annotations)
        ./ldak5.2.mac --sum-hers snpher-BLD_BPi --summary pre-CCGWAS_LDAK/BPi_maf_gt_0.05.raw --tagfile bld.ldak.hapmap.gbr.tagging --fixed-n 792643 --check-sums NO
        ./ldak5.2.mac --sum-hers snpher-BLD_BPii --summary pre-CCGWAS_LDAK/BPii_maf_gt_0.05.raw --tagfile bld.ldak.hapmap.gbr.tagging --fixed-n 792643 --check-sums NO


# estimating the number of independent loci 
    ### Mixer 
    cd $MIX 
    module load singularity 
    sbatch $MIX/mixer.job

        ###mixer.job
        #!/bin/bash
        #SBATCH --job-name=mix_BPivii
        #SBATCH --time=48:00:00
        #SBATCH --cpus-per-task=16
        #SBATCH --array=1-20
        #SBATCH --account=scw2173
        #SBATCH --partition compute
        #SBATCH --mem-per-cpu=8000M
        #SBATCH --mail-user=PearsonO2@cardiff.ac.uk
        #SBATCH --mail-type=BEGIN,END,FAIL
        #SBATCH --output=mixer_4.out
        #SBATCH --error=mixer_4.err

        cd /home/$USER/containers/mixer

        export COMORMENT=/home/$USER/containers
        export SINGULARITY_BIND=$COMORMENT/mixer/reference:/REF
        export SIF=$COMORMENT/mixer/singularity
        export MIXER_COMMON_ARGS="--ld-file /REF/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld --bim-file /REF/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.@..bim --threads 16"
        export PYTHON="singularity exec --home=$PWD:/home $SIF/mixer.sif python"
        export BPI=/scratch/c.c23045409/dissertation/ccgwas_input/LDSR/BPI.sumstats.gz
        export BPII=/scratch/c.c23045409/dissertation/ccgwas_input/LDSR/BPII.sumstats.gz
        export REP="rep${SLURM_ARRAY_TASK_ID}"
        export EXTRACT="--extract /REF/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.$REP.snps"

        #univariate analysis
        $PYTHON /tools/mixer/precimed/mixer.py fit1 $MIXER_COMMON_ARGS $EXTRACT --trait1-file $BPI --out BPi2.fit.$REP
        $PYTHON /tools/mixer/precimed/mixer.py fit1 $MIXER_COMMON_ARGS $EXTRACT --trait1-file $BPII --out BPii2.fit.$REP

        $PYTHON /tools/mixer/precimed/mixer.py test1 $MIXER_COMMON_ARGS --trait1-file $BPI --load-params /nfshome/store04/users/c.c23045409/containers/mixer/c.c23045409/containers/mixer/BPi2.fit.$REP.json --out BPI2.test.$REP
        $PYTHON /tools/mixer/precimed/mixer.py test1 $MIXER_COMMON_ARGS --trait1-file $BPII --load-params /nfshome/store04/users/c.c23045409/containers/mixer/c.c23045409/containers/mixer/BPii2.fit.$REP.json --out BPII2.test.$REP


#once job finished, create the final file
        cd $MIX_out
        module load singularity 
        singularity shell --bind /home/c.c23045409/containers/mixer:/mnt /home/c.c23045409/containers/mixer/singularity/mixer.sif

    # creating the file
    python /tools/mixer/precimed/mixer_figures.py combine --json $MIX_out/BPi2.fit.rep@.json  --out $MIX_out/BPi.fit
    python /tools/mixer/precimed/mixer_figures.py combine --json $MIX_out/BPii2.fit.rep@.json  --out $MIX_out/BPii.fit 
    python /tools/mixer/precimed/mixer_figures.py combine --json $MIX_out/BPivBPii.fit.rep@.json  --out $MIX_out/BPivBPii.fit
    python /tools/mixer/precimed/mixer_figures.py one --json $MIX_out/BPi.fit.json $MIX_out/BPii.fit.json --out $MIX_dir/BPi_and_BPii.fit --trait1 BPi BPii --statistic mean std --ext png 
