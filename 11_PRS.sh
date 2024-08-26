#PRS_prep 

export PRS_DATA=/scratch/c.c23045409/dissertation/postGWAS/PRS/DATA
export PRS_dir=/scratch/c.c23045409/dissertation/postGWAS/PRS/dir
export BDRN=/scratch/c.c23045409/dissertation/postGWAS/PRS/test/raw
LD=/scratch/c.mpmlh/MET588_h2_rg/MET588_LSH/
export PRS=/scratch/c.c23045409/dissertation/postGWAS/PRS
PRScs=/scratch/c.c23045409/dissertation/postGWAS/PRS/PRScs
export PRS_LDSR=/scratch/c.c23045409/dissertation/postGWAS/PRS/LDSR
export QC=/scratch/c.c23045409/dissertation/postGWAS/PRS/QC/targetdata
export PRS_R=/scratch/c.c23045409/dissertation/postGWAS/PRS/Ranalysis
export PRS_res=/scratch/c.c23045409/dissertation/postGWAS/PRS/dir


##########################################################################################
# generating pcs
cd $PRS

module purge
module load plink/2.0

plink2 --bfile ${QC}/M_BDRN.qc --pca 10 --out ${QC}/M_BDRN.PCA

##########################################################################################
# create sumstats file 

  #GWAS summary statistics file must have the following column headers: SNP, A1, A1, OR/BETA, SE
    gunzip $PRS_DATA/PRS.out.results.gz
    awk 'NR==1 {print $1, "A1", "A2", "BETA", "SE"} NR>1 {print $1, $4, $5, $6, $7}' $PRS_DATA/PRS.out.results > sumstats.txt
    gzip $PRS_DATA/PRS.out.results


##########################################################################################
cd $PRS
###PRS.job

#SBATCH --job-name=PRS_CCGWAS
#SBATCH --array=1-22
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --array=1-20
#SBATCH --account=scw2173
#SBATCH --partition compute 
#SBATCH --mem-per-cpu=8000M
#SBATCH --mail-user=PearsonO2@cardiff.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=PRS.out
#SBATCH --error=PRS.err

module load anaconda/2020.02
cd /scratch/c.c23045409/dissertation/postGWAS/PRS

export REF=/scratch/c.c23045409/dissertation/postGWAS/PRS/REF
export DATA=/scratch/c.c23045409/dissertation/postGWAS/PRS/DATA
export TEST=/scratch/c.c23045409/dissertation/postGWAS/PRS/test
export DIR=/scratch/c.c23045409/dissertation/postGWAS/PRS/dir
export QC=/scratch/c.c23045409/dissertation/postGWAS/PRS/QC/targetdata
export PYTHON="python PRScs/PRScs.py"

$PYTHON \
  --ref_dir=$REF/ldblk_1kg_eur \
  --bim_prefix=$QC/M_BDRN.qc \
  --sst_file=$DATA/sumstats.txt \
  --n_gwas=207365 \
  --n_iter=100 --n_burnin=50 \
  --out_dir=$DIR/cs

module purge
module load plink/1.9

plink \
 --bfile $QC/M_BDRN.qc \
 --maf 0.1 \ 
 --score $DIR/cs_pst_eff_a1_b0.5_phiauto_chr${SLURM_ARRAY_TASK_ID}.txt 2 4 6 sum \
 --out $DIR/chr${SLURM_ARRAY_TASK_ID}

##########################################################################################

# extract PRS scores using R

cd $PRS_R

module purge 
module load R/4.4.0
R 

#setup 
library(readr)
library(tidyverse)
library(datawizard)
library(ggplot2)
library(performance)
library(gtools)

#list all .profile files 
list_file <- list.files(path = Sys.getenv("PRS_res"), pattern = "*.profile", full.name=TRUE) %>% 
  gtools::mixedsort()

# create function to sum across chromosomes 
analyze <- function(threshold) {
  #prepare
  dat <- lapply(threshold,read.table, header = TRUE, sep= "")
  my_cols <- c("FID", "SCORESUM")
  dat <- lapply(dat, "[", , my_cols)
  #sum
  sum <- bind_rows(dat) %>%
    group_by(FID) %>%
    summarise(ALL_CHR_SUM = sum(SCORESUM))
}

# apply function to the list of *.profile files
prs_cs <- analyze(list_file)

# rename columns so nice and clean
colnames(prs_cs) <- c("IID", "PRScs")

write.table(prs_cs, "sumPRScs.txt", col.names=TRUE, row.names=FALSE, quote=FALSE)
quit()



