
## Comparative Genetic Analysis of bipolar type 1 and type 2

To make full use of the suggested environment variables, set out own directory as such: 
- ccgwas_input
	- data (original GWAS summary stats)
	- LDAK
	- LDSR
	- mixer 
	- ref
- ccgwas_output
- post_ccgwas
	- annot
	- LDAK
	- LDSR
	- LD_clumping
	- FINEMAP
		- ref
	- postGWAS_xtLDSR
		- MDD
		- SCZ
		- CDG
	- sex_strat
	- SMR
		- ref 
	- PRS
		- DATA (for base data)
		- test (for test data)
		- dir
		- LDSR
		- mixer
		- PRScs
		- QC
			- basedata
			- targetdata
		- Ranalysis
		- REF


1) CC-GWAS_prep.sh
	- This is for determining the input parameters for the CC-GWAS. 
	- software required: 
		- LDSR ([https://github.com/bulik/ldsc](https://github.com/bulik/ldsc))
		- LDAK SumHer ([https://dougspeed.com/sumher/](https://dougspeed.com/sumher/))
		- Mixer ([https://github.com/precimed/mixer](https://github.com/precimed/mixer))
		- R studio (loaded into the shell environment)

2) CC-GWAS.sh
	- CC-GWAS is an R package but is loaded into the shell environment
		- [https://github.com/wouterpeyrot/CCGWAS](https://github.com/wouterpeyrot/CCGWAS)

3) interrogate_CCGWAS_output.R
	- an R script for lambda calculations, QQ and Manhattan plots. 
	- can be run in /ccgwas_output

4) post_CC-GWAS.sh
	- This is for annotating any significant loci using PLINK and SMR, LD clumping and fine mapping
	- Software required: 
		- SMR: [https://yanglab.westlake.edu.cn/software/smr/#Overview](https://yanglab.westlake.edu.cn/software/smr/)
		- FINEMAP: [http://www.christianbenner.com](http://www.christianbenner.com/)
		- PLINK: [https://www.cog-genomics.org/plink/1.9/](https://www.cog-genomics.org/plink/1.9/)

5) postCC-GWAS_H2_genetic_correlations.sh
	- This is for the cross trait and sex-stratified genetic correlations

6) preCCGWAS_PRS.sh
	- same software required as file 1

7) PRS_CC-GWAS.sh
	- preparing base data files for CC-GWAS in R 
	- performing CC-GWAS in R

8) interrogate_PRS_CC-GWAS_output.R
	- generate QQ plots, lambda values in R

9) post_PRS_CC-GWAS_QC.sh
	- prepare PRS CC-GWAS file for PRS analysis in R
	- check heritability (LDSR)
	- MAF > 0.01 and INFO > 0.8
	- remove duplicate SNPs
	- remove ambiguous SNPs
	- remove SNPs in extended MHC locus
	- as per https://choishingwan.github.io/PRS-Tutorial/

10) BDRN_QC.sh
	- preparing BDRN sample for PRS 
		eg adding Sex info
	- LD prunining (PLINK)
	- interrogate heterozygosity
	- recode mismatching SNPs 
	- sex check
	- remove individuals with relatedness f coefficient greater than 0.125
	- as per https://choishingwan.github.io/PRS-Tutorial/

11) PRS.sh
	- generate 10 PCs from QC'd base data
	- create the final file for PRS 
	- run PRScs 
	- extract the scores using R 
	- software: [https://github.com/getian107/PRScs](https://github.com/getian107/PRScs)

12) generating_phenofile_for_regression.R
	- adding the PCs and phenotypes to PRS scores 

13) PRS_regression.R
	- plot eigenvalues
	- scale PRS values 
	- for running the regressions, comparing the models and obtaining the beta, p-value, standard errors and confidence intervals
