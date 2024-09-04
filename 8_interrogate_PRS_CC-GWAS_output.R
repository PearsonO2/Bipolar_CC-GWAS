# interrogate PRS CC-GWAS output

#prepare R environment
library(dplyr)
library(data.table)
library(qqman)

# read in CC-GWAS results 
PRS_CCGWAS <- fread("../PRS.out.results.gz")

##QQ plot

PRS_CCGWAS_QQ <- PRS_CCGWAS %>%
  select(SNP, CHR, BP, EA, OLS_beta, OLS_se, OLS_pval, Exact_beta, Exact_se, Exact_pval)

#OLS values 
tiff("QQ_plot_OLS_PRS.tiff", width = 10, height = 9, units = 'cm', res = 300)
qq(PRS_CCGWAS_QQ$OLS_pval)
dev.off()

#Exact values 
tiff("QQ_plot_Exact_PRS.tiff", width = 10, height = 9, units = 'cm', res = 300)
qq(PRS_CCGWAS_QQ$Exact_pval)
dev.off()

## calculation of lambda and lambda_1000 
cases<-24086 #combined 19578+4508
controls<-183279 #155439+73009-45169
gclambda<-1

#OLS
PRS_CCGWAS_QQ$CHISQ_OLS<-qchisq(1-PRS_CCGWAS_QQ$OLS_pval,1)
PRS_CCGWAS_QQ$CHISQ_OLS<-PRS_CCGWAS_QQ$CHISQ_OLS/gclambda

lambda_OLS<-median(PRS_CCGWAS_QQ$CHISQ_OLS)/0.456

lambda1000_OLS<-1+(lambda_OLS-1)*(1/cases +1/controls)/(1/1000+1/1000) #UCL BD vs New Half of 9999 controls
print(lambda_OLS) 
print(lambda1000_OLS) 

#Exact 
PRS_CCGWAS_QQ$CHISQ_Exact<-qchisq(1-PRS_CCGWAS_QQ$Exact_pval,1)
PRS_CCGWAS_QQ$CHISQ_Exact<-PRS_CCGWAS_QQ$CHISQ_Exact/gclambda

lambda_Exact<-median(PRS_CCGWAS_QQ$CHISQ_Exact)/0.456

lambda1000_Exact<-1+(lambda_Exact-1)*(1/cases +1/controls)/(1/1000+1/1000) #UCL BD vs New Half of 9999 controls
print(lambda_Exact)
print(lambda1000_Exact) 


