#interrogating CC-GWAS output 


# prepare R environment
library(dplyr)
library(data.table)
library(qqman)

CCGWAS <- fread("test_final.out.results.gz")


##QQ plot

CCGWAS_QQ <- CCGWAS %>%
  select(SNP, CHR, BP, EA, OLS_beta, OLS_se, OLS_pval, Exact_beta, Exact_se, Exact_pval)

#OLS values 
tiff("QQ_plot_OLS.tiff", width = 10, height = 9, units = 'cm', res = 300)
qq(CCGWAS_QQ$OLS_pval)
dev.off()

#Exact values 
tiff("QQ_plot_Exact.tiff", width = 10, height = 9, units = 'cm', res = 300)
qq(CCGWAS_QQ$Exact_pval)
dev.off()

## calculation of lambda and lambda_1000 
cases<-31841 #combined 6781+25060
controls<-760802 #364075+449978-53251
gclambda<-1

#OLS
CCGWAS_QQ$CHISQ_OLS<-qchisq(1-CCGWAS_QQ$OLS_pval,1)
CCGWAS_QQ$CHISQ_OLS<-CCGWAS_QQ$CHISQ_OLS/gclambda

lambda_OLS<-median(CCGWAS_QQ$CHISQ_OLS)/0.456

lambda1000_OLS<-1+(lambda_OLS-1)*(1/cases +1/controls)/(1/1000+1/1000) #UCL BD vs New Half of 9999 controls
print(lambda_OLS) #1.144684
print(lambda1000_OLS) #1.002367

#Exact 
CCGWAS_QQ$CHISQ_Exact<-qchisq(1-CCGWAS_QQ$Exact_pval,1)
CCGWAS_QQ$CHISQ_Exact<-CCGWAS_QQ$CHISQ_Exact/gclambda

lambda_Exact<-median(CCGWAS_QQ$CHISQ_Exact)/0.456

lambda1000_Exact<-1+(lambda_Exact-1)*(1/cases +1/controls)/(1/1000+1/1000) #UCL BD vs New Half of 9999 controls
print(lambda_Exact) #0.9120177
print(lambda1000_Exact) #0.9985606

##manhattan plot 

#OLS
tiff("ManhattanOLS.tiff", width = 20, height = 8, units = 'cm', res = 300)
manhattan(CCGWAS, p = "OLS_pval", suggestiveline = -log10(5e-8),
          col = c("blue", "lightblue"))
dev.off()

#Exact
tiff("ManhattanExact.tiff", width = 20, height = 8, units = 'cm', res = 300)
manhattan(CCGWAS, p = "Exact_pval", suggestiveline = -log10(1e-4),
          col = c("blue", "lightblue"))
dev.off()


