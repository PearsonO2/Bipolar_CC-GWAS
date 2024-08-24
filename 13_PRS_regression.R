# perfrom regression analysis on PRS score and multiple variables

library(dplyr)
library(data.table)
library(MASS)
library(ggplot2)
library(performance)

#loading data
setwd("~/Documents/Masters/Dissertation/code /final_code")
PRS <- fread("PRS_complete.txt")
eigenval <- fread("../../week8/M_BDRN.PCA.eigenval")

#scaling PRS
PRS$PRSadj <- scale(PRS$PRScs)

### scree plot ###
eigenval <- fread("../../week8/M_BDRN.PCA.eigenval")
eigenval <- eigenval %>%
  mutate(Principal_Component = rownames(eigenval))

# Convert the Principal_Component column to a numeric type if it's not already
eigenval$Principal_Component <- as.numeric(eigenval$Principal_Component)

# Plot the scree plot
ggplot(eigenval, aes(x = Principal_Component, y = V1)) +
  geom_point(size = 2) +
  geom_line() +
  labs(title = "Scree Plot",
       x = "Principal Component",
       y = "Eigenvalue") +
  theme_minimal()

############################################################################
#all phenotypes are ordinal except for age at interview


####### no episodes mania LE as a continuous variable ########
no_mania <- PRS %>%
  filter(Number_of_episodes_mania_LE !=999) # lost 390 individuals 

no_mania$Number_of_episodes_mania_LE <- factor(no_mania$Number_of_episodes_mania_LE, ordered = TRUE)

# as an ordinal variable 
m1a_ordinal <- polr(Number_of_episodes_mania_LE ~ PRSadj + Age_at_Interview + Sex + ARRAY + PC1+ PC2+ PC3+ PC4+ PC5, data = no_mania, Hess = TRUE)
#without PRS
m1b_ordinal <- polr(Number_of_episodes_mania_LE ~ Age_at_Interview + Sex++ ARRAY + PC1+ PC2+ PC3+ PC4+ PC5, data = no_mania, Hess = TRUE)


compare_performance(m1a_ordinal, m1b_ordinal)

#obtaining the p-value 
summary(m1a_ordinal) #tval = 1.0339
p_values <- 2 * (1 - pnorm(abs(-1.0339)))
print(p_values)

####### no episodes depression LE #######
no_depression <- PRS %>%
  filter(Number_of_episodes_depression_LE !=999) # lost 460 individuals

no_depression$Number_of_episodes_depression_LE <- factor(no_depression$Number_of_episodes_depression_LE, ordered = TRUE)

#with PRS
m2a_ordinal <- polr(Number_of_episodes_depression_LE ~ PRSadj + Age_at_Interview + Sex + ARRAY + PC1 + PC2 + PC3 + PC4 + PC5, data = no_depression, Hess = TRUE)

# without PRS 
m2b_ordinal <- polr(Number_of_episodes_depression_LE ~ Age_at_Interview + Sex + ARRAY + PC1 + PC2 + PC3 + PC4 + PC5, data = no_depression, Hess = TRUE)

compare_performance(m2a_ordinal, m2b_ordinal)

#obtaining beta and pval
summary(m2a_ordinal)
p_values <- 2 * (1 - pnorm(abs(-3.00651)))
print(p_values)

####### BADDS psychosis as a count variable #######
BADDS_P <- PRS %>%
  filter(BADDS_P !=999) # lost 351 individuals 

BADDS_P$BADDS_P <- factor(BADDS_P$BADDS_P, ordered = TRUE)

# with PRS
m3a_ordinal <- polr(BADDS_P ~ PRSadj + Age_at_Interview + Sex + ARRAY + PC1 + PC2 + PC3 + PC4 + PC5, data = BADDS_P, Hess = TRUE)

# without PRS
m3b_ordinal <- polr(BADDS_P ~ Age_at_Interview + Sex + ARRAY + PC1 + PC2 + PC3 + PC4 + PC5, data = BADDS_P, Hess = TRUE)

compare_performance(m3b_ordinal, m3a_ordinal)

#obtaining beta and pval
summary(m3a_ordinal)
p_values <- 2 * (1 - pnorm(abs(4.256)))
print(p_values)

####### BADDS mood incongruent #######
BADDS_I <- PRS %>%
  filter(BADDS_I !=999) # lost 1797 individuals 

BADDS_I$BADDS_I <- factor(BADDS_I$BADDS_I, ordered = TRUE)

# with PRS
m4a_ordinal <- polr(BADDS_I ~ PRSadj + Age_at_Interview + Sex + ARRAY + PC1 + PC2 + PC3 + PC4 + PC5, data = BADDS_I, Hess = TRUE)

# without PRS
m4b_ordinal <- polr(BADDS_I ~ Age_at_Interview + Sex + ARRAY + PC1 + PC2 + PC3 + PC4 + PC5, data = BADDS_I, Hess = TRUE)

compare_performance(m4a_ordinal, m4b_ordinal)

#obtaining beta and pval
summary(m4a_ordinal)
p_values <- 2 * (1 - pnorm(abs(0.62116)))
print(p_values)

####### BADDS mania ####### 
BADDS_M <- PRS %>%
  filter(BADDS_M != 999) # lost 290 individuals 

BADDS_M$BADDS_M <- factor(BADDS_M$BADDS_M, ordered = TRUE)

# with PRS
m5a_ordinal <- polr(BADDS_M ~ PRSadj + Age_at_Interview + Sex + ARRAY + PC1 + PC2 + PC3 + PC4 + PC5, data = BADDS_M, Hess = TRUE)

#without PRS
m5b_ordinal <- polr(BADDS_M ~ Age_at_Interview + Sex + ARRAY + PC1 + PC2 + PC3 + PC4 + PC5, data = BADDS_M, Hess = TRUE)

compare_performance(m5a_ordinal, m5b_ordinal)

#obtaining beta and pval
summary(m5a_ordinal)
p_values <- 2 * (1 - pnorm(abs(7.6829)))
print(p_values)

####### BADDS depression as a continuous variable ####### 
BADDS_D <- PRS %>%
  filter(BADDS_D != 999) # lost 418 individuals

BADDS_D$BADDS_D <- factor(BADDS_D$BADDS_D, ordered = TRUE)

# with PRS
m6a_ordinal <- polr(BADDS_D ~ PRSadj + Age_at_Interview + Sex + ARRAY + PC1 + PC2 + PC3 + PC4 + PC5, data = BADDS_D, Hess = TRUE)

#without PRS
m6b_ordinal <- polr(BADDS_D ~ Age_at_Interview + Sex + ARRAY + PC1 + PC2 + PC3 + PC4 + PC5, data = BADDS_D, Hess = TRUE)

compare_performance(m6a_ordinal, m6b_ordinal)

#obtaining beta and pval
summary(m6a_ordinal)
p_values <- 2 * (1 - pnorm(abs(-3.1631)))
print(p_values)

####### Age of onset of impairment as a continuous variable ####### 
age_onset_impairment <- PRS %>%
  filter(Age_onset_impairment != 999) # 190 lost 

model7a <- lm(Age_onset_impairment ~ PRSadj + Age_at_Interview + Sex + ARRAY+ PC1+ PC2+ PC3+ PC4+ PC5, data = age_onset_impairment)
summary(model7a)

model7b <- lm(Age_onset_impairment ~ Age_at_Interview + Sex + ARRAY+ PC1+ PC2+ PC3+ PC4+ PC5, data = age_onset_impairment)
summary(model7b)

compare_performance(model7a, model7b)

coef(model7a)





