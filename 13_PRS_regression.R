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
# creating new co-variable

PRS_filt <- PRS %>%
  filter(Age_at_Interview !=999) # lost 34 individuals 


####### no episodes mania LE as a continuous variable ########
no_mania <- PRS_filt %>%
  filter(Number_of_episodes_mania_LE !=999) %>%
  filter(Age_onset_impairment !=999) %>%
  mutate(length_of_disorder = (Age_at_Interview - Age_onset_impairment)) # lost 497 individuals 

no_mania$Number_of_episodes_mania_LE <- factor(no_mania$Number_of_episodes_mania_LE, ordered = TRUE)

# as an ordinal variable 
m1a <- polr(Number_of_episodes_mania_LE ~ PRSadj + Age_at_Interview + Sex + ARRAY + PC1+ PC2+ PC3+ PC4+ PC5 + length_of_disorder, data = no_mania, Hess = TRUE)
#without PRS
m1b <- polr(Number_of_episodes_mania_LE ~ Age_at_Interview + Sex++ ARRAY + PC1+ PC2+ PC3+ PC4+ PC5 + length_of_disorder, data = no_mania, Hess = TRUE)

compare_performance(m1a, m1b)

#obtaining the p-value 
summary(m1a_ordinal) #tval = -1.6329
p_values <- 2 * (1 - pnorm(abs(-1.6329)))
print(p_values)

####### no episodes depression LE #######
no_depression <- PRS_filt %>%
  filter(Number_of_episodes_depression_LE !=999) %>%
  filter(Age_onset_impairment !=999) %>%
  mutate(length_of_disorder = (Age_at_Interview - Age_onset_impairment)) # lost 563 individuals

no_depression$Number_of_episodes_depression_LE <- factor(no_depression$Number_of_episodes_depression_LE, ordered = TRUE)

#with PRS
m2a<- polr(Number_of_episodes_depression_LE ~ PRSadj + Age_at_Interview + Sex + ARRAY + PC1 + PC2 + PC3 + PC4 + PC5 + length_of_disorder, data = no_depression, Hess = TRUE)

# without PRS 
m2b <- polr(Number_of_episodes_depression_LE ~ Age_at_Interview + Sex + ARRAY + PC1 + PC2 + PC3 + PC4 + PC5 + length_of_disorder, data = no_depression, Hess = TRUE)

compare_performance(m2a, m2b)

#obtaining beta and pval
summary(m2a)
p_values <- 2 * (1 - pnorm(abs(-3.4983)))
print(p_values)

####### BADDS psychosis as a count variable #######
BADDS_P <- PRS_filt %>%
  filter(BADDS_P !=999) # lost 319 individuals 

BADDS_P$BADDS_P <- factor(BADDS_P$BADDS_P, ordered = TRUE)

# with PRS
m3a <- polr(BADDS_P ~ PRSadj + Age_at_Interview + Sex + ARRAY + PC1 + PC2 + PC3 + PC4 + PC5, data = BADDS_P, Hess = TRUE)

# without PRS
m3b <- polr(BADDS_P ~ Age_at_Interview + Sex + ARRAY + PC1 + PC2 + PC3 + PC4 + PC5, data = BADDS_P, Hess = TRUE)

compare_performance(m3b, m3a)

#obtaining beta and pval
summary(m3a)
p_values <- 2 * (1 - pnorm(abs(4.255)))
print(p_values)

####### BADDS mood incongruent #######
BADDS_I <- PRS_filt %>%
  filter(BADDS_I !=999) # lost 1765 individuals 

BADDS_I$BADDS_I <- factor(BADDS_I$BADDS_I, ordered = TRUE)

# with PRS
m4a <- polr(BADDS_I ~ PRSadj + Age_at_Interview + Sex + ARRAY + PC1 + PC2 + PC3 + PC4 + PC5, data = BADDS_I, Hess = TRUE)

# without PRS
m4b <- polr(BADDS_I ~ Age_at_Interview + Sex + ARRAY + PC1 + PC2 + PC3 + PC4 + PC5, data = BADDS_I, Hess = TRUE)

compare_performance(m4a, m4b)

#obtaining beta and pval
summary(m4a)
p_values <- 2 * (1 - pnorm(abs(0.62130)))
print(p_values)

####### BADDS mania ####### 
BADDS_M <- PRS_filt %>%
  filter(BADDS_M != 999) # lost 258 individuals 

BADDS_M$BADDS_M <- factor(BADDS_M$BADDS_M, ordered = TRUE)

# with PRS
m5a <- polr(BADDS_M ~ PRSadj + Age_at_Interview + Sex + ARRAY + PC1 + PC2 + PC3 + PC4 + PC5, data = BADDS_M, Hess = TRUE)

#without PRS
m5b <- polr(BADDS_M ~ Age_at_Interview + Sex + ARRAY + PC1 + PC2 + PC3 + PC4 + PC5, data = BADDS_M, Hess = TRUE)

compare_performance(m5a, m5b)

#obtaining beta and pval
summary(m5a)
p_values <- 2 * (1 - pnorm(abs(7.6938)))
print(p_values)

####### BADDS depression ####### 
BADDS_D <- PRS_filt %>%
  filter(BADDS_D != 999) # lost 386 individuals

BADDS_D$BADDS_D <- factor(BADDS_D$BADDS_D, ordered = TRUE)

# with PRS
m6a <- polr(BADDS_D ~ PRSadj + Age_at_Interview + Sex + ARRAY + PC1 + PC2 + PC3 + PC4 + PC5, data = BADDS_D, Hess = TRUE)

#without PRS
m6b <- polr(BADDS_D ~ Age_at_Interview + Sex + ARRAY + PC1 + PC2 + PC3 + PC4 + PC5, data = BADDS_D, Hess = TRUE)

compare_performance(m6a, m6b)

#obtaining beta and pval
summary(m6a)
p_values <- 2 * (1 - pnorm(abs(-3.1682)))
print(p_values)

####### Age of onset of impairment ####### 
age_onset_impairment <- PRS_filt %>%
  filter(Age_onset_impairment != 999) # 185 lost 

m7a <- lm(Age_onset_impairment ~ PRSadj + Age_at_Interview + Sex + ARRAY+ PC1+ PC2+ PC3+ PC4+ PC5, data = age_onset_impairment)
summary(model7a)

m7b <- lm(Age_onset_impairment ~ Age_at_Interview + Sex + ARRAY+ PC1+ PC2+ PC3+ PC4+ PC5, data = age_onset_impairment)
summary(model7b)

compare_performance(m7a,m7b)

#obtaining beta and p-value
summary(m7a)

####### Rapid cycling  ####### 
rapid_cycling <- PRS_filt %>%
  filter(Rapid_cycling != 9) # 1255 individuals lost 

rapid_cycling$Rapid_cycling <- factor(rapid_cycling$Rapid_cycling, ordered = TRUE)

m8a <- polr(Rapid_cycling ~ PRSadj + Age_at_Interview + Sex + ARRAY+ PC1+ PC2+ PC3+ PC4+ PC5, data = rapid_cycling, Hess = T)
m8b <- polr(Rapid_cycling ~ Age_at_Interview + Sex + ARRAY+ PC1+ PC2+ PC3+ PC4+ PC5, data = rapid_cycling, Hess=T)

compare_performance(m8a,m8b)

# obtaining beta and p-value
summary(m8a)
p_values <- 2 * (1 - pnorm(abs(-3.93512)))
print(p_values)


####### Suicidal Ideation ####### 
suicidal_ideation <- PRS_filt %>%
  filter(Suicidal_ideation_LE != 9) # lost 239 individuals

suicidal_ideation$Suicidal_ideation_LE <- factor(suicidal_ideation$Suicidal_ideation_LE, ordered = TRUE)

m9a <- polr(Suicidal_ideation_LE ~ PRSadj + Age_at_Interview + Sex + ARRAY+ PC1+ PC2+ PC3+ PC4+ PC5, data = suicidal_ideation, Hess = T)
m9b <- polr(Suicidal_ideation_LE ~ Age_at_Interview + Sex + ARRAY+ PC1+ PC2+ PC3+ PC4+ PC5, data = suicidal_ideation, Hess=T)

compare_performance(m9a, m9b)

#obtaining beta and p-value
summary(m9a)
p_values <- 2 * (1 - pnorm(abs(-2.0635)))
print(p_values)

