# perfrom regression analysis on PRS score and multiple variables

library(dplyr)
library(data.table)
library(MASS)
library(ggplot2)
library(performance)

#loading data
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