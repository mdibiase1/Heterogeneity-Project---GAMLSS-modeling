
# GAMLSS Z-score Calculation Script for Biochemical Data (SHASH Distribution)
# Author: Maria Di Biase
# Description: This script fits a GAMLSS model using the SHASH distribution
# and computes z-scores for all subjects in the dataset based on predicted
# distribution parameters (mu, sigma, nu, tau).

rm(list = ls())

# Load necessary libraries
library(gamlss)
library(R.matlab)
library(pracma)

# Set working directory
setwd("GAMLSS_input_biochem/")

# Define output path
PATH_OUT <- "GAMLSS_output/"

# Define the biochemical measure
biochem_measure <- 'Eosinophill_percentage'
FN <- paste("GAMLSSinput_Biochem_", biochem_measure, ".mat", sep="")

# Load input data
M <- readMat(FN)
biochem_measure_data <- as.data.frame(M$biochem.measure)
mydata_tmp <- as.data.frame(cbind(M$AGE, M$SEX, M$data))
names(mydata_tmp) <- c("age", "sex", "phenotype")

# Filter healthy controls
ind <- which(M$DX[,1] %in% c(1))
mydata <- mydata_tmp[ind, ]

# Fit the GAMLSS model
fam_dist <- 'SHASH'
mdl <- gamlss(phenotype ~ fp(age, npoly=1) + sex,
              sigma.fo = ~fp(age, npoly=1),
              family = fam_dist,
              data = mydata,
              robust = TRUE)

# Generate z-scores for all subjects
n <- nrow(mydata_tmp)
z_scores <- rep(NA, n)
p <- rep(NA, n)

for (sub in 1:n) {
  sub_data <- data.frame(age = mydata_tmp$age[sub],
                         sex = mydata_tmp$sex[sub],
                         phenotype = mydata_tmp$phenotype[sub])

  # Predict parameters
  params <- predictAll(mdl, newdata = sub_data)

  # Calculate CDF value for observed phenotype
  p[sub] <- pSHASH(mydata_tmp$phenotype[sub],
                   mu = params$mu,
                   sigma = params$sigma,
                   nu = params$nu,
                   tau = params$tau)

  # Convert CDF to z-score using standard normal inverse
  z_scores[sub] <- qNO(p[sub])
}

# Export results
output_file <- paste(PATH_OUT, "GAMLSSout_main_ZSCORES", biochem_measure, ".mat", sep = "")
writeMat(output_file, z_scores_unseen_subs = z_scores)
