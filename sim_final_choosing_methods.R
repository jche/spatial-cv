
# Simple script to pick optimal values for:
#
# K: K-fold CV
# 
# NROW/NCOL: grid CV
# NROW/NCOL: buffered grid CV
# RADIUS: SLOO
# 
# NROW/NCOL: fixed-tile bootstrap
# SIZE: moving tile bootstrap (no toroidal)
# SIZE: moving tile bootstrap (with toroidal)
# RADIUS/NBOOT: moving circle bootstrap
#
# and also to choose between the methods (once they've been optimized:
#
# CV: grid, buffered grid, SLOO (too computationally intensive?), SKCV
# Bootstrap: fixed-tile, moving-tile (no toroidal), moving-tile (with toroidal), moving-circle

# GOAL: find optimal balance of simplicity and good performance


# Setup
library(tidyverse)
library(tictoc)
source("Documents/Thesis/Simulations/sim_final_generate_sample.R")
source("Documents/Thesis/Simulations/cv_procedures.R")
source("Documents/Thesis/Simulations/bootstrap_procedures.R")
tic()

####################################
# Define parameters for simulation #
####################################
set.seed(20)
NUMSAMP <- 100   # Number of simulation samples to generate
NPOINTS <- 500   # Number of points per sample
# Boundaries of data to use
MAX_X <- 100
MAX_Y <- 100

######################################### EDIT THIS  ######################################### 
param_vec <- c(1,5,10,15,20,25)

# Define true model
f <- function(x, y, z){
  b1 <- 2
  b2 <- 2
  b3 <- 4
  return(b1*sin(pi*x) + b2*(y) + b3*(z>0))
}
model_formula <- as.formula("growth ~ X1 + X2 + X3")

# Generate random samples
all_samps <- vector(mode = "list", length = NUMSAMP)
all_samps <- lapply(all_samps, function(x) {generate_sample(f)})

# Compute true error and error estimates for each sample
err_df <- data.frame()
for (i in 1:NUMSAMP) {
  # Get sample i
  df_train <- all_samps[[i]]
  
  # Fit lm to df_train
  model_lm <- lm(model_formula, data=df_train)
  # Compute test errors on all other training sets
  get_err <- function(df) {
    rmse(predict(model_lm, df), df$growth)
  }
  errs_test <- sapply(all_samps, get_err)
  # Compute true error (need to remove error for df_train)
  err_true <- (sum(errs_test)-errs_test[i]) / (NUMSAMP-1)
  
  ######################################### EDIT THIS  ######################################### 
  err_vec <- lapply(param_vec, function(x){
    cv_SLOO_lm(df_train, x, model_formula=model_formula)
  })

  # Output results
  err_df <- rbind(err_df, c(i, err_true, unlist(err_vec)))
}
colnames(err_df) <- c("set_number", "err_true", 
                      sapply(param_vec, function(x) {
                        paste("err", x, sep="")
                      }))

# Write output
output <- err_df
######################################### EDIT THIS  ######################################### 
file_name <- paste("Documents/Thesis/Simulations/Results/final_results/cv-SLOO_larger",
                   ".csv", sep="")
write_csv(output, file_name)
toc()

