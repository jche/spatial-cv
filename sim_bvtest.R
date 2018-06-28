
# Simulation study of different methods for estimating errors

# GOAL: Show bias-variance tradeoff in K for K-fold CV

# METHODS:
# Data: simulated data
#  - 500 observations on 100x100 grid
#  - Spatially correlated variables X1...X10
#  - Spatially correlated noise
# True model: function of X1...X10
# Fit model: linear regression
# Blocking strategies:
#  1) 2-fold CV
#  2) 10-fold CV
#  3) 50-fold CV
#  4) LOO CV
# Error metric: RMSE

# Setup
library(tidyverse)
library(tictoc)
source("gen_sample.R")
source("cv_procedures.R")
tic()

####################################
# Define parameters for simulation #
####################################
SEED <- 2
set.seed(SEED)
MAX_X <- 100   # Boundaries of data to use
MAX_Y <- 100
NUMSAMP <- 100   # Number of simulation samples to generate
NPOINTS <- 500   # Number of points per sample

NFOLDS <- c(2,3,4,5,6,7,8,9,10,25,50,100,NPOINTS)   # Number of folds for K-fold cv

SP_VARS <- FALSE   # Whether vars are spat. correlated
SP_NOISE <- FALSE   # Whether noise is spat. correlated

# Define true model
f <- function(a,b,c,d,e,f){
  # return(a+b+c+d+e+f)
  return(2*sin(pi*a)+2*b+4*(c>0)+2*sin(pi*d)+2*e+4*(f>0))
}
model_formula <- as.formula("val ~ X1 + X2 + X3 + X4 + X5 + X6")

# Generate random samples
all_samps <- vector(mode = "list", length = NUMSAMP)
all_samps <- lapply(all_samps, function(x) {
  generate_sample(f, spat_vars=SP_VARS, spat_noise=SP_NOISE, npoints=NPOINTS)})

# Compute true error and error estimates for each sample
err_df <- data.frame()
for (i in 1:NUMSAMP) {
  # Get sample i to use for training
  df_train <- all_samps[[i]]
  
  # Fit lm to df_train
  model_lm <- lm(model_formula, data=df_train)
  
  # Compute test errors on all other training sets
  get_err <- function(df) {
    rmse(predict(model_lm, df), df$val)
  }
  errs_test <- unlist(lapply(all_samps, get_err))
  # Compute true error (need to remove error for df_train)
  err_true <- (sum(errs_test)-errs_test[i]) / (NUMSAMP-1)
  
  # Compute training error (1)
  err_train <- rmse(predict(model_lm, df_train), df_train$val)
  
  # Compute K-fold CV errors (2)
  err_cv <- sapply(NFOLDS, function(x) {
    cv_lm(df_train, x, model_formula=model_formula)})

  # Output results
  err_df <- rbind(err_df, c(i, err_true, err_train, err_cv))
}
colnames(err_df) <- c("set_number", "err_true", "err_train", "err_cv2",
                      "err_cv3", "err_cv4", "err_cv5", "err_cv6",
                      "err_cv7", "err_cv8", "err_cv9", "err_cv10",
                      "err_cv25", "err_cv50", "err_cv100", "err_loo")

# Write output
output <- err_df
file_name <- paste("Results/bvtest/sim-",
                   "spvars-", SP_VARS,
                   "_spnoise-", SP_NOISE,
                   "_", MAX_X, "x", MAX_Y, "-", NPOINTS,
                   "-k", NFOLDS[1],
                   "thru", NFOLDS[13],
                   "_seed", SEED,
                   "_complexsignal",
                   ".csv", sep="")
write_csv(output, file_name)
toc()

