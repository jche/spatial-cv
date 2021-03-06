
# Simulation study of different methods for estimating errors

# GOAL: Show that spatial blocking matters, i.e. sampling trees rather
#       than blocks optimistically biases error predictions

# METHODS:
# Data: simulated data
#  - 500 observations on 100x100 grid
#  - Spatially correlated variables X1...X10
#  - Spatially correlated noise
# True model: function of X1...X10
# Fit model: linear regression
# Blocking strategies:
#  1) Training error
#  2) K-fold CV
#  3) LOO CV
#  4) Buffered Grid CV
#  5) Spatial LOO CV
# Error metric: RMSE

# Setup
library(tidyverse)
library(tictoc)
source("gen_sample.R")
source("cv_procedures.R")

####################################
# Define parameters for simulation #
####################################
SEED <- 100
set.seed(SEED)
MAX_X <- 100   # Boundaries of data to use
MAX_Y <- 100
NUMSAMP <- 100   # Number of simulation samples to generate
NPOINTS <- 500   # Number of points per sample

NFOLDS <- 16   # Number of folds for K-fold cv
NROW <- 4   # Number of rows for grid CV
NCOL <- 4   # Number of columns for grid CV
BUFFER <- 15   # Buffer for SLOO CV

for (tuple in list(c(F,F), c(T,T))) {
tic()
  
SP_VARS <- tuple[1]   # Whether vars are spat. correlated
SP_NOISE <- tuple[2]   # Whether noise is spat. correlated
# SP_VARS <- TRUE   # Whether vars are spat. correlated
# SP_NOISE <- TRUE   # Whether noise is spat. correlated

# Define true model
f <- function(a,b,c){
  # original function
  return(2*sin(pi*a)+a+2*b+4*(c>0))
}
model_formula <- as.formula("val ~ X1 + X2 + X3")

# Generate random samples
all_samps <- vector(mode = "list", length = NUMSAMP)
all_samps <- lapply(all_samps, function(x) {
  generate_sample(f, spat_vars=SP_VARS, spat_noise=SP_NOISE, npoints=NPOINTS)})

# Compute true error and error estimates (1-5) for each sample
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
  
  # Compute K-fold CV error (2)
  err_cv <- cv_lm(df_train, NFOLDS, model_formula=model_formula)
  
  # Compute LOO CV error (3)
  err_loo <- cv_lm(df_train, NPOINTS, model_formula=model_formula)
  
  # Compute buffered grid CV error (4)
  err_grid <- cv_grid_buffer_lm(df_train, NROW, NCOL, model_formula=model_formula)
  
  # Compute SLOO CV error (5)
  err_sloo <- cv_SLOO_lm(df_train, BUFFER, model_formula=model_formula)
  
  # Output results
  err_df <- rbind(err_df, c(i, err_true, err_train, err_cv, err_loo, err_grid, err_sloo))
}
colnames(err_df) <- c("set_number", "err_true", "err_train", "err_cv",
                      "err_loo", "err_grid_buffer", "err_sloo")

# Write output
output <- err_df
file_name <- paste("Results/assessment/sim-",
                   "spvars-", SP_VARS,
                   "_spnoise-", SP_NOISE,
                   "_", MAX_X, "x", MAX_Y, "-", NPOINTS,
                   "_cv", NCOL, "x", NROW,
                   "-k", NFOLDS,
                   "-buffer", BUFFER,
                   "_seed", SEED,
                   # "_noisier",
                   "_plusA_noFlip_decay10",
                   ".csv", sep="")
write_csv(output, file_name)
toc()

}

#######
# EDA #
#######

# Plot learning curve for linear model
if (FALSE) {
  # Define custom f
  f <- function(a,b,c){
    return(2*sin(pi*a)+2*b+4*(c>0))
  }
  temp_train <- generate_sample(f, spat_vars=SP_VARS, spat_noise=SP_NOISE, npoints=NPOINTS)
  temp_test <- generate_sample(f, spat_vars=SP_VARS, spat_noise=SP_NOISE, npoints=NPOINTS)
  
  err_df <- data.frame(matrix(ncol = 2, nrow = 0))
  
  nrow <- nrow(temp_train)
  npoints <- 20   # number of learning curve points
  inc <- floor(nrow/npoints)   # increments between learning curve points
  for (i in 5:npoints) {
    train <- temp_train[1:(i*inc),]
    
    # Train lm
    model_formula <- "val ~ X1 + X2 + X3"
    model_lm <- lm(model_formula, data=train)
    
    # Get training error (on whole training set)
    err_train <- sqrt(mean(model_lm$residuals^2))
    # Get test error
    err_test <- rmse(
      predict(model_lm, temp_test),
      temp_test$val)
    # Store errors
    err_df <- rbind(err_df, c(err_train, err_test))
  }
  err_df <- cbind(err_df, seq(5*inc, nrow, by=inc))
  colnames(err_df) <- c("err_train", "err_test", "num")
  
  # err_df %>%
  #   ggplot(aes(x=num)) +
  #   geom_point(aes(y=err_train))
  err_df %>%
    ggplot(aes(x=num)) +
    geom_point(aes(y=err_test)) +
    geom_point(aes(y=err_train), color="red")
}

# Plot residuals in space for linear model
if (FALSE) {
  require(broom)
  
  temp_train <- all_samps[[1]]
  temp_lm <- lm(model_formula, data=temp_train)
  
  augment(temp_lm, temp_train) %>%
    ggplot(aes(x=x, y=y)) +
    geom_point(aes(color=.resid), size=2) +
    scale_color_continuous(low="blue", high="orange")
}

