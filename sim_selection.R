
# Simulation study of different methods for model selection

# GOAL: Show that model selection based on error estimates
#       from spatial blocking leads to more accurate model
#       selection (i.e. choosing the 'correct' model)

# METHODS:
# Data: simulated data
#  - 500 observations on 100x100 grid
#  - Spatially correlated variables X1...X6
#  - Spatially correlated noise
# True model: function of X1...X6
# Fit model: linear regression
# Blocking strategies:
#  1) Training error
#  2) K-fold CV
#  3) Flipped K-fold CV
#  4) LOO CV
#  5) Buffered Grid CV
#  6) Flipped Buffered Grid CV
#  7) Spatial LOO CV
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
SEED <- 630
set.seed(SEED)
MAX_X <- 100   # Boundaries of data to use
MAX_Y <- 100
NUMSAMP <- 100   # Number of simulation samples to generate
NPOINTS <- 500   # Number of points per sample

NFOLDS <- 10   # Number of folds for K-fold cv
NROW <- 4   # Number of rows for grid CV
NCOL <- 4   # Number of columns for grid CV
BUFFER <- 10   # Buffer for SLOO CV

for (tuple in list(c(F,F), c(F,T), c(T,F), c(T,T))) {
tic()
  
SP_VARS <- tuple[1]   # Whether vars are spat. correlated
SP_NOISE <- tuple[2]   # Whether noise is spat. correlated
# SP_VARS <- TRUE   # Whether vars are spat. correlated
# SP_NOISE <- TRUE   # Whether noise is spat. correlated

# Define true model
f <- function(a,b,c,d,e,f){
  return(a+b+c)
}

# Generate random samples
all_samps <- vector(mode = "list", length = NUMSAMP)
all_samps <- lapply(all_samps, function(x) {
  generate_sample(f, spat_vars=SP_VARS, spat_noise=SP_NOISE, npoints=NPOINTS)})

# Compute true error and error estimates (1-7) for each sample, model
# --> Select model with least error
sel_df <- data.frame()
errs_df <- data.frame()
for (i in 1:NUMSAMP) {
  # Get sample i
  df_train <- all_samps[[i]]
  model_formulas <- list(
    as.formula("val ~ X1 + X2"),
    as.formula("val ~ X1 + X3"),
    as.formula("val ~ X2 + X3"),
    as.formula("val ~ X1 + X2 + X3"),
    as.formula("val ~ X1 + X2 + X3 + X4"),
    as.formula("val ~ X1 + X2 + X3 + X5"),
    as.formula("val ~ X1 + X2 + X3 + X6"))
  
  # Compute training error (1) for each model
  get_err_train <- function(formula) {
    model_lm <- lm(formula, data=df_train)
    return(rmse(predict(model_lm, df_train), df_train$val))
  }
  errs_train <- lapply(model_formulas, get_err_train)
  # Select model with least training error
  sel_train <- which.min(errs_train)
  
  # Compute K-fold CV error (2) for each model
  get_err_cv <- function(formula) {
    return(cv_lm(df_train, NFOLDS, model_formula=formula))
  }
  errs_cv <- lapply(model_formulas, get_err_cv)
  # Select model with least K-fold CV error
  sel_cv <- which.min(errs_cv)
    
  # Compute flipped K-fold CV error (3) for each model
  get_err_cv_flip <- function(formula) {
    return(cv_lm_flip(df_train, NFOLDS, model_formula=formula))
  }
  errs_cv_flip <- lapply(model_formulas, get_err_cv_flip)
  # Select model with least K-fold CV error
  sel_cv_flip <- which.min(errs_cv_flip)
  
  # Compute LOO CV error (4) for each model
  get_err_loo <- function(formula) {
    return(cv_lm(df_train, NPOINTS, model_formula=formula))
  }
  errs_loo <- lapply(model_formulas, get_err_loo)
  # Select model with least K-fold CV error
  sel_loo <- which.min(errs_loo)
  
  # Compute buffered grid CV error (5) for each model
  get_err_grid <- function(formula) {
    return(cv_grid_buffer_lm(df_train, NROW, NCOL, model_formula=formula))
  }
  errs_grid <- lapply(model_formulas, get_err_grid)
  # Select model with least buffered grid CV error
  sel_grid <- which.min(errs_grid)
  
  # Compute flipped buffered grid CV error (6) for each model
  get_err_grid_flip <- function(formula) {
    return(cv_grid_buffer_lm_flip(df_train, NROW, NCOL, model_formula=formula))
  }
  errs_grid_flip <- lapply(model_formulas, get_err_grid_flip)
  # Select model with least buffered grid CV error
  sel_grid_flip <- which.min(errs_grid_flip)
  
  # Compute SLOO CV error (7) for each model
  get_err_sloo <- function(formula) {
    return(cv_SLOO_lm(df_train, BUFFER, model_formula=formula))
  }
  errs_sloo <- lapply(model_formulas, get_err_sloo)
  # Select model with least buffered grid CV error
  sel_sloo <- which.min(errs_sloo)

  # Output results of both selection and assessment
  sel_df <- rbind(sel_df, c(i, sel_train, sel_cv, sel_cv_flip, sel_loo, sel_grid, sel_grid_flip, sel_sloo))
}
colnames(sel_df) <- c("set_number", "sel_train", "sel_cv", "sel_cv_flip", "sel_loo", 
                      "sel_grid", "sel_grid_flip", "sel_sloo")

# Output results from model selection
output <- sel_df
file_name <- paste("Results/selection/sim-",
                   "spvars-", SP_VARS,
                   "_spnoise-", SP_NOISE,
                   "_", MAX_X, "x", MAX_Y, "-", NPOINTS,
                   "_cv", NCOL, "x", NROW,
                   "-k", NFOLDS,
                   "-buffer", BUFFER,
                   "_seed", SEED,
                   ".csv", sep="")
write_csv(output, file_name)

toc()
}

