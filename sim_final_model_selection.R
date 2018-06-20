
# Simulation study of different methods for model selection

# GOAL: Show that model selection based on error estimates
#       from spatial blocking leads to more accurate model
#       selection (i.e. choosing the 'correct' model)

# METHODS:
# Data: simulated data
#  - 500 observations on 100x100 grid
#  - Spatially correlated variables X1, X2, and X3
#  - Spatially correlated noise
# True model: function of X1, X2
# Fit model: linear regression
# Blocking strategies:
#  1) Training error
#  2) K-fold CV
#  3) LOO Bootstrap
#  4) Buffered Grid CV
#  5) Spatial LOO CV
#  6) Moving Circle Bootstrap
# Error metric: RMSE

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
set.seed(600)
MAX_X <- 100   # Boundaries of data to use
MAX_Y <- 100
NUMSAMP <- 100   # Number of simulation samples to generate
NPOINTS <- 500   # Number of points per sample

NFOLDS <- 10   # Number of folds for K-fold cv
NROW <- 4   # Number of rows for grid CV
NCOL <- 4   # Number of columns for grid CV
RADIUS_SLOO <- 15   # Size of buffer for SLOO
RADIUS <- 5   # Size of each moving circle bootstrap sample
NBOOT <- 30   # Number of moving circle bootstrap samples to take

# Define true model
f <- function(x, y, z){
  b1 <- 2
  b2 <- 2
  return(b1*sin(pi*x) + b2*(y))
}

# Generate random samples
all_samps <- vector(mode = "list", length = NUMSAMP)
all_samps <- lapply(all_samps, function(x) {generate_sample(f)})

# Compute true error and error estimates (1-5) for each sample, model
# --> Select model with least error
sel_df <- data.frame()
errs_df <- data.frame()
for (i in 1:NUMSAMP) {
  # Get sample i
  df_train <- all_samps[[i]]
  model_formulas <- list(
    as.formula("growth ~ X1"),
    as.formula("growth ~ X2"),
    as.formula("growth ~ X1 + X2"),
    as.formula("growth ~ X1 + X2 + X3"))
  
  # Compute training error (1) for each model
  get_err_train <- function(formula) {
    model_lm <- lm(formula, data=df_train)
    return(rmse(predict(model_lm, df_train), df_train$growth))
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
    
  # Compute LOO bootstrap error (3) for each model
  get_err_LOO <- function(formula) {
    return(boot_LOO_lm(df_train, model_formula=formula))
  }
  errs_LOO <- lapply(model_formulas, get_err_LOO)
  # Select model with least LOO bootstrap error
  sel_LOO <- which.min(errs_LOO)
  
  # Compute buffered grid CV error (4) for each model
  get_err_grid_buffer <- function(formula) {
    return(cv_grid_buffer_lm(df_train, NROW, NCOL, model_formula=formula))
  }
  errs_grid_buffer <- lapply(model_formulas, get_err_grid_buffer)
  # Select model with least buffered grid CV error
  sel_grid_buffer <- which.min(errs_grid_buffer)
  
  # Compute spatial LOO bootstrap error (5) for each model
  get_err_SLOO <- function(formula) {
    return(cv_SLOO_lm(df_train, RADIUS_SLOO, model_formula=formula))
  }
  errs_SLOO <- lapply(model_formulas, get_err_SLOO)
  # Select model with least LOO bootstrap error
  sel_SLOO <- which.min(errs_SLOO)
  
  # Compute moving circle bootstrap error (6) for each model
  get_err_moving <- function(formula) {
    return(boot_movingcircle_lm(df_train, RADIUS, NBOOT, model_formula=formula))
  }
  errs_moving <- lapply(model_formulas, get_err_moving)
  # Select model with least moving circle bootstrap error
  sel_moving <- which.min(errs_moving)

  # Output results of both selection and assessment
  sel_df <- rbind(sel_df, c(i, sel_train, sel_cv, sel_LOO, sel_grid_buffer, sel_SLOO, sel_moving))
  errs_df <- rbind(errs_df, data.frame(
    rep(i, 6),
    c("train", "cv", "LOO", "grid_buffer", "SLOO", "moving"),
    t(matrix(c(errs_train, errs_cv, errs_LOO, errs_grid_buffer, errs_SLOO, errs_moving), nrow=4, ncol=6))))
}
colnames(sel_df) <- c("set_number", "sel_train", "sel_cv", 
                      "sel_LOO", "sel_grid_buffer", "sel_SLOO", "sel_moving")
colnames(errs_df) <- c("set_number", "method", "err1", "err2", "err3", "err4")

# Output results from model selection
output <- sel_df
file_name <- paste("Documents/Thesis/Simulations/Results/final_results/selection-",
                   MAX_X, "x", MAX_Y, "-", NPOINTS,
                   "_cv", NCOL, "x", NROW, 
                   "-k", NFOLDS, "-sloo", RADIUS_SLOO,
                   "_boot", NBOOT, "x", RADIUS,
                   ".csv", sep="")
write_csv(output, file_name)

# Output results from model assessment
output <- errs_df %>%
  mutate(
    err1 = as.numeric(err1),
    err2 = as.numeric(err2),
    err3 = as.numeric(err3),
    err4 = as.numeric(err4))
file_name <- paste("Documents/Thesis/Simulations/Results/final_results/assessment-",
                   MAX_X, "x", MAX_Y, "-", NPOINTS,
                   "_cv", NCOL, "x", NROW, 
                   "-k", NFOLDS, "-sloo", RADIUS_SLOO,
                   "_boot", NBOOT, "x", RADIUS,
                   ".csv", sep="")
write_csv(output, file_name)
toc()

