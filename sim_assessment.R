
# Simulation study of different methods for estimating errors

# GOAL: Show that spatial blocking matters, i.e. sampling trees rather
#       than blocks optimistically biases error predictions

# METHODS:
# Data: simulated data
#  - 500 observations on 100x100 grid
#  - Spatially correlated variables X1, X2, and X3
#  - Spatially correlated noise
# True model: function of X1, X2, X3
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
set.seed(300)
MAX_X <- 100   # Boundaries of data to use
MAX_Y <- 100
NUMSAMP <- 100   # Number of simulation samples to generate
NPOINTS <- 5000   # Number of points per sample

NFOLDS <- 10   # Number of folds for K-fold cv
NROW <- 4   # Number of rows for grid CV
NCOL <- 4   # Number of columns for grid CV

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
all_samps <- lapply(all_samps, function(x) {generate_sample(f, npoints=NPOINTS)})

# Compute true error and error estimates (1-5) for each sample
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
  errs_test <- unlist(lapply(all_samps, get_err))
  # Compute true error (need to remove error for df_train)
  err_true <- (sum(errs_test)-errs_test[i]) / (NUMSAMP-1)
  
  # Compute training error (1)
  err_train <- rmse(predict(model_lm, df_train), df_train$growth)
  
  # Compute K-fold CV error (2)
  err_cv <- cv_lm(df_train, NFOLDS, model_formula=model_formula)
  
  # Compute flipped K-fold CV error (3)
  err_cv_flip <- cv_lm_flip(df_train, NFOLDS, model_formula=model_formula)
  
  # Compute buffered grid CV error (4)
  err_grid_buffer <- cv_grid_buffer_lm(df_train, NROW, NCOL, model_formula=model_formula)
  
  # Compute flipped buffered grid CV error (5)
  err_grid_buffer_flip <- cv_grid_buffer_lm_flip(df_train, NROW, NCOL, model_formula=model_formula)
  
  # Output results
  err_df <- rbind(err_df, c(i, err_true, err_train, err_cv, err_cv_flip, err_grid_buffer, err_grid_buffer_flip))
}
colnames(err_df) <- c("set_number", "err_true", "err_train", 
                      "err_cv", "err_cv_flip", "err_grid_buffer", "err_grid_buffer_flip")

# Write output
output <- err_df
file_name <- paste("Documents/Thesis/Simulations/Results/final_results/sim-",
                   MAX_X, "x", MAX_Y, "-", NPOINTS,
                   "_cv", NCOL, "x", NROW, 
                   "-k", NFOLDS,
                   "_flipped",
                   ".csv", sep="")
write_csv(output, file_name)
toc()



#######
# EDA #
#######

# Plot learning curve for linear model
if (FALSE) {
  err_df <- data.frame(matrix(ncol = 2, nrow = 0))
  temp_train <- all_samps[[90]]
  temp_test <- all_samps[[43]]
  
  nrow <- nrow(temp_train)
  npoints <- 20   # number of learning curve points
  inc <- floor(nrow/npoints)   # increments between learning curve points
  for (i in 5:npoints) {
    train <- temp_train[1:(i*inc),]
    
    # Train lm
    model_formula <- "growth ~ X1 + X2 + X3"
    model_lm <- lm(model_formula, data=train)
    
    # Get training error (on whole training set)
    err_train <- sqrt(mean(model_lm$residuals^2))
    # Get test error
    err_test <- rmse(
      predict(model_lm, temp_test),
      temp_test$growth)
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
