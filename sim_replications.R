
# Replications of simulations in:
#  - Shao 1993: Linear Model Selection by Cross-Validation
#    (http://www.libpls.net/publication/MCCV_Shao_1993.pdf)
#  - Zhang 1993: Model Selection via Multifold Cross-Validation
#    (https://projecteuclid.org/download/pdf_1/euclid.aos/1176349027)

require(tidyverse)
source("cv_procedures.R")

# Shao
set.seed(1)

# Data from paper
df <- read_csv("shao_data.csv")

# Data-generating function from paper
f <- function(a,b,c,d,e){
  b1 <- 2
  b2 <- 0
  b3 <- 0
  b4 <- 4
  b5 <- 0
  return(b1*a + b2*b + b3*c + b4*d + b5*e)
}
# Models to consider from paper
model_formulas <- list(
  as.formula("val ~ X1 + X4"),
  as.formula("val ~ X1 + X2 + X4"),
  as.formula("val ~ X1 + X3 + X4"),
  as.formula("val ~ X1 + X4 + X5"),
  as.formula("val ~ X1 + X2 + X3 + X4"),
  as.formula("val ~ X1 + X2 + X4 + X5"),
  as.formula("val ~ X1 + X3 + X4 + X5"),
  as.formula("val ~ X1 + X2 + X3 + X4 + X5"))

results_vec <- vector(1000, mode="numeric")
for (i in 1:1000) {
  # Generate noisy df
  noise <- rnorm(40, 0, 1)
  df <- df %>%
    mutate(val = f(X1,X2,X3,X4,X5) + noise)

  # Compute LOOCV error for each model
  get_err_cv <- function(formula) {
    return(cv_lm(df, 40, model_formula=formula))
  }
  errs_cv <- lapply(model_formulas, get_err_cv)
  # Select model with least K-fold CV error
  sel_cv <- which.min(errs_cv)
  
  # Store results
  results_vec[i] <- sel_cv
}

table(results_vec)/1000
# 1     2     3     4     5     6     7     8 
# 0.489 0.147 0.108 0.134 0.050 0.030 0.029 0.013

# Same results as paper!


# Zhang
set.seed(1)

# Data-generating function from paper
f <- function(a, b){
  b1 <- .6
  return(b1*a)
}
# Models to consider from paper
model_formulas <- list(
  as.formula("val ~ 1"),
  as.formula("val ~ X1"),
  as.formula("val ~ X1 + X2"))
# K for K-fold CV from paper
k_vals <- c(2, 5, 10, 20)

results_mat <- matrix(nrow = 4, ncol = 500)
for (i in 1:500) {
  # Generate random df
  df <- data.frame(rnorm(20, 0, 1), rnorm(20, 0, 1))
  names(df) <- c("X1", "X2")
  df <- df %>%
    mutate(val = f(X1, X2))
  
  # For each K, for each model, choose least-error model
  err_mat <- sapply(k_vals, function(k_val) {
    lapply(model_formulas, function(formula) {
      return(cv_lm(df, k_val, model_formula=formula))
    })})
  sel_vec <- apply(err_mat, 2, which.min)
  
  # Store results
  results_mat[,i] <- sel_vec
}

apply(results_mat, 1, table)/500
#    [,1] [,2]  [,3]  [,4]
# 2 0.554  0.5 0.576 0.556
# 3 0.446  0.5 0.424 0.444

# Different results than paper...

