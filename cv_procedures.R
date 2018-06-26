
# Cross-validation procedures
# K-fold CV code from https://stats.stackexchange.com/questions/61090/how-to-split-a-data-set-to-do-10-fold-cross-validation

# INPUT:
# - data frame
# - # folds
# - (model formula)
# - (model parameters)
#
# OUTPUT:
# - CV error estimate (rmse)

require(tidyverse)
require(rpart)
require(ModelMetrics)
require(sperrorest)

########################
# K-fold CV procedures #
########################

# K-fold CV for tree (using rpart)
cv_tree <- function(df, nfolds, cp, model_formula=as.formula("val ~ X1")) {
  cv_errors <- vector(mode="numeric", length=nfolds)

  # Shuffle data and partition folds
  temp <- df[sample(nrow(df)),]
  folds <- cut(seq(1,nrow(temp)),
               breaks=nfolds,
               labels=FALSE)
  for(i in 1:nfolds){
    # Split training/testing sets
    test_index <- which(folds==i, arr.ind=TRUE)
    test <- temp[test_index, ]
    train <- temp[-test_index, ]
    # Fit tree
    tree_parameters <- rpart.control(cp=cp)
    fitted_model <- rpart(model_formula, 
                          data = train, 
                          control=tree_parameters)
    # Store error
    diffs <- predict(fitted_model, test) - test$val
    cv_errors[i] <- mean(diffs^2)
  }
  return(sqrt(mean(cv_errors)))
}

# K-fold CV for lm
cv_lm <- function(df, nfolds, model_formula) {
  cv_errors <- vector(mode="numeric", length=nfolds)
  
  # Shuffle data and partition folds
  temp <- df[sample(nrow(df)),]
  folds <- cut(seq(1,nrow(temp)),
               breaks=nfolds,
               labels=FALSE)
  for(i in 1:nfolds){
    # Split training/testing sets
    test_index <- which(folds==i, arr.ind=TRUE)
    test <- temp[test_index, ]
    train <- temp[-test_index, ]
    # Fit lm
    fitted_model <- lm(model_formula, data = train)
    # Store error
    diffs <- predict(fitted_model, test) - test$val
    cv_errors[i] <- mean(diffs^2)
  }
  return(sqrt(mean(cv_errors)))
}

# K-fold CV for spline
cv_spline <- function(df, nfolds, deg_free) {
  cv_errors <- vector(mode="numeric", length=nfolds)
  
  # Shuffle data and partition folds
  temp <- df[sample(nrow(df)),]
  folds <- cut(seq(1,nrow(temp)),
               breaks=nfolds,
               labels=FALSE)
  for(i in 1:nfolds){
    # Split training/testing sets
    test_index <- which(folds==i, arr.ind=TRUE)
    test <- temp[test_index, ]
    train <- temp[-test_index, ]
    # Fit spline
    model_spline <- smooth.spline(x=train$X1,
                                  y=train$val, 
                                  df=deg_free)
    # Store error
    diffs <- predict(model_spline, x=test$X1)$y - test$val
    cv_errors[i] <- mean(diffs^2)
  }
  return(sqrt(mean(cv_errors)))
}

# Flipped K-fold CV for lm
cv_lm_flip <- function(df, nfolds, model_formula) {
  # Empty dataframe with 1 column per observation in df
  obs_errors <- matrix(nrow=0, ncol=nrow(df))
  
  # Shuffle data and partition folds
  temp <- df[sample(nrow(df)),]
  folds <- cut(seq(1,nrow(temp)),
               breaks=nfolds,
               labels=FALSE)
  for(i in 1:nfolds){
    # Split training/testing sets
    train_index <- which(folds==i, arr.ind=TRUE)
    train <- temp[train_index, ]
    test <- temp[-train_index, ]
    # Fit lm
    fitted_model <- lm(model_formula, data = train)
    # Store error for each observation (0 for obs in train)
    err_vec <- vector(mode="numeric", length=nrow(df))
    diffs <- predict(fitted_model, test) - test$val
    err_vec[-train_index] <- diffs
    obs_errors <- rbind(obs_errors, err_vec)
  }
  # Average error for each observation, then rmse
  errs <- (rowSums(t(obs_errors))/(nfolds-1))^2
  return(sqrt(mean(errs)))
}

#########################
# Spatial CV procedures #
#########################

# Buffered grid CV
cv_grid_buffer_lm <- function(df, nrow, ncol, model_formula) {
  # Use sperrorest to partition data into spatial folds
  nfolds_spat <- nrow*ncol
  partition <- partition_tiles(
    as.data.frame(df),   # See testing_sperrorest_package.Rmd
    nsplit = c(ncol,nrow),   # nrow x ncol grid
    reassign = FALSE)   # Don't reassign observations in small folds
  
  # Hacky method to number which fold each observation is in
  all_folds <- partition$`1`
  count = 1
  fold_df <- df %>%
    select(x, y)
  for (fold in all_folds){
    varname <- paste("fold", count, sep="")
    fold_df <- fold_df %>%
      mutate(!!varname := count*(row_number() %in% fold$test))
    count <- count + 1
  }
  folds <- rowSums(fold_df[,c(3:ncol(fold_df))])
  df_spatial <- df %>%
    mutate(fold = folds)
  
  # Find fold neighbors of each fold (again, hacky)
  fold_neighbors <- vector(mode="list", length=nfolds_spat)
  for (i in 1:nfolds_spat) {
    list1 <- c(i, i-nrow, i+nrow)
    list2 <- switch(as.character(i %% nrow),
                    "0" = c(i-1, i-1-nrow, i-1+nrow),
                    "1" = c(i+1, i+1-nrow, i+1+nrow),
                    c(i-1, i-1-nrow, i-1+nrow, i+1, i+1-nrow, i+1+nrow))
    fold_neighbors[[i]] <- c(list1, list2)
  }
  
  # Compute buffered grid CV error
  results_list <- vector(mode="numeric", length=nfolds_spat)
  num_empty_folds <- 0
  for (j in 1:nfolds_spat) {
    train <- df_spatial %>% filter(!(fold %in% fold_neighbors[[j]]))
    test <- df_spatial %>% filter(fold == j)
    if (dim(test)[1] == 0) {   # If test set is empty
      num_empty_folds <- num_empty_folds + 1
      next
    }
    # Train model & make predictions
    model_lm <- lm(model_formula, data=train)
    model_preds <- predict(model_lm, test)
    # Store results
    results_list[[j]] <- rmse(model_preds, test$val)
  }
  return(sum(unlist(results_list)) / (nfolds_spat-num_empty_folds))
}

# Flipped buffered grid CV for lm
cv_grid_buffer_lm_flip <- function(df, nrow, ncol, model_formula) {
  # Use sperrorest to partition data into spatial folds
  nfolds_spat <- nrow*ncol
  partition <- partition_tiles(
    as.data.frame(df),   # See testing_sperrorest_package.Rmd
    nsplit = c(ncol,nrow),   # nrow x ncol grid
    reassign = FALSE)   # Don't reassign observations in small folds
  
  # Hacky method to number which fold each observation is in
  all_folds <- partition$`1`
  count = 1
  fold_df <- df %>%
    select(x, y)
  for (fold in all_folds){
    varname <- paste("fold", count, sep="")
    fold_df <- fold_df %>%
      mutate(!!varname := count*(row_number() %in% fold$test))
    count <- count + 1
  }
  folds <- rowSums(fold_df[,c(3:ncol(fold_df))])
  
  # Find fold neighbors of each fold (again, hacky)
  fold_neighbors <- vector(mode="list", length=nfolds_spat)
  for (i in 1:nfolds_spat) {
    list1 <- c(i, i-nrow, i+nrow)
    list2 <- switch(as.character(i %% nrow),
                    "0" = c(i-1, i-1-nrow, i-1+nrow),
                    "1" = c(i+1, i+1-nrow, i+1+nrow),
                    c(i-1, i-1-nrow, i-1+nrow, i+1, i+1-nrow, i+1+nrow))
    fold_neighbors[[i]] <- c(list1, list2)
  }
  
  # Compute flipped buffered grid CV error
  # TODO: include check for empty training sets
  obs_errors <- matrix(nrow=0, ncol=nrow(df))
  for (j in 1:nfolds_spat) {
    train_index <- which(folds == j, arr.ind=TRUE)
    test_index <- which(!(folds %in% fold_neighbors[[j]]), arr.ind=TRUE)
    train <- df[train_index,]
    test <- df[test_index,]
    # Train model & make predictions
    model_lm <- lm(model_formula, data=train)
    # Store error for each observation (0 for obs in train)
    err_vec <- vector(mode="numeric", length=nrow(df))
    diffs <- predict(model_lm, test) - test$val
    err_vec[test_index] <- diffs
    obs_errors <- rbind(obs_errors, err_vec)
  }
  # Average errors for each observation, accounting for fact that
  # each observation may have a different number of estimated errors
  errs <- vector(mode="numeric", length=nrow(df))
  for (k in 1:nrow(df)) {
    temp <- obs_errors[,k]
    num_zero <- sum(temp == 0)
    errs[k] <- (sum(temp)/(nfolds_spat-num_zero))^2
  }
  return(sqrt(mean(errs)))
}

# Unbuffered grid CV
cv_grid_nobuff_lm <- function(df, nrow, ncol, model_formula) {
  model <- sperrorest(
    model_formula,
    data=df,
    coords=c("x", "y"),   # variables that contain x/y coordinates
    model_fun=glm,
    pred_fun=predict,
    smp_fun=partition_tiles,   # spatial cv by rectangular grids
    smp_args = list(repetition = 1,
                    nsplit = c(nrow, ncol),
                    reassign = FALSE,
                    seed1 = 123),   # 1 repetition
    progress=FALSE)
  return(summary(model$error_fold)["test.rmse",]$mean)
}

# SLOO
cv_SLOO_lm <- function(df, buffer, model_formula) {
  model <- sperrorest(
    model_formula,
    data=df,
    coords=c("x", "y"),   # variables that contain x/y coordinates
    model_fun=glm,
    pred_fun=predict,
    smp_fun=partition_loo,   # SLOO
    smp_args = list(repetition = 1,
                    buffer = buffer,
                    seed1 = 123),   # 1 repetition
    progress=FALSE)
  return(summary(model$error_rep)["test_rmse",]$mean)
}
