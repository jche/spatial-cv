
# Bootstrap procedures

# INPUT:
# - data frame
# - 
# - (model formula)
# - (model parameters)
#
# OUTPUT:
# - bootstrap error estimate (rmse)

require(tidyverse)
require(rpart)
require(ModelMetrics)
require(sperrorest)
# require(boot)

#################################
# Standard bootstrap procedures #
#################################

# LOO Bootstrap for lm
boot_LOO_lm <- function(df, nboot=500, model_formula) {
  N <- nrow(df)
  boot_errors <- vector(mode="numeric", length=nboot)
  for (i in 1:nboot) {
    # Partition training/validation sets
    samp_indices <- sample(1:N, N, replace=TRUE)
    train <- df[unique(samp_indices),]
    test <- df[-samp_indices,]
    # Fit lm
    fitted_model <- lm(model_formula, data = train)
    # Return error
    diffs <- predict(fitted_model, test) - test$growth
    boot_errors[i] <- mean(diffs^2)
  }
  return(sqrt(mean(boot_errors)))
}

# Check stability of bootstrap
# err_LOO <- vector(mode = "numeric", length = 25)
# for (i in seq(20,500,20)) {
#   err_LOO[i/20] <- boot_LOO_lm(df_train, nboot=i, model_formula=model_formula)
# }
# temp <- data.frame(1:500, err_LOO) %>%
#   filter(!is.na(err_LOO)) %>%
#   filter(err_LOO!=0)
# ggplot(temp) +
#   geom_point(aes(x=X1.500,y=err_LOO))

################################
# Spatial bootstrap procedures #
################################

# Fixed tile bootstrap
boot_fixedtile_lm <- function(df, nrow, ncol, model_formula) {
  model <- sperrorest(
    model_formula,
    data=as.data.frame(df),
    coords=c("x", "y"),   # variables that contain x/y coordinates
    model_fun=glm,
    pred_fun=predict,
    smp_fun=represampling_tile_bootstrap,   # fixed tile bootstrap
    smp_args = list(repetition = 1,
                    nsplit = c(nrow, ncol),
                    reassign = FALSE,
                    seed1 = 44),   # 1 repetition
    progress=FALSE)
  return(summary(model$error_fold)["test.rmse",]$mean)
}

# Moving circle bootstrap
boot_movingcircle_lm <- function(df, radius, nboot, model_formula) {
  # Note: sperrorest does not do a LOO moving-tile bootstrap,
  # so I construct it manually
  parti <- represampling_disc_bootstrap(
    as.data.frame(df),
    radius=radius,
    nboot=nboot)
  samp_indices <- parti[[1]][1]$`1`$train
  
  # Partition training/validation sets
  train <- df[unique(samp_indices),]
  test <- df[-samp_indices,]
  # print(nrow(unique(test))); print(nrow(unique(train)))
  # Fit lm
  fitted_model <- lm(model_formula, data=train)
  # Return error
  diffs <- predict(fitted_model, test) - test$growth
  return(sqrt(mean(diffs^2)))
}

boot_movingtiletw_lm <- function(df, max_x, max_y, size, nboot, model_formula) {
  # Function to randomly sample tile of width `size`, w/ toroidal wrapping
  sample_tile_tw <- function(df, size=size) {
    # Randomly generate center of tile
    x_coord <- runif(1, 0, max_x)
    y_coord <- runif(1, 0, max_y)
    # Sample tile, w/ toroidal wrapping
    temp <- as.data.frame(df) %>%
      mutate(index = row_number()) %>%
      filter(
        x > x_coord-(size/2) | x < x_coord-max_x+(size/2),
        y > y_coord-(size/2) | y < y_coord-max_y+(size/2),
        x < x_coord+(size/2) | x > x_coord+max_x-(size/2),
        y < y_coord+(size/2) | y > y_coord+max_y-(size/2)
      )
    return(temp$index)
  }
  
  # Sample nboot tiles
  train_indices <- c()
  for (i in 1:nboot) {
    train_indices <- c(train_indices, sample_tile_tw(df, size))
  }
  train_indices <- unique(train_indices)
  
  # Partition training/validation sets
  train <- df[train_indices,]
  test <- df[-train_indices,]
  # Fit lm
  fitted_model <- lm(model_formula, data=train)
  # Return error
  diffs <- predict(fitted_model, test) - test$growth
  return(sqrt(mean(diffs^2)))
}

