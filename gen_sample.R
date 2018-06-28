
# Build function to randomly generate spatially correlated data
#
# INPUT: 
#  - true function f:(X1...X6) -> val
#  - whether variables & noise are independent or spatially correlated
#  - (maximum x-coordinate)
#  - (maximum y-coordinate)
#  - (number of points)
#  - (SD of noise)
# OUTPUT: df with `npoints` points, each with:
#  - x/y coordinates
#  - X1, ..., X6 values, where Xi are spatially correlated
#  - Y values from f(X1...X6) + e, where e is noise
#
# Two kinds of noise are added:
#  1) Independent normal noise N(0, sigma)
#  2) Spatially autocorrelated noise

require(mvnfast)

generate_sample <- function(f, spat_vars, spat_noise, max_x=100, max_y=100, npoints=500, sigma=0.5) {
  # Generate x/y coordinates, uniformly across space
  #  - TODO: generate x/y coordinates with spatial clustering
  #          - can use mixture of normals in space to do so
  x <- runif(npoints, 0, max_x)
  y <- runif(npoints, 0, max_y)
  df <- data.frame(x,y)
  
  if(spat_vars){
    # Spatially dependent covariance matrix
    mat <- as.matrix(proxy::dist(df))   # distance matrix
    sigma_mat <- 1 - mat/(mat+1)   # decaying function to simulate distance decay
    sigma_mat <- (sigma^2)*sigma_mat   # diagonals equal variance
  } else{
    # Independent covariance matrix
    sigma_mat <- diag(npoints)
  }
  
  vars <- rmvn(
    n=3,   # number of variables to generate (without noise)
    mu=rep(0, npoints),   # vector of variable means
    sigma=sigma_mat, 
    ncores=4)
  
  if(spat_noise){
    if(!spat_vars){
      # Generate sp. dep. cov. matrix if haven't done so yet
      mat <- as.matrix(proxy::dist(df))   # distance matrix
      sigma_mat <- 1 - mat/(mat+1)   # decaying function to simulate distance decay
      sigma_mat <- (sigma^2)*sigma_mat   # diagonals equal variance
    }
    # Spatially dependent noise 
    noise <- as.vector(rmvn(
      n=1,
      mu=rep(0, npoints),
      sigma=sigma_mat, 
      ncores=4))
  } else{
    # Spatially independent noise
    noise <- rnorm(npoints, 0, sigma)
  }
  
  # Add another layer of nonspatial noise
  extra_noise <- rnorm(npoints, 0, sigma)
  noise <- noise+extra_noise
  
  df <- cbind(df, t(vars), noise)
  colnames(df) <- c("x","y","X1","X2","X3","eps")
  df <- df %>%
    mutate(val = f(X1,X2,X3) + eps)
  return(df)
}


#######
# EDA #
#######

if (FALSE) {
  require(tidyverse)
  
  f <- function(a,b,c){
    return(2*sin(pi*a)+2*b+4*(c>0))
  }
  foo <- generate_sample(f, spat_vars=FALSE, spat_noise=FALSE, sigma=0.8, npoints=1000)
  
  # Plot spatial structure of variables
  foo %>%
    ggplot(aes(x=x,y=y)) +
    geom_point(aes(color=X1), size=2) +
    scale_color_continuous(low="blue", high="orange")
  # Plot spatial structure of noise
  foo %>%
    ggplot(aes(x=x,y=y)) +
    geom_point(aes(color=eps), size=2) +
    scale_color_continuous(low="blue", high="orange")
  # Plot spatial structure of response
  foo %>%
    ggplot(aes(x=x,y=y)) +
    geom_point(aes(color=val), size=2) +
    scale_color_continuous(low="blue", high="orange")
  
  # See how Xi relate to f
  foo %>%
    ggplot(aes(x=X1,y=val)) +
    geom_point()
}


