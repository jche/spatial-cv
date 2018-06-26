
# Build function to randomly generate spatially correlated data
#
# INPUT: 
#  - true function f:(X1...X6) -> val
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

generate_sample <- function(f, max_x=100, max_y=100, npoints=500, sigma=0.5, independent=FALSE) {
  # Generate x/y coordinates, uniformly across space
  #  - TODO: generate x/y coordinates with spatial clustering
  #          - can use mixture of normals in space to do so
  x <- runif(npoints, 0, max_x)
  y <- runif(npoints, 0, max_y)
  df <- data.frame(x,y)
  
  if(!independent){
    # Generate spatially correlated X variables, noise
    mat <- as.matrix(proxy::dist(df))   # distance matrix
    sigma_mat <- 1 - mat/(mat+1)   # decaying function to simulate distance decay
    sigma_mat <- (sigma^2)*sigma_mat   # diagonals equal variance
  } else{
    # Generate independent X variables
    sigma_mat <- diag(npoints)
  }
  # Generate nonspatial normal noise
  noise <- rnorm(npoints, 0, sigma)
  
  vars <- rmvn(
    n=7,   # number of variables to generate (including noise)
    mu=rep(0, npoints),   # vector of variable means
    sigma=sigma_mat, 
    ncores=4)
  df <- cbind(df, t(vars))
  colnames(df) <- c("x","y","X1","X2","X3","X4","X5","X6","eps")
  df <- df %>%
    mutate(val = f(X1,X2,X3,X4,X5,X6) + eps + noise)
  return(df)
}


#######
# EDA #
#######

if (FALSE) {
  require(tidyverse)
  
  f <- function(a,b,c,d,e,f){
    return(a+b+c)
  }
  foo <- generate_sample(f, sigma=0.8, npoints=1000)
  
  # Plot spatial structure of variables
  foo %>%
    ggplot(aes(x=x,y=y)) +
    geom_point(aes(color=val), size=2) +
    scale_color_continuous(low="blue", high="orange")
  
  # See how Xi relate to f
  foo %>%
    ggplot(aes(x=X5,y=val)) +
    geom_point()
}


