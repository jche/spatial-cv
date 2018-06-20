
# Build functions to randomly generate noisy growth values
# and spatially correlated variables
# INPUT: 
#  - (number of points)
#  - true function f:(X1,X2,X3) -> growth
# OUTPUT: df with npoints points, each with:
#  - x/y coordinates
#  - X1, X2, X3 values, where Xi are spatially correlated
#  - Y values from f(X1,X2,X3) + e, where e is spatially correlated
#
# Two kinds of noise are added:
#  1) Independent normal noise N(0, .15)
#  2) Spatially autocorrelated noise, using mvnfast package

require(mvnfast)

generate_sample <- function(f, max_x=100, max_y=100, npoints=500, sigma=0.5) {
  # Generate x/y coordinates, uniformly across space
  #  - TODO: generate x/y coordinates with spatial clustering
  x <- runif(npoints, 0, max_x)
  y <- runif(npoints, 0, max_y)
  df <- data.frame(x,y)
  
  # Generate spatially correlated X1, X2, X3, noise
  mat <- as.matrix(proxy::dist(df))   # distance matrix
  sigma_mat <- 1 - mat/(mat+10)   # decaying function to simulate distance decay
  sigma_mat <- (sigma^2)*sigma_mat   # diagonals equal variance
  # Generate nonspatial normal noise
  noise <- rnorm(npoints, 0, sigma)
  
  vars <- rmvn(
    n=4,   # number of variables to generate
    mu=rep(0, npoints),   # vector of variable means
    sigma=sigma_mat, 
    ncores=4)
  df <- cbind(df, vars[1,], vars[2,], vars[3,], vars[4,])
  colnames(df) <- c("x","y","X1","X2","X3","eps")
  df <- df %>%
    mutate(growth = f(X1,X2,X3) + eps + noise)
  return(df)
}


#######
# EDA #
#######

if (FALSE) {
  f <- function(x, y, z){
    b1 <- 2
    b2 <- 2
    b3 <- 4
    return(b1*sin(pi*x) + b2*(y) + b3*(z>0))
  }
  foo <- generate_sample(f, sigma=0.8, npoints=1000)
  
  # Plot spatial structure of variables
  foo %>%
    ggplot(aes(x=x,y=y)) +
    geom_point(aes(color=X3), size=2) +
    scale_color_continuous(low="blue", high="orange")
  
  # See how Xi relate to f
  foo %>%
    ggplot(aes(x=X1,y=growth)) +
    geom_point()
}


