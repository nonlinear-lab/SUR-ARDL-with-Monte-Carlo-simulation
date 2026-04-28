#+++++++++++++++++++++++++++++++++
# Step 0: set a dorking directory
#+++++++++++++++++++++++++++++++
setwd("C:/R/ardl")

#++++++++++++++++++++++++++++++++++++++++++++++++
# Step 1: load data
#++++++++++++++++++++++++++++++++++++++++++++++++
my_data <- read.csv("data(2).csv")
my_data <- as.matrix(my_data)
T_obs   <- nrow(my_data)
k       <- 1
mout    <- matrix(NA,nrow=5,ncol=1)

#++++++++++++++++++++++++++++++++++++++++++++++++
# Step 2: Bounds test function
#++++++++++++++++++++++++++++++++++++++++++++++++
ardl_bounds <- function(y,x){
  T_obs   <- nrow(y)
  k       <- 1    # k is number of regressor
  y_lag  <- y[-T_obs]
  dy     <- diff(y)
  x_lag  <- x[-T_obs]
  dx     <- diff(x)
  # restricted model
  X_r <- cbind(1, dx)
  beta_r <- solve(t(X_r) %*% X_r) %*% t(X_r) %*% dy
  res_r <- dy - (X_r %*% beta_r)
  SSR_R <- sum(res_r^2)
  # unrestricted model
  X_ur <- cbind(1, dx, y_lag, x_lag)
  beta_ur <- solve(t(X_ur) %*% X_ur) %*% t(X_ur) %*% dy
  res_ur <- dy - (X_ur %*% beta_ur)
  SSR_ur <- sum(res_ur^2)
  # degree of freedom
  df_num <- k + 1
  df_den <- length(dy) - ncol(X_ur)

  bounds_stat <- ((SSR_R - SSR_ur) / df_num) / (SSR_ur / df_den)
  return(bounds_stat)
}
#++++++++++++++++++++++++++++++++++++++++++++++++
# Step 3: Run the analysis
#++++++++++++++++++++++++++++++++++++++++++++++++

for (j in 1:5) {
  if (j == 1) {
    y <- as.matrix(my_data[, 3])
    x <- as.matrix(my_data[, 2])
  } else if (j == 2) {
    y <- as.matrix(my_data[, 4])
    x <- as.matrix(my_data[, 2])
  } else if (j == 3) {
    y <- as.matrix(my_data[, 5])
    x <- as.matrix(my_data[, 2])
  } else if (j == 4) {
    y <- as.matrix(my_data[, 6])
    x <- as.matrix(my_data[, 2])
  } else {
    y <- as.matrix(my_data[, 7])
    x <- as.matrix(my_data[, 2])
  }
  results <- ardl_bounds (y,x)
  mout[j,1]=results
}

print(" bounds test (k=1, Case III):")
print(mout)
