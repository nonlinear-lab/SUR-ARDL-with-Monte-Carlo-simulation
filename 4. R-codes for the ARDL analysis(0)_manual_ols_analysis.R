#+++++++++++++++++++++++++++++++++
# Step 0: set a working directory
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
# Step 2: assign data
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

  #++++++++++++++++++++++++++++++++++++++++++++++++
  # Step 3: data analysis
  #++++++++++++++++++++++++++++++++++++++++++++++++
  y_lag  <- y[-T_obs]
  dy     <- diff(y)
  x_lag  <- x[-T_obs]
  dx     <- diff(x)
  y_t    <- y[-1]
  data_set <- data.frame(dy = dy, dx = dx, x_lag = x_lag, y_lag = y_lag)
  # restricted model
  model_r  <- lm(dy ~ dx, data = data_set)
  ssr_r <- sum(residuals(model_r)^2)
  # unrestricted model
  model_ur <- lm(dy ~ dx + y_lag + x_lag, data = data_set)
  ssr_ur <- sum(residuals(model_ur)^2)
  # bounds test
  df_num <- k + 1
  # Denominator df is Sample Size minus Number of Parameters in UR model
  df_den <- nrow(data_set) - length(coef(model_ur))
  bounds_stat <- ((ssr_r - ssr_ur) / df_num) / (ssr_ur / df_den)
  mout[j,1]=bounds_stat
}
print(" bounds test (k=1, Case III):")
print(mout)
