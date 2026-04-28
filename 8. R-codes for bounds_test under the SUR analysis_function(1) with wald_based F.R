#+++++++++++++++++++++++++++++++++
# Step 0: set a working directory
#+++++++++++++++++++++++++++++++
setwd("C:/R/ardl")
library(systemfit)

#+++++++++++++++++++++++++++++++++
# Step 1: Load Data
#+++++++++++++++++++++++++++++++
data_set <- read.csv("data(5).csv")
data_set <- data_set[, -1]
colnames(data_set) <- c("x1", "y1", "x2", "y2", "x3", "y3", "x4", "y4", "x5", "y5")
T_obs   <- nrow(data_set)
k       <- 1
mout    <- matrix(NA,nrow=5,ncol=1)

#++++++++++++++++++++++++++++++++++++++++++++++
# Step 2: Wald test function
#++++++++++++++++++++++++++++++++++++++++++++++
wald_sur_eq <- function(model, restriction) {
  b <- coef(model)
  V <- vcov(model)
  idx <- match(restriction, names(b))
  if (any(is.na(idx))) {
    cat("\nAvailable coefficient names are:\n")
    print(names(b))
    cat("\nRequested restrictions were:\n")
    print(restriction)
    stop("Some restriction names were not found.")
  }
  R <- matrix(0, nrow = length(idx), ncol = length(b))
  for (i in seq_along(idx)) {
    R[i, idx[i]] <- 1
  }
  Rb <- R %*% b
  RVRT <- R %*% V %*% t(R)
  # Pure Wald statistic
  W <- as.numeric(t(Rb) %*% solve(RVRT) %*% Rb)
  F_version <- W / length(idx)
  return(F_version)
  return(W)
}

#++++++++++++++++++++++++++++++++++++++++++++++
# Step 3: SUR_ARDL bounds test function
#++++++++++++++++++++++++++++++++++++++++++++++
sur_ardl_bounds <- function(x1,x2,x3,x4,x5,y1,y2,y3,y4,y5){
  # creating difference and lag value
  for (j in 1:5) {
    temp_x <- get(paste0("x", j))
    temp_y <- get(paste0("y", j))
    assign(paste0("x_lag", j), temp_x[-T_obs])
    assign(paste0("dx",    j), diff(temp_x))
    assign(paste0("y_lag", j), temp_y[-T_obs])
    assign(paste0("dy",    j), diff(temp_y))
  }
  # creating difference and lag value
  eq11 <- dy1 ~ dx1
  eq12 <- dy2 ~ dx2
  eq13 <- dy3 ~ dx3
  eq14 <- dy4 ~ dx4
  eq15 <- dy5 ~ dx5
  eq21 <- dy1 ~ dx1 + x_lag1 + y_lag1
  eq22 <- dy2 ~ dx2 + x_lag2 + y_lag2
  eq23 <- dy3 ~ dx3 + x_lag3 + y_lag3
  eq24 <- dy4 ~ dx4 + x_lag4 + y_lag4
  eq25 <- dy5 ~ dx5 + x_lag5 + y_lag5

  # Combine them into a list for the systemfit function
  sur_data <- data.frame(dy1 = dy1, dx1 = dx1, x_lag1 = x_lag1, y_lag1 = y_lag1,
                        dy2 = dy2, dx2 = dx2, x_lag2 = x_lag2, y_lag2 = y_lag2,
                        dy3 = dy3, dx3 = dx3, x_lag3 = x_lag3, y_lag3 = y_lag3,
                        dy4 = dy4, dx4 = dx4, x_lag4 = x_lag4, y_lag4 = y_lag4,
                        dy5 = dy5, dx5 = dx5, x_lag5 = x_lag5, y_lag5 = y_lag5)
  system_eq_r <-  list(eq1 = eq11, eq2 = eq12, eq3 = eq13, eq4 = eq14, eq5 = eq15)
  system_eq_ur <- list(eq1 = eq21, eq2 = eq22, eq3 = eq23, eq4 = eq24, eq5 = eq25)

  # method = OLS
  ols_model <- systemfit(system_eq_ur, method = "OLS", data = sur_data)
  res_matrix <- residuals(ols_model)
  cor_matrix <- cor(res_matrix)
  upper_triangle   <- cor_matrix[upper.tri(cor_matrix)]
  mean_cor    <- mean(upper_triangle)
  abs_mean_cor <- mean(abs(upper_triangle))
  cat(sprintf("  Mean |correlation|     : %.4f\n", abs_mean_cor))

  # unrestricted model
  sur_model_ur <- systemfit(system_eq_ur, method = "SUR", data = sur_data,control = systemfit.control(methodResidCov = "noDfCor"))
  # Wald tests for H0: x_lag_j = 0 and y_lag_j = 0
  bounds1 <- wald_sur_eq(sur_model_ur, c("eq1_x_lag1", "eq1_y_lag1"))
  bounds2 <- wald_sur_eq(sur_model_ur, c("eq2_x_lag2", "eq2_y_lag2"))
  bounds3 <- wald_sur_eq(sur_model_ur, c("eq3_x_lag3", "eq3_y_lag3"))
  bounds4 <- wald_sur_eq(sur_model_ur, c("eq4_x_lag4", "eq4_y_lag4"))
  bounds5 <- wald_sur_eq(sur_model_ur, c("eq5_x_lag5", "eq5_y_lag5"))

  return(list(bounds1=bounds1,
              bounds2=bounds2,
              bounds3=bounds3,
              bounds4=bounds4,
              bounds5=bounds5))
}

#++++++++++++++++++++++++++++++++++++++++++++++++
# Step 3: assign data
#++++++++++++++++++++++++++++++++++++++++++++++++
for (j in 1:5) {
  # paste0("x", j) creates the name "x1", "x2", etc.
  # assign() takes that name and assigns the data to it in your environment
  # Grab the data for this specific loop iteration first
  temp_x <- data_set[, (j-1)*2+1]
  temp_y <- data_set[, (j-1)*2+2]

  # Now generate all the names and assign them!
  assign(paste0("x",     j), temp_x)
  assign(paste0("y",     j), temp_y)
  assign(paste0("x_lag", j), temp_x[-T_obs])
  assign(paste0("dx",    j), diff(temp_x))
  assign(paste0("y_lag", j), temp_y[-T_obs])
  assign(paste0("dy",    j), diff(temp_y))
}

#++++++++++++++++++++++++++++++++++++++++++++++++
# Step 4: Run the function
#++++++++++++++++++++++++++++++++++++++++++++++++

results   <- sur_ardl_bounds(x1,x2,x3,x4,x5,y1,y2,y3,y4,y5)
mout[1,1] <- results$bounds1
mout[2,1] <- results$bounds2
mout[3,1] <- results$bounds3
mout[4,1] <- results$bounds4
mout[5,1] <- results$bounds5

print(mout)
