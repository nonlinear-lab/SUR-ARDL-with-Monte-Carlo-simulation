#+++++++++++++++++++++++++++++++++++
# Step 1: set monte carlo simulation
#+++++++++++++++++++++++++++++++++++
# Monte Carlo demonstration of SUR efficiency gain
library(systemfit)
set.seed(123)
reps      <- 1000
T         <- 200
rho       <- 0.8
true_beta <- 2.5
T_obs     <- T
mout <- matrix(NA, nrow = reps, ncol = 5)
colnames(mout) <- c("AR=0","AR=0.5","AR=0,8","AR=0.95","AR=1")

#++++++++++++++++++++++++++++++++++++++++++++++
# Step 2: wald function
#++++++++++++++++++++++++++++++++++++++++++++++
wald_sur_eq <- function(model, restriction) {
  b <- coef(model)
  V <- vcov(model)
  idx <- match(restriction, names(b))
  R <- matrix(0, nrow = length(idx), ncol = length(b))
  for (i in seq_along(idx)) {
    R[i, idx[i]] <- 1
  }

  Rb <- R %*% b
  RVRT <- R %*% V %*% t(R)
  W <- as.numeric(t(Rb) %*% solve(RVRT) %*% Rb)
  F_version <- W / length(idx)
  return(F_version)
}

#++++++++++++++++++++++++++++++++++++++++++++++
# Step 2: SUR_ARDL bounds test function
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

  # unrestricted model
  sur_model_ur <- systemfit(system_eq_ur, method = "SUR", data = sur_data,control = systemfit.control(methodResidCov = "noDfCor"))

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

#+++++++++++++++++++++++++++++++++++
# Step 3 simulation loop (1000 reps)
#+++++++++++++++++++++++++++++++++++
for (s in 1:reps) {
  if (s %% 100 == 0) {
    cat("Now replication is", s, "out of", reps, "\n")
    flush.console()
  }
  common_y <- rnorm(T)
  for (j in 1:5) {
    # Generate a different x for each equation
    u_x <- rnorm(T)
    x <- cumsum(u_x)

    # Generate a different y for each equation
    # but allow cross-equation correlation through common_y
    u_y <- sqrt(rho) * common_y + sqrt(1 - rho) * rnorm(T)
    y <- cumsum(u_y)

    assign(paste0("x", j), x)
    assign(paste0("y", j), y)
  }

  results <- sur_ardl_bounds(x1, x2, x3, x4, x5,
                             y1, y2, y3, y4, y5)

  mout[s, 1] <- results$bounds1
  mout[s, 2] <- results$bounds2
  mout[s, 3] <- results$bounds3
  mout[s, 4] <- results$bounds4
  mout[s, 5] <- results$bounds5
}

#++++++++++++++++++++++++++++++++++++++++++++++++
# Step 4: Report the Findings (FIXED DIMNAMES)
#++++++++++++++++++++++++++++++++++++++++++++++++
# Define upper-tail significance levels
alpha_levels <- c(0.01, 0.025, 0.05, 0.10, 0.50, 0.90, 0.95, 0.975, 0.99)

# Calculate the quantiles for all 5 equations at once
# apply(matrix, 2, ...) means "calculate this for every column"
crit_vals <- apply(mout, 2, quantile, probs = alpha_levels, na.rm = TRUE)

# Convert to a clean data frame
results_table <- data.frame(crit_vals)

# FIX: Exactly 5 column names for the 5 equations
colnames(results_table) <- c("S1", "S2", "S3", "S4", "S5")

# Add row-wise mean of the displayed critical values
results_table$Mean <- rowMeans(results_table[, c("S1", "S2", "S3", "S4", "S5")],
                               na.rm = TRUE)
# Exactly 4 row names for the 4 significance levels
rownames(results_table) <- c("1%", "2.5%", "5%", "10%", "50%",
                             "90%", "95%", "97.5%", "99%")



print("Monte Carlo Simulation Complete! SUR-ARDL Critical Values:")
print(round(results_table, 3))
