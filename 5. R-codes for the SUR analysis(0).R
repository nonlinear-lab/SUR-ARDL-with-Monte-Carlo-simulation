#+++++++++++++++++++++++++++++++++
# Step 0: set a working directory
#+++++++++++++++++++++++++++++++
setwd("C:/R/ardl")
library(systemfit)

#+++++++++++++++++++++++++++++++++
# Step 1: Load Data
#+++++++++++++++++++++++++++++++
data_set <- read.csv("data(3).csv")
data_set <- data_set[, -1]
colnames(data_set) <- c("x1", "y1", "x2", "y2", "x3", "y3", "x4", "y4", "x5", "y5")
T_obs   <- nrow(data_set)
k       <- 1
mout    <- matrix(NA,nrow=5,ncol=1)

#++++++++++++++++++++++++++++++++++++++++++++++++
# Step 2: assign data
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

#++++++++++++++++++++++++++++++++++++++++++++++
# Step 2: Define the Equations
#++++++++++++++++++++++++++++++++++++++++++++++
# We define a separate formula for each pair
eq1 <- dy1 ~ dx1
eq2 <- dy2 ~ dx2
eq3 <- dy3 ~ dx3
eq4 <- dy4 ~ dx4
eq5 <- dy5 ~ dx5

# Combine them into a list for the systemfit function
sur_data <- data.frame(dy1 = dy1, dx1 = dx1, x_lag1 = x_lag1, y_lag1 = y_lag1,
                       dy2 = dy2, dx2 = dx2, x_lag2 = x_lag2, y_lag2 = y_lag2,
                       dy3 = dy3, dx3 = dx3, x_lag3 = x_lag3, y_lag3 = y_lag3,
                       dy4 = dy4, dx4 = dx4, x_lag4 = x_lag4, y_lag4 = y_lag4,
                       dy5 = dy5, dx5 = dx5, x_lag5 = x_lag5, y_lag5 = y_lag5)

system_eq <- list(eq1 = eq1, eq2 = eq2, eq3 = eq3, eq4 = eq4, eq5 = eq5)

#++++++++++++++++++++++++++++++++++++++++++++++
# Step 3: Diagnostic Estimation (OLS)
#++++++++++++++++++++++++++++++++++++++++++++++
# method = "SUR" applies Zellner's seemingly unrelated regressions approach
ols_model <- systemfit(system_eq, method = "OLS", data = sur_data, control = systemfit.control(methodResidCov = "noDfCor"))
res_matrix <- residuals(ols_model)
cor_matrix <- cor(res_matrix)
upper_triangle   <- cor_matrix[upper.tri(cor_matrix)]
mean_cor    <- mean(upper_triangle)
abs_mean_cor <- mean(abs(upper_triangle))

cat("\n--- Mean Cross-Equation Error Correlation ---\n")
cat(sprintf("  Number of unique pairs : %d\n",   length(upper_triangle)))
cat(sprintf("  Mean correlation       : %.4f\n", mean_cor))
cat(sprintf("  Mean |correlation|     : %.4f\n", abs_mean_cor))
cat(sprintf("  Min  correlation       : %.4f\n", min(upper_triangle)))
cat(sprintf("  Max  correlation       : %.4f\n", max(upper_triangle)))



# Print the comprehensive summary
print(ols_model)

#++++++++++++++++++++++++++++++++++++++++++++++
# Step 4: Estimate the SUR Model
#++++++++++++++++++++++++++++++++++++++++++++++
# Temporarily run standard OLS instead of SUR to bypass the Cholesky crash
sur_model <- systemfit(system_eq, method = "SUR", data = sur_data)
res1      <- residuals(sur_model)$eq1
ssr_r1     <- sum(res1^2)
res2      <- residuals(sur_model)$eq2
ssr_r2     <- sum(res2^2)
res3      <- residuals(sur_model)$eq3
ssr_r3     <- sum(res3^2)
res4      <- residuals(sur_model)$eq4
ssr_r4     <- sum(res1^4)
res5      <- residuals(sur_model)$eq5
ssr_r5     <- sum(res1^5)

# Print the comprehensive summary
print(sur_model)
