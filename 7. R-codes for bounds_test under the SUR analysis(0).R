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

#++++++++++++++++++++++++++++++++++++++++++++++
# Step 3: Diagnostic Estimation (OLS)
#++++++++++++++++++++++++++++++++++++++++++++++
# method = "SUR" applies Zellner's seemingly unrelated regressions approach
ols_model <- systemfit(system_eq_ur, method = "OLS", data = sur_data)
res_matrix <- residuals(ols_model)
cor_matrix <- cor(res_matrix)
upper_triangle   <- cor_matrix[upper.tri(cor_matrix)]
mean_cor    <- mean(upper_triangle)
abs_mean_cor <- mean(abs(upper_triangle))
cat(sprintf("  Mean |correlation|     : %.4f\n", abs_mean_cor))

# Print the comprehensive summary
print(ols_model)

#++++++++++++++++++++++++++++++++++++++++++++++
# Step 4: Estimate the SUR Model
#++++++++++++++++++++++++++++++++++++++++++++++
# Temporarily run standard OLS instead of SUR to bypass the Cholesky crash
sur_model_r <- systemfit(system_eq_r, method = "SUR", data = sur_data, control = systemfit.control(methodResidCov = "noDfCor"))
res11      <- residuals(sur_model_r)$eq1
ssr_r1     <- sum(res11^2)
res12      <- residuals(sur_model_r)$eq2
ssr_r2     <- sum(res12^2)
res13      <- residuals(sur_model_r)$eq3
ssr_r3     <- sum(res13^2)
res14      <- residuals(sur_model_r)$eq4
ssr_r4     <- sum(res14^2)
res15      <- residuals(sur_model_r)$eq5
ssr_r5     <- sum(res15^2)

sur_model_ur <- systemfit(system_eq_ur, method = "SUR", data = sur_data,control = systemfit.control(methodResidCov = "noDfCor"))
res21      <- residuals(sur_model_ur)$eq1
ssr_ur1     <- sum(res21^2)
res22      <- residuals(sur_model_ur)$eq2
ssr_ur2     <- sum(res22^2)
res23      <- residuals(sur_model_ur)$eq3
ssr_ur3     <- sum(res23^2)
res24      <- residuals(sur_model_ur)$eq4
ssr_ur4     <- sum(res24^2)
res25      <- residuals(sur_model_ur)$eq5
ssr_ur5     <- sum(res25^2)

df_num <- 2 # this is number of restriction, x_lag, y_lag
df_den <- nrow(sur_data) - (length(coef(sur_model_ur))/5)#Sample Size minus Number of Parameters (4)

bounds_stat1 <- ((ssr_r1 - ssr_ur1) / df_num) / (ssr_ur1 / df_den)
bounds_stat2 <- ((ssr_r2 - ssr_ur2) / df_num) / (ssr_ur2 / df_den)
bounds_stat3 <- ((ssr_r3 - ssr_ur3) / df_num) / (ssr_ur3 / df_den)
bounds_stat4 <- ((ssr_r4 - ssr_ur4) / df_num) / (ssr_ur4 / df_den)
bounds_stat5 <- ((ssr_r5 - ssr_ur5) / df_num) / (ssr_ur5 / df_den)

mout[1,1]=bounds_stat1
mout[2,1]=bounds_stat2
mout[3,1]=bounds_stat3
mout[4,1]=bounds_stat4
mout[5,1]=bounds_stat5

print(mout)
