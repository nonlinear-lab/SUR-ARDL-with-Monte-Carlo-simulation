#+++++++++++++++++++++++++++++++++
# Step 0: Load data
#+++++++++++++++++++++++++++++++
setwd("C:/R/ardl")
data_set <- read.csv("data(3).csv")
data_set <- data_set[, -1]                          # ADD THIS LINE
colnames(data_set) <- c("x1", "y1", "x2", "y2", "x3", "y3", "x4", "y4", "x5", "y5")

T_obs0 <- nrow(data_set)
T_obs <- T_obs0 - 1
M <- 5 # Number of equations
k <- 2 # Intercept + 1 slope variable

colnames(data_set) <- c("x1", "y1", "x2", "y2", "x3", "y3", "x4", "y4", "x5", "y5")
for (j in 1:5) {
  temp_x <- data_set[, (j-1)*2+1]
  temp_y <- data_set[, (j-1)*2+2]
  assign(paste0("x",     j), temp_x)
  assign(paste0("y",     j), temp_y)
  assign(paste0("x_lag", j), temp_x[-T_obs])
  assign(paste0("dx",    j), diff(temp_x))
  assign(paste0("y_lag", j), temp_y[-T_obs])
  assign(paste0("dy",    j), diff(temp_y))
}
#++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Step 1: First Stage OLS (To get residual covariance)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Updated formulas to match EXACT column names from the CSV
eq1 <- lm(dy1 ~ dx1)
eq2 <- lm(dy2 ~ dx2)
eq3 <- lm(dy3 ~ dx3)
eq4 <- lm(dy4 ~ dx4)
eq5 <- lm(dy5 ~ dx5)

# Stack residuals into a (T_obs x M) matrix
res_matrix <- cbind(residuals(eq1), residuals(eq2), residuals(eq3),
                    residuals(eq4), residuals(eq5))

# Calculate the Variance-Covariance matrix of the residuals (Sigma)
Sigma_e <- t(res_matrix) %*% res_matrix / (T_obs - k)

# Get the inverse of the covariance matrix
Sigma_inv <- solve(Sigma_e)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Step 2: STRICTLY formatted Matrices for SUR
#++++++++++++++++++++++++++++++++++++++++++++++++++++++
# X matrices use x1, x2, x3...
X_list <- list(
  matrix(c(rep(1, T_obs), dx1), nrow = T_obs, ncol = 2),
  matrix(c(rep(1, T_obs), dx2), nrow = T_obs, ncol = 2),
  matrix(c(rep(1, T_obs), dx3), nrow = T_obs, ncol = 2),
  matrix(c(rep(1, T_obs), dx4), nrow = T_obs, ncol = 2),
  matrix(c(rep(1, T_obs), dx5), nrow = T_obs, ncol = 2)
)

# Y matrices use the exact verbose names from your data
Y_list <- list(
  matrix(dy1, nrow = T_obs, ncol = 1),
  matrix(dy2, nrow = T_obs, ncol = 1),
  matrix(dy3, nrow = T_obs, ncol = 1),
  matrix(dy4, nrow = T_obs, ncol = 1),
  matrix(dy5, nrow = T_obs, ncol = 1)
)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Step 3: Zellner's GLS Calculation
#++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Initialize the 10x10 matrix and 10x1 vector with zeros
XT_Omega_X <- matrix(0, nrow = M * k, ncol = M * k)
XT_Omega_Y <- matrix(0, nrow = M * k, ncol = 1)

for (i in 1:M) {
  # 1. Fill the XT_Omega_Y vector
  sum_Y <- matrix(0, nrow = k, ncol = 1)
  for (j in 1:M) {
    sum_Y <- sum_Y + Sigma_inv[i, j] * (t(X_list[[i]]) %*% Y_list[[j]])
  }
  row_idx_Y <- ((i - 1) * k + 1):(i * k)
  XT_Omega_Y[row_idx_Y, 1] <- sum_Y

  # 2. Fill the XT_Omega_X matrix
  for (j in 1:M) {
    block_val <- Sigma_inv[i, j] * (t(X_list[[i]]) %*% X_list[[j]])
    row_idx_X <- ((i - 1) * k + 1):(i * k)
    col_idx_X <- ((j - 1) * k + 1):(j * k)
    XT_Omega_X[row_idx_X, col_idx_X] <- block_val
  }
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Step 4: Final Estimation (Updated Output)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Variance-Covariance matrix of the SUR estimators: (X' \Omega^-1 X)^-1
V_sur <- solve(XT_Omega_X)

# Final SUR Coefficients: (X' \Omega^-1 X)^-1 * (X' \Omega^-1 Y)
beta_sur <- V_sur %*% XT_Omega_Y

# Extract Standard Errors
se_sur <- sqrt(diag(V_sur))

# Compute t-statistics
t_stats <- beta_sur / se_sur

# Print the raw coefficients just to guarantee the math worked!
print("Raw Coefficients (beta_sur):")
print(beta_sur)

# Force the matrices into flat numeric vectors so data.frame() doesn't choke
beta_vec <- as.numeric(beta_sur)
se_vec   <- as.numeric(se_sur)
t_vec    <- as.numeric(t_stats)

# Format output nicely
results <- data.frame(
  Equation  = rep(paste0("Eq", 1:5), each = 2),
  Term      = rep(c("Intercept", "Slope"), 5),
  Estimate  = round(beta_vec, 4),
  Std_Error = round(se_vec, 4),
  t_value   = round(t_vec, 4)
)

print("--- Formatted Manual SUR Estimation Results ---")
print(results)
