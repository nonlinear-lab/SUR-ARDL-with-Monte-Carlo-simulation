#+++++++++++++++++++++++++++++++++++++++
# Step 0: dorking directory and library
#+++++++++++++++++++++++++++++++++++++++
setwd("C:/R/ardl")
# Install the package if not already installed
# install.packages("ARDL")

# Load the library
library(ARDL)

#++++++++++++++++++++++++++++++++++++++++++++++++
# Step 1: load data
#++++++++++++++++++++++++++++++++++++++++++++++++
my_data <- read.csv("data(2).csv")
my_data <- as.matrix(my_data)
T_obs   <- nrow(my_data)
k       <- 1
mout    <- matrix(NA,nrow=5,ncol=1)


#++++++++++++++++++++++++++++++++++++++++++++++++
# Step 3: Run the program
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
  temp_data   <- data.frame(y = y, x = x)
  ardl_model  <- ardl(y ~ x, data=temp_data, order = c(1, 1))
  print(ardl_model)
  bounds_test <- bounds_f_test(ardl_model, case = 3)
  print(bounds_test)
  mout[j,1]   <- bounds_test$statistic
}
print(mout)
