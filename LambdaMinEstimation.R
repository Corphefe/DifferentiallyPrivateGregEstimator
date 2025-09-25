library(tidyverse)

set.seed(144)
# -------------------- Parameters --------------------
N <- 10000
n <- 500
B_x <- 1
B_y <- 3
dist_choice <- "exponential"  # "normal", "uniform", or "exponential"
beta_true <- c(-1.44, 0.42)  # choosing betas that would keep most Y in [-B_y, B_y] given sigma_e
sigma_e <- 0.44
simulations <- 10000

# -------------------- Population --------------------
# Generate X according to chosen distribution
if (dist_choice == "normal") {
  x_pop <- rnorm(N, mean = 0, sd = 0.44)
} else if (dist_choice == "uniform") {
  x_pop <- runif(N, min = -B_x, max = B_x)
} else if (dist_choice == "exponential") {
  x_raw <- rexp(N, rate = 1)
  x_pop <- x_raw - mean(x_raw)  # center before clipping
} else {
  stop("Unknown dist_choice")
}

# Clip and generate Y
x_pop <- pmin(pmax(x_pop, -B_x), B_x)
y_pop <- beta_true[1] + beta_true[2] * x_pop + rnorm(N, mean = 0, sd = sigma_e)
y_pop <- pmin(pmax(y_pop, -B_y), B_y)

# Build augmented X with intercept
X_aug <- cbind(1, x_pop)  # [n x d] with intercept

mu_x <- colMeans(X_aug)

# --------------------- lambda_min estimation -------------------
lambda_mins <- numeric(simulations)

for (i in 1:simulations) {
  sample_indices <- sample(1:N, n, replace = FALSE)
  X_sample <- X_aug[sample_indices, , drop = FALSE]
  
  A <- t(X_sample) %*% X_sample
  eig_vals <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
  lambda_mins[i] <- min(eig_vals)
}

lambda_min_empirical <- min(lambda_mins)
lambda_min_empirical 



