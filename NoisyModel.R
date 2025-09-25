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
  lambda_min <- 72
} else if (dist_choice == "uniform") {
  x_pop <- runif(N, min = -B_x, max = B_x)
  lambda_min <- 140
} else if (dist_choice == "exponential") {
  x_raw <- rexp(N, rate = 1)
  x_pop <- x_raw - mean(x_raw)  # center before clipping
  lambda_min <- 175
} else {
  stop("Unknown dist_choice")
}

# Clip and generate Y
x_pop <- pmin(pmax(x_pop, -B_x), B_x)
y_pop <- beta_true[1] + beta_true[2] * x_pop + rnorm(N, mean = 0, sd = sigma_e)
y_pop <- pmin(pmax(y_pop, -B_y), B_y)

X_aug <- cbind(1, x_pop)  # [n x d] with intercept
mu_x <- colMeans(X_aug)

# ------------------------ Privacy ---------------------------
# Target privacy
rho <- 0.045
rho_stats <- rho / 2 

# -------------------- Sensitivities & sigmas (zCDP) --------------------
R_x <- sqrt(1 + B_x^2)
B_e <- 2*B_y

Delta_beta <- 2 * R_x * B_e / lambda_min
term1 <- ( (B_x^2 * B_e) / (lambda_min - B_x^2) )^2
term2 <- (4 * B_e^2 * B_x^2) / (lambda_min - B_x^2)
term3 <- (1 / (n - 1)) * B_e^2
Delta_V <- (1 / n) * (1 - n / N) * (term1 + term2 + term3)

sigma2_beta  <- Delta_beta^2  / (2 * rho_stats)
sigma2_V <- Delta_V^2 / (2 * rho_stats)

# -------------------- Helper: variance correction (delta method) --------------------
variance_correction <- function(sigma2_beta) {
   (sum(mu_x^2) * sigma2_beta)
}

# -------------------- DP GREG NoisyModel (A2) --------------------
dp_greg_run_A2 <- function(X_aug, y_pop, n, sigma2_beta, sigma2_V, lambda_min) {
  N <- length(x_pop)
  idx <- sample(1:N, n, replace = FALSE)
  X_s <- X_aug[idx, , drop = FALSE]
  y_s <- y_pop[idx]
  
  xx_bar <- crossprod(X_s) / n       # (d x d) = (X^T X)/n
  xy_bar <- (t(X_s) %*% y_s) / n     # (d x 1) = (X^T y)/n
  
  beta_hat <- solve(xx_bar, xy_bar)  # solves (xx_bar) * beta = xy_bar
  beta_tilde <- beta_hat + rnorm(length(beta_hat), mean = 0, sd = sqrt(sigma2_beta))
  
  mu_hat <- t(mu_x) %*% beta_hat
  mu_hat <- as.numeric(mu_hat)
  
  mu_tilde <- t(mu_x) %*% beta_tilde
  mu_tilde <- as.numeric(mu_tilde)
  
  e_s_hat <- y_s - (X_s) %*% beta_hat
  
  V_hat <- (1/n) * (1 - n/N) * (1/ (n-1)) * sum(e_s_hat^2)
  V_tilde <- V_hat + rnorm(1, mean = 0, sd = sqrt(sigma2_V)) + sum(mu_x^2) * sigma2_beta
  
  list(
    mu_hat = mu_hat, 
    mu_tilde = mu_tilde,
    V_hat = V_hat, 
    V_tilde = V_tilde
  )
  
}

# -------------------- Monte Carlo --------------------
mu_hat_vals      <- numeric(simulations)
mu_tilde_vals    <- numeric(simulations)
V_hat_vals       <- numeric(simulations)
V_tilde_vals     <- numeric(simulations)

for (i in seq_len(simulations)) {
  out <- dp_greg_run_A2(X_aug = X_aug, y_pop = y_pop, n = n, sigma2_beta = sigma2_beta, sigma2_V = sigma2_V, lambda_min = lambda_min)
  mu_hat_vals[i]      <- out$mu_hat
  mu_tilde_vals[i]    <- out$mu_tilde
  V_hat_vals[i]       <- out$V_hat
  V_tilde_vals[i]     <- out$V_tilde
}

results_df <- tibble(
  mu_hats   = mu_hat_vals,
  mu_tildes = mu_tilde_vals,
  V_hats    = V_hat_vals,
  V_tildes  = V_tilde_vals
)

# -------- True mean of Y --------
mu_y <- mean(y_pop)

# ---------- 1) POINT ESTIMATOR SUMMARY ----------
# --- Empirical moments for point estimators ---
emp_var_np <- var(results_df$mu_hats,    na.rm = TRUE)
emp_var_dp <- var(results_df$mu_tildes,  na.rm = TRUE)

bias_np <- mean(results_df$mu_hats   - mu_y, na.rm = TRUE)
bias_dp <- mean(results_df$mu_tildes - mu_y, na.rm = TRUE)

mse_np <- bias_np^2 + emp_var_np
mse_dp <- bias_dp^2 + emp_var_dp

point_tbl <- tibble::tibble(
  estimator         = c("NP", "DP"),
  mean_estimate     = c(mean(results_df$mu_hats,   na.rm = TRUE),
                        mean(results_df$mu_tildes, na.rm = TRUE)),
  bias              = c(bias_np, bias_dp),
  emp_var_estimates = c(emp_var_np, emp_var_dp),
  mse               = c(mse_np, mse_dp),
  # Relative-to-NP metrics (NP row is 1 or 0 by definition)
  percent_rel_bias  = c(bias_np / (abs(mu_y)) * 100, bias_dp / (abs(mu_y)) * 100),
  rvar_vs_np        = c(1, emp_var_dp / emp_var_np),
  rmse_vs_np        = c(1, mse_dp / mse_np),
)

# ---------- 2) VARIANCE ESTIMATOR SUMMARY ----------
# Empirical variance of the point estimators (targets for calibration)
emp_var_mu_np <- emp_var_np
emp_var_mu_dp <- emp_var_dp

mean_V_np <- mean(results_df$V_hats,   na.rm = TRUE)
mean_V_dp <- mean(results_df$V_tildes, na.rm = TRUE)

calib_np <- mean_V_np / emp_var_mu_np   # should be ~1 if well-calibrated
calib_dp <- mean_V_dp / emp_var_mu_dp

biasV_np <- mean_V_np - emp_var_mu_np
biasV_dp <- mean_V_dp - emp_var_mu_dp

varV_np  <- var(results_df$V_hats,   na.rm = TRUE)
varV_dp  <- var(results_df$V_tildes, na.rm = TRUE)

var_tbl <- tibble::tibble(
  estimator        = c("NP", "DP"),
  mean_V           = c(mean_V_np, mean_V_dp),
  emp_var_mu       = c(emp_var_mu_np, emp_var_mu_dp),
  calibration_kappa= c(calib_np, calib_dp),
  bias_V           = c(biasV_np, biasV_dp),
  var_of_V         = c(varV_np, varV_dp),
)

# ---------- 3) CONFIDENCE INTERVAL SUMMARY ----------
z <- qnorm(0.975)

# Guard against nonpositive/NA variances
valid_np <- is.finite(results_df$V_hats)   & (results_df$V_hats   > 0)
valid_dp <- is.finite(results_df$V_tildes) & (results_df$V_tildes > 0)

# NP CIs
ci_lower_np <- results_df$mu_hats[valid_np] - z * sqrt(results_df$V_hats[valid_np])
ci_upper_np <- results_df$mu_hats[valid_np] + z * sqrt(results_df$V_hats[valid_np])
coverage_np <- mean(ci_lower_np <= mu_y & mu_y <= ci_upper_np)
width_np    <- mean(ci_upper_np - ci_lower_np)

Z_np   <- (results_df$mu_hats[valid_np] - mu_y) / sqrt(results_df$V_hats[valid_np])
meanZ_np <- mean(Z_np)
sdZ_np   <- sd(Z_np)
tail_np  <- mean(abs(Z_np) > z)
n_np     <- sum(valid_np)

# DP CIs
ci_lower_dp <- results_df$mu_tildes[valid_dp] - z * sqrt(results_df$V_tildes[valid_dp])
ci_upper_dp <- results_df$mu_tildes[valid_dp] + z * sqrt(results_df$V_tildes[valid_dp])
coverage_dp <- mean(ci_lower_dp <= mu_y & mu_y <= ci_upper_dp)
width_dp    <- mean(ci_upper_dp - ci_lower_dp)

Z_dp   <- (results_df$mu_tildes[valid_dp] - mu_y) / sqrt(results_df$V_tildes[valid_dp])
meanZ_dp <- mean(Z_dp)
sdZ_dp   <- sd(Z_dp)
tail_dp  <- mean(abs(Z_dp) > z)
n_dp     <- sum(valid_dp)

ci_tbl <- tibble::tibble(
  estimator        = c("NP", "DP"),
  coverage         = c(coverage_np, coverage_dp),
  avg_width        = c(width_np, width_dp),
  rwidth_vs_np     = c(1, width_dp / width_np),
  mean_Z           = c(meanZ_np, meanZ_dp),
  sd_Z             = c(sdZ_np, sdZ_dp),
  tail_rate_gt_1.96= c(tail_np, tail_dp),
  n_used           = c(n_np, n_dp)
)

# View
point_tbl_disp <- point_tbl %>% mutate(across(where(is.numeric), ~round(., 6)))
var_tbl_disp   <- var_tbl   %>% mutate(across(where(is.numeric), ~round(., 6)))
ci_tbl_disp    <- ci_tbl    %>% mutate(across(where(is.numeric), ~round(., 6)))
point_tbl_disp
var_tbl_disp
ci_tbl_disp