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

mu_x <- mean(x_pop)  # known population mean (public)

# ------------------------ Privacy ---------------------------
# Target privacy
rho <- 0.045
rho_stats <- rho / 5  

# -------------------- Sensitivities & sigmas (zCDP) --------------------
Delta_x   <- 2*B_x
Delta_y   <- 2*B_y
Delta_xx  <- B_x^2
Delta_xy  <- 2*B_x*B_y
Delta_yy  <- B_y^2

sigma2_x  <- (Delta_x^2)  / (2 * rho_stats * n^2)
sigma2_y  <- (Delta_y^2)  / (2 * rho_stats * n^2)
sigma2_xx <- (Delta_xx^2) / (2 * rho_stats * n^2)
sigma2_xy <- (Delta_xy^2) / (2 * rho_stats * n^2)
sigma2_yy <- (Delta_yy^2) / (2 * rho_stats * n^2)

# -------------------- Helper: variance correction (delta method) --------------------
var_correction <- function(a, b, c, d) {
  D <- c - a^2
  A <- d - a*b
  beta1_hat_tilde <- A / D
  
  ### 2. gradients
  g_a <- -beta1_hat_tilde + (mu_x - a) * ((-b*D + 2*a*A) / D^2)
  g_b <-  1 - (mu_x - a) * (a / D)
  g_c <- -(mu_x - a) * (A / D^2)
  g_d <-  (mu_x - a) / D
  
  ### 3. DP-noise variances (these are what you used in rnorm())
  var_a <- sigma2_x
  var_b <- sigma2_y
  var_c <- sigma2_xx
  var_d <- sigma2_xy
  
  var_dp_mu <- g_a^2 * var_a +
    g_b^2 * var_b +
    g_c^2 * var_c +
    g_d^2 * var_d
  
  (var_dp_mu)
}

# -------------------- DP GREG NoisyStats (A1) --------------------
dp_greg_run_A1 <- function(x_pop, y_pop, n, sigma2_x, sigma2_y, sigma2_xx, sigma2_xy, sigma2_yy) {
  N <- length(x_pop)
  idx <- sample(1:N, n, replace = FALSE)
  x_s <- x_pop[idx]
  y_s <- y_pop[idx]
  
  # --- True sufficient statistics (means) ---
  x_bar  <- mean(x_s)
  y_bar  <- mean(y_s)
  xx_bar <- mean(x_s^2)
  xy_bar <- mean(x_s * y_s)
  yy_bar <- mean(y_s^2)
  
  # --- Noisy sufficient statistics ---
  x_bar_tilde  <- x_bar  + rnorm(1, 0, sqrt(sigma2_x))
  y_bar_tilde  <- y_bar  + rnorm(1, 0, sqrt(sigma2_y))
  xx_bar_tilde <- xx_bar + rnorm(1, 0, sqrt(sigma2_xx))
  xy_bar_tilde <- xy_bar + rnorm(1, 0, sqrt(sigma2_xy))
  yy_bar_tilde <- yy_bar + rnorm(1, 0, sqrt(sigma2_yy))
  
  # --- Compute betas without privacy ---
  beta1_hat <- (xy_bar - x_bar * y_bar) / 
    (xx_bar - x_bar^2)
  beta0_hat <- y_bar - beta1_hat * x_bar
  
  # --- Compute betas with privacy ---
  beta1_tilde <- (xy_bar_tilde - x_bar_tilde * y_bar_tilde) /
    (xx_bar_tilde - x_bar_tilde^2)
  beta0_tilde <- y_bar_tilde - beta1_tilde * x_bar_tilde
  
  greg_hat <- beta0_hat + mu_x * beta1_hat
  greg_tilde <- beta0_tilde + mu_x * beta1_tilde
  
  V_hat <- (1/n) * (1 - (n/N)) * (1/(n-1)) * (
    n * (yy_bar - 2*beta0_hat*y_bar - 2*beta1_hat*xy_bar + beta0_hat^2 + 2*beta1_hat*beta0_hat*x_bar + beta1_hat*xx_bar*beta1_hat ) - 
      n * (y_bar - beta0_hat - beta1_hat * x_bar)^2
  )
  V_hat_tilde <- (1/n) * (1 - (n/N)) * (1/(n-1)) * (
    n * (yy_bar_tilde - 2*beta0_tilde*y_bar_tilde - 2*beta1_tilde*xy_bar_tilde + beta0_tilde^2 + 2*beta1_tilde*beta0_tilde*x_bar_tilde + beta1_tilde*xx_bar_tilde*beta1_tilde ) - 
      n * (y_bar_tilde - beta0_tilde - beta1_tilde * x_bar_tilde)^2
  )
  
  # Additive DP-noise variance via delta method
  var_dp_mu <- var_correction(
    a = x_bar_tilde, b = y_bar_tilde, c = xx_bar_tilde, d = xy_bar_tilde
  )
  
  V_tilde <- V_hat_tilde + var_dp_mu
  
  list(
    mu_hat = greg_hat, 
    mu_tilde = greg_tilde,
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
  out <- dp_greg_run_A1(
    x_pop = x_pop, y_pop = y_pop, n = n,
    sigma2_x = sigma2_x, sigma2_y = sigma2_y,
    sigma2_xx = sigma2_xx, sigma2_xy = sigma2_xy, sigma2_yy = sigma2_yy
  )
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