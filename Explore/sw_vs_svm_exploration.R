## ============================================================#
### WORKING DIRECTORY
## ============================================================#
#setwd("~/Desktop/OSU/Research/Pretest-Simulation/User_framework_Rpkg/ML_application/SW_vs_ML_Application_eploration")
#' ============================================================
#' REPRODUCIBLE RNG SETUP (cross-platform, cross-version safe)
#' ============================================================
RNGversion("4.2.0")
set.seed(
  12345,
  kind        = "Mersenne-Twister",
  normal.kind = "Inversion",
  sample.kind = "Rejection"
)

## `%||%` null-coalescing operator — returns rhs when lhs is NULL
`%||%` <- function(lhs, rhs) if (!is.null(lhs)) lhs else rhs

##' ============================================================#
##' PACKAGE INSTALLATION & LOADING
##' ============================================================#

#' Use a stable CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
#' Detect OS
is_windows <- .Platform$OS.type == "windows"
#' Core package list 
pkgs <- c("MASS", "tidyverse","foreach", "doParallel", "doRNG", "doSNOW","nortest", "dgof", 
          "goftest","DescTools", "moments", "tseries","LaplacesDemon", "VGAM", "evd","Rfit",
          "lawstat", "coin","ineq","caret", "glmnet", "randomForest", "gbm","nnet", "mlbench", 
          "kernlab", "e1071","pROC", "ROCR","fractaldim","openxlsx", "Lmoments"
)

#' Add doMC only on Unix-like systems (Linux / macOS)
if (!is_windows) pkgs <- c(pkgs, "doMC")
#' Install any missing packages and attach all
pacman::p_load(char = pkgs)


#' ============================================================
#' SOURCE EXTERNAL ML NORMALITY FRAMEWORK
#' ============================================================
#' The ML framework is sourced into its own environment to avoid
#' overwriting user-framework functions with the same names, such as:
#' ============================================================

ml_env <- new.env(parent = globalenv())

source(
  file  = "~/Desktop/OSU/Research/Pretest-Simulation/User_framework_Rpkg/ML_application/ML_Approach_application_v2.R",
  local = ml_env
)


#' ============================================================================
#' 1. Generate standardized data with flexible distribution options
#' ============================================================================

generate_data <- function(n, dist, par = NULL, center_by = "median", standardize = TRUE) {
  
  # Input validation
  if (!is.numeric(n) || length(n) != 1 || n <= 0 || n != round(n)) {
    stop("n must be a single positive integer.")
  }
  
  center_by <- tolower(trimws(center_by))
  if (!center_by %in% c("median", "mean")) {
    stop("center_by must be 'median' or 'mean'.")
  }
  
  # Initialize storage:These will be populated below based on the distribution chosen
  samples <- NULL
  theoretical_mean <- theoretical_median <- sigma <- quartile_deviation <- NULL
  dist_label <- NULL
  
  # CASE 1: dist is a function
  if (is.function(dist)) {
    samples <- dist(n)
    theoretical_mean <- mean(samples)
    theoretical_median <- median(samples)
    sigma <- sd(samples)
    quartile_deviation <- IQR(samples)
    dist_label <- "custom_function"
  }
  
  # CASE 2: dist is a numeric vector (empirical)
  else if (is.numeric(dist) && length(dist) > 0) {
    samples <- sample(dist, size = n, replace = TRUE)
    theoretical_mean <- mean(dist)
    theoretical_median <- median(dist)
    sigma <- sd(dist)
    quartile_deviation <- IQR(dist)
    dist_label <- "empirical"
  }
  
  # CASE 3: dist is a character string (built-in distribution)
  else if (is.character(dist)) {
    dist <- tolower(trimws(dist))
    
    # Built-in distributions
    if (dist == "normal") {
      if (is.null(par)) par <- c(0, 1)
      samples <- rnorm(n, mean = par[1], sd = par[2])
      theoretical_mean <- par[1]
      theoretical_median <- par[1]
      sigma <- par[2]
      quartile_deviation <- qnorm(0.75, par[1], par[2]) - qnorm(0.25, par[1], par[2])
      
    } else if (dist == "chi_square") {
      if (is.null(par)) par <- 3
      samples <- rchisq(n, df = par)
      theoretical_mean <- par
      theoretical_median <- qchisq(0.5, df = par)
      sigma <- sqrt(2 * par)
      quartile_deviation <- qchisq(0.75, df = par) - qchisq(0.25, df = par)
      
    } else if (dist == "gamma") {
      if (is.null(par)) par <- c(3, 0.1)
      samples <- stats::rgamma(n, shape = par[1], rate = par[2])
      theoretical_mean <- par[1] / par[2]
      theoretical_median <- stats::qgamma(0.5, shape = par[1], rate = par[2])
      sigma <- sqrt(par[1]) / par[2]
      quartile_deviation <- qgamma(0.75, shape = par[1], rate = par[2]) - qgamma(0.25, shape = par[1], rate = par[2])
      
    } else if (dist == "exponential") {
      if (is.null(par)) par <- 1
      samples <- rexp(n, rate = par)
      theoretical_mean <- 1 / par
      theoretical_median <- log(2) / par
      sigma <- 1 / par
      quartile_deviation <- qexp(0.75, rate = par) - qexp(0.25, rate = par)
      
    } else if (dist == "f") {
      if (is.null(par)) par <- c(6, 15)
      samples <- rf(n, df1 = par[1], df2 = par[2])
      theoretical_mean <- if (par[2] > 2) par[2] / (par[2] - 2) else NA
      theoretical_median <- qf(0.5, df1 = par[1], df2 = par[2])
      sigma <- if (par[2] > 4) {
        sqrt(2 * par[2]^2 * (par[1] + par[2] - 2) / (par[1] * (par[2] - 2)^2 * (par[2] - 4)))
      } else NA
      quartile_deviation <- qf(0.75, df1 = par[1], df2 = par[2]) - qf(0.25, df1 = par[1], df2 = par[2])
      
    } else if (dist == "t") {
      if (is.null(par)) par <- 3
      samples <- rt(n, df = par)
      theoretical_mean <- if (par > 1) 0 else NA
      theoretical_median <- 0
      sigma <- if (par > 2) sqrt(par / (par - 2)) else NA
      quartile_deviation <- qt(0.75, df = par) - qt(0.25, df = par)
      
    } else if (dist == "uniform") {
      if (is.null(par)) par <- c(0, 1)
      samples <- runif(n, min = par[1], max = par[2])
      theoretical_mean <- (par[1] + par[2]) / 2
      theoretical_median <- (par[1] + par[2]) / 2
      sigma <- (par[2] - par[1]) / sqrt(12)
      quartile_deviation <- (par[2] - par[1]) / 2
      
    } else if (dist == "laplace") {
      if (is.null(par)) par <- c(2, 7)
      samples <- LaplacesDemon::rlaplace(n, location = par[1], scale = par[2])
      theoretical_mean <- par[1]
      theoretical_median <- par[1]
      sigma <- sqrt(2) * par[2]
      quartile_deviation <- par[2] * log(4)
      
    } else if (dist == "cauchy") {
      if (is.null(par)) par <- c(0, 1)
      samples <- rcauchy(n, location = par[1], scale = par[2])
      theoretical_mean <- NA
      theoretical_median <- par[1]
      sigma <- NA
      quartile_deviation <- 2 * par[2]
      
    } else if (dist == "gumbel") {
      if (is.null(par)) par <- c(0, 1)
      samples <- evd::rgumbel(n, loc = par[1], scale = par[2])
      theoretical_mean <- par[1] + par[2] * 0.5772156649 # Euler-Mascheroni constant
      theoretical_median <- par[1] - par[2] * log(log(2))
      sigma <- (pi * par[2]) / sqrt(6)
      quartile_deviation <- par[2] * (log(log(4/3)) - log(log(4)))
      
    } else if (dist == "weibull") {
      if (is.null(par)) par <- c(1, 2)
      samples <- rweibull(n, shape = par[1], scale = par[2])
      theoretical_mean <- par[2] * gamma(1 + 1 / par[1])
      theoretical_median <- par[2] * log(2)^(1 / par[1])
      sigma <- par[2] * sqrt(gamma(1 + 2 / par[1]) - gamma(1 + 1 / par[1])^2)
      quartile_deviation <- par[2] * ((-log(0.25))^(1 / par[1]) - (-log(0.75))^(1 / par[1]))
      
    } else if (dist == "lognormal") {
      if (is.null(par)) par <- c(0, 1)
      samples <- rlnorm(n, meanlog = par[1], sdlog = par[2])
      theoretical_mean <- exp(par[1] + par[2]^2 / 2)
      theoretical_median <- exp(par[1])
      sigma <- sqrt((exp(par[2]^2) - 1) * exp(2 * par[1] + par[2]^2))
      quartile_deviation <- exp(par[1]) * (exp(qnorm(0.75) * par[2]) - exp(qnorm(0.25) * par[2]))
      
    } else if (dist == "beta") {
      if (is.null(par)) par <- c(2, 5)
      samples <- rbeta(n, shape1 = par[1], shape2 = par[2])
      theoretical_mean <- par[1] / (par[1] + par[2])
      theoretical_median <- qbeta(0.5, shape1 = par[1], shape2 = par[2])
      sigma <- sqrt((par[1] * par[2]) /((par[1] + par[2])^2 * (par[1] + par[2] + 1)))
      quartile_deviation <- qbeta(0.75, shape1 = par[1], shape2 = par[2]) - qbeta(0.25, shape1 = par[1], shape2 = par[2])
      
    } else if (dist == "logistic") {
      if (is.null(par)) par <- c(0, 1)
      samples <- rlogis(n, location = par[1], scale = par[2])
      theoretical_mean <- par[1]
      theoretical_median <- par[1]
      sigma <- par[2] * pi / sqrt(3)
      quartile_deviation <- par[2] * log(3)
      
    } else if (dist == "pareto") {
      if (is.null(par)) par <- c(3, 1)  # shape (alpha), scale (xm)
      alpha <- par[1]; xm <- par[2]
      samples            <- VGAM::rpareto(n, scale = xm, shape = alpha)
      theoretical_mean   <- if (alpha > 1) xm * alpha / (alpha - 1) else NA
      theoretical_median <- xm * 2^(1/alpha)
      sigma              <- if (alpha > 2) xm * sqrt(alpha / ((alpha-1)^2 * (alpha-2))) else NA
      # Pareto quantiles: Q(p) = xm * (1-p)^(-1/alpha)
      quartile_deviation <- xm * (0.25)^(-1/alpha) - xm * (0.75)^(-1/alpha)
      
    } else if (dist == "contaminated") {
      if (is.null(par)) par <- c(0.65, 0, 1, 5)
      component_indicator <- rbinom(n, size = 1, prob = par[1])
      sd_draw <- ifelse(component_indicator == 1, par[3], par[4])
      samples <- rnorm(n, mean = par[2], sd = sd_draw)
      theoretical_mean <- par[2]
      theoretical_median <- par[2]
      sigma <- sqrt(par[1] * par[3]^2 + (1 - par[1]) * par[4]^2)
      quartile_deviation <- if (is.finite(sigma) && sigma > 0) {
        qnorm(0.75, 0, sigma) - qnorm(0.25, 0, sigma)
      } else NA
      
    } else {
      stop("Unsupported distribution: '", dist, "'")
    }
    
    dist_label <- dist
  }
  
  else {
    stop("dist must be a function, numeric vector, or character string")
  }
  
  # Attach basic attributes
  attr(samples, "distribution") <- dist_label
  
  # Return early if no standardization
  if (!standardize) {
    attr(samples, "standardized") <- FALSE
    return(samples)
  }
  
  # Standardization
  center_value <- if (center_by == "median") theoretical_median else theoretical_mean
  
  if (is.null(center_value) || !is.finite(center_value)) {
    warning("Cannot center by ", center_by, ". Falling back to median.")
    center_value <- theoretical_median
  }
  
  # center by center_value
  samples <- samples - center_value
  
  # Scaling
  if (is.finite(sigma) && sigma > 0) {
    scale_factor <- sigma
    scale_method <- "standard_deviation"
  } else if (is.finite(quartile_deviation) && quartile_deviation > 0) {
    scale_factor <- quartile_deviation
    scale_method <- "quartile_deviation"
    warning("Scaling by quartile deviation (SD not available/finite)")
  } else {
    scale_factor <- 1
    scale_method <- "none"
    warning("No valid scale found - returning centered data only")
  }
  
  # scale by the scale_factor
  samples <- samples / scale_factor
  
  # Attach metadata as attributes
  attr(samples, "center_by") <- center_by
  attr(samples, "center_value") <- center_value
  attr(samples, "sigma") <- sigma
  attr(samples, "quartile_deviation") <- quartile_deviation
  attr(samples, "standardized") <- TRUE
  attr(samples, "scale_factor") <- scale_factor
  attr(samples, "scale_method") <- scale_method
  
  return(samples)
}


#' ============================================================================
#' 2. Perform normality/diagnostic tests on a numeric sample
#' ============================================================================

generate_tests <- function(x, test, mu = 0, sigma = 1) {
  
  # Input validation
  if (!is.numeric(x)) stop("x must be a numeric vector.")
  if (length(x) < 3) stop("x must have at least 3 observations.")
  
  # CASE 1: test is a function
  if (is.function(test)) {
    result <- test(x)
    if (!is.numeric(result) || length(result) != 1) {
      stop("Test function must return a single named numeric scalar")
    }
    return(result)
  }
  
  # CASE 2: test is a character string
  if (!is.character(test)) {
    stop("test must be a character string or function")
  }
  
  # enforce upper case test names
  test <- toupper(trimws(test))
  
  # Built-in tests
  result <- switch(test,
                   # Shapiro-Wilk
                   SW = c(p.value = shapiro.test(x)$p.value),
                   # Shapiro-Francia
                   SF = c(p.value = nortest::sf.test(x)$p.value),
                   # Lilliefors
                   LF = c(p.value = nortest::lillie.test(x)$p.value),
                   # Kolmogorov-Smirnov
                   KS = c(p.value = ks.test(x, "pnorm", mean = mu, sd = sigma)$p.value),
                   # Anderson-Darling
                   AD = c(p.value = nortest::ad.test(x)$p.value),
                   # Anderson-Darling with specified parameters
                   AD2 = c(p.value = DescTools::AndersonDarlingTest(x, null = "pnorm", mean = mu, sd = sigma)$p.value),
                   # Cramér-von Mises
                   CVM = c(p.value = nortest::cvm.test(x)$p.value),
                   # Jarque-Bera
                   JB = c(p.value = tseries::jarque.bera.test(x)$p.value),
                   # D'Agostino-Pearson
                   DAP = c(p.value = moments::agostino.test(x)$p.value),
                   # Anscombe-Glynn
                   ANS = c(p.value = moments::anscombe.test(x)$p.value),
                   # Skewness z-test
                   SKEW = {
                     s <- moments::skewness(x)
                     z <- s / sqrt(6 / length(x))
                     c(p.value = 2 * (1 - pnorm(abs(z))))
                   },
                   # Kurtosis z-test
                   KURT = {
                     k <- moments::kurtosis(x)
                     z <- (k - 3) / sqrt(24 / length(x))
                     c(p.value = 2 * (1 - pnorm(abs(z))))
                   },
                   # Third central moment
                   MOM3 = c(moment3 = mean((x - mean(x))^3)),
                   # Fourth central moment
                   MOM4 = c(moment4 = mean((x - mean(x))^4)),
                   # Default: unrecognized test
                   stop(paste("Unknown test:", test, "\nBuilt-in codes: SW, SF, LF, KS, AD, AD2, CVM, JB, DAP, ANS, SKEW, KURT, MOM3, MOM4"))
  )
  
  return(result)
}


#' =========================================================================== #
#' 3. General function for computing area under curves for both normality
#'  tests and downstream test such as power and type I error curves
#'  --------------------------------------------------------------------------
compute_auc <- function(fpr, tpr, ensure_endpoints = TRUE, normalize = FALSE) {
  # Filter finite values
  valid_points <- is.finite(fpr) & is.finite(tpr)
  x <- fpr[valid_points]
  y <- tpr[valid_points]
  
  # Sort by FPR (and TPR for ties)
  sort_order <- order(x, y)
  fpr <- x[sort_order]
  tpr <- y[sort_order]
  
  if (length(fpr) < 2) return(NA)
  
  # Ensure (0,0) and (1,1) endpoints for ROC curves
  if (ensure_endpoints) {
    if (fpr[1] > 0 || tpr[1] > 0) { 
      fpr <- c(0, fpr) 
      tpr <- c(0, tpr) 
    }
    if (fpr[length(fpr)] < 1 || tpr[length(tpr)] < 1) { 
      fpr <- c(fpr, 1) 
      tpr <- c(tpr, 1) 
    }
  }
  
  # Trapezoidal integration
  auc <- sum(diff(fpr) * (head(tpr, -1) + tail(tpr, -1)) / 2)
  
  # Normalize for partial curves 
  if (normalize) {
    fpr_range <- max(fpr) - min(fpr)
    auc <- if (fpr_range > 0) auc / fpr_range else NA
  }
  
  return(auc)
}
#' ------------------------------------------------------------------------#

#' =============================================================================#
# 4. Define Utility Functions for Effect Size Specification and Processing
# -----------------------------------------------------------------------------#
zero_like <- function(x) {
  if (is.null(x)) return(0)
  if (!is.numeric(x)) stop("effect_size must be numeric (scalar or vector) or list(H0=..., H1=...).")
  z <- rep(0, length(x))
  names(z) <- names(x)
  z
}

#' -----------------------------------------------------------------------------#
split_effect_size <- function(effect_size) {
  # Case 1: User supplies explicit H0 and H1 as a named list
  if (is.list(effect_size) && !is.data.frame(effect_size)) {
    if (!all(c("H0","H1") %in% names(effect_size))) {
      stop("If effect_size is a list, it must be list(H0=..., H1=...).")
    }
    H0 <- effect_size$H0
    H1 <- effect_size$H1
    if (!is.numeric(H0) || !is.numeric(H1)) stop("effect_size$H0 and effect_size$H1 must be numeric.")
    if (length(H0) != length(H1)) stop("effect_size$H0 and effect_size$H1 must have the same length.")
    return(list(H0 = H0, H1 = H1))
  }
  
  # Case 2: User supplies only H1 values (H0 defaults to zero)
  if (!is.numeric(effect_size)) stop("effect_size must be numeric or list(H0=..., H1=...).")
  list(H0 = zero_like(effect_size), H1 = effect_size)
}
#' -----------------------------------------------------------------------------#


#' =====================================================================
#' 5. USER-DEFINED FUNCTIONS FOR DOWNSTREAM STATISTICAL TESTS
#' =====================================================================
#' ----------------------------------------------------------------------
#' A) ONE-SAMPLE TEST FUNCTIONS
#' ----------------------------------------------------------------------
#' (i) Generate Data for One-Sample Tests
onesample_data <- function(n = 20, effect_size = 0.0, sd = 1, dist = "Exponential", par = NULL, center_by = "median",...) {
  dist <- tolower(dist)
  x <- effect_size + sd * generate_data(n, dist, par = par, center_by = center_by)
  return(x)
}

#' -----------------------------------------------------------------
#' (ii) One-Sample Test Parameters. 
#' Returns list of parameters for one-sample tests
onesample_parameters <- function(n = 10, effect_size = 0.5, sd = 1, dist = "exponential", par = NULL, center_by = "median",...) {
  list(
    n = n,
    effect_size = effect_size,
    sd = sd,
    dist = dist,
    par = par,
    center_by = center_by 
  )
}

#' -----------------------------------------------------------------
#'(iii) Normality test object. 
#' Identity function returning input data unchanged
raw_data <- function(data) {
  return(data)
}

#' -----------------------------------------------------------------
#'(iv) One-Sample t-Test
one_sample_t_test <- function(data) {
  test_result <- t.test(data, mu = 0)
  return(list(p.value = test_result$p.value))
}

#' -----------------------------------------------------------------
#'(v) Sign Test (Binomial Test)
sign_test <- function(data, mu0 = 0) {
  x0 <- data - mu0
  n_valid <- sum(x0 != 0)
  if (n_valid == 0) return(list(p.value = 1))
  signs <- sum(x0 > 0)
  return(list(p.value = binom.test(signs, n_valid, p = 0.5)$p.value))
}

#' -----------------------------------------------------------------
#' (vi) One-Sample Permutation Test using t-test statistic
one_sample_perm_test <- function(data, mu0 = 0, nresample = 1000) {
  n <- length(data)
  observe_stat <- sqrt(n) * (mean(data) - mu0) / sd(data)
  permuted_stat <- replicate(nresample, {
    index <- sample(c(-1, 1), n, replace = TRUE)
    sample_data <- mu0 + index * (data - mu0)
    sqrt(n) * (mean(sample_data) - mu0) / sd(sample_data)
  })
  p <- (sum(abs(permuted_stat) >= abs(observe_stat)) + 1) / (nresample + 1)
  return(list(p.value = p))
}

#' -----------------------------------------------------------------
#'(vii) One-Sample Bootstrap-t Test using t-test statistic
one_sample_boot_test <- function(data, mu0 = 0, nresample = 1000) {
  data <- as.numeric(data)
  n <- length(data)
  t_obs <- sqrt(n) * (mean(data) - mu0) / sd(data)
  x0 <- data - mean(data) + mu0
  boot_stat <- replicate(nresample, {
    x_star <- sample(x0, size = n, replace = TRUE)
    sqrt(n) * (mean(x_star) - mu0) / sd(x_star)
  })
  p_value <- mean(abs(boot_stat) >= abs(t_obs), na.rm = TRUE)
  
  return(list(p.value = p_value))
}

#' -----------------------------------------------------------------
#'(viii) Bootstrap Median Test
bootstrap_median_test <- function(data, mu = 0, nresample = 1000) {
  data <- data[!is.na(data)]
  n <- length(data)
  obs_median <- median(data)
  centered_data <- data - obs_median + mu
  bootstrap_medians <- replicate(nresample, {
    median(sample(centered_data, n, replace = TRUE))
  })
  p_value <- mean(abs(bootstrap_medians - mu) >= abs(obs_median - mu))
  
  return(list(p.value = p_value))
}

#' ----------------------------------------------------------------------
#' B) TWO-SAMPLE TEST FUNCTIONS
#' ----------------------------------------------------------------------
#' (i) Generate Data for Two-Sample Tests
two_sample_data <- function(n1 = 10, n2 = 10, mean1 = 0.0, effect_size = 0.0, sd1 = 1, sd2 = 1, dist = "exponential", par = NULL, center_by = "median",...) {
  
  group1 <- mean1 + sd1 * generate_data(n = n1, dist, par = par, center_by = center_by)
  group2 <- (mean1 + effect_size) + sd2 * generate_data(n = n2, dist, par = par, center_by = center_by)
  
  return(data.frame(
    group = factor(rep(c("x", "y"), c(n1, n2))),
    value = c(group1, group2)
  ))
}
#' -----------------------------------------------------------------
#' (ii) Two-Sample Test Parameters
twosample_parameters <- function(n = 10, effect_size = 0.5, sd = 1, dist = "exponential", par = NULL, center_by = "median", ...) {
  
  #' Handle flexible n specification
  if (length(n) == 1) {
    n1 <- n
    n2 <- n
  } else if (length(n) == 2) {
    n1 <- n[1]
    n2 <- n[2]
  } else {
    stop("n should be length 1 or 2 for two-sample case")
  }
  
  list(
    n1 = n1,
    n2 = n2,
    mean1 = 0,  
    effect_size = effect_size,  
    sd1 = sd,
    sd2 = sd,  
    dist = dist,
    par = par,
    center_by = center_by
  )
}

#' -----------------------------------------------------------------
#' (iii) Two-Sample t-Test
twosample_t_test <- function(data) {
  x_data <- data$value[data$group == "x"]
  y_data <- data$value[data$group == "y"]
  test_result <- t.test(x = x_data, y = y_data)
  return(list(p.value = test_result$p.value))
}

#' -----------------------------------------------------------------
#' (iv) Mann-Whitney U Test (Wilcoxon Rank-Sum Test)
Mann_whitney_U_test <- function(data) {
  x_data <- data$value[data$group == "x"]
  y_data <- data$value[data$group == "y"]
  test_result <- wilcox.test(x_data, y_data)
  return(list(p.value = test_result$p.value))
}

#' -----------------------------------------------------------------
#'(v) Two-Sample Permutation Test
two_sample_perm_test <- function(data, nresample = 1e3) {
  stopifnot(all(c("group", "value") %in% names(data)))
  data <- data[!is.na(data$value) & !is.na(data$group), ]
  
  x <- data$value[data$group == "x"]
  y <- data$value[data$group == "y"]
  observed <- (mean(x) - mean(y))/sqrt(var(x)/length(x) + var(y)/length(y))
  
  perm_stats <- replicate(nresample, {
    g_star <- sample(data$group)
    x_star <- data$value[g_star == "x"]
    y_star <- data$value[g_star == "y"]
    (mean(x_star) - mean(y_star))/sqrt(var(x_star)/length(x_star) + var(y_star)/length(y_star))
  })
  
  return(list(p.value = mean(abs(perm_stats) >= abs(observed))))
}

#' ----------------------------------------------------------------------
#' C) ANOVA TEST FUNCTIONS
#' ----------------------------------------------------------------------
#' (i) Generate Data for One-Way ANOVA
#' 
#' Creates k independent groups with different means specified by effect_size.
#' Each element of effect_size is the mean for that group.
#' ----------------------------------------------------------------------------
anova_gen_data <- function(n = 10, effect_size = c(0, 0, 0), sd = 1, dist = "normal", par = NULL, center_by = "median", ...) {
  
  #' Validate inputs
  if (!is.numeric(effect_size) || length(effect_size) < 2) {
    stop("effect_size must be a numeric vector with at least 2 elements")
  }
  
  #' # Detect if all groups are identical 
  #' if (length(unique(effect_size)) == 1) {
  #'   warning("All effect sizes are identical - ANOVA F-test may not be appropriate")
  #' }
  
  k <- length(effect_size)
  group_labels <- LETTERS[1:k]
  
  #' Generate data for each group
  values <- unlist(lapply(seq_along(effect_size), function(i) {
    mean_i <- effect_size[i]
    mean_i + sd * generate_data(n, dist, par = par, center_by = center_by)
  }))
  
  #' Create data frame with proper factor ordering
  data.frame(
    group = factor(rep(group_labels, each = n), levels = group_labels),
    value = values
  )
}

#' -----------------------------------------------------------------
#' (ii) ANOVA Parameters
#' 
#' Create Parameter List for ANOVA Data Generation
#' -----------------------------------------------------------------
anova_parameters <- function(n = 10,effect_size = c(0, 0, 0), sd = 1, dist = "exponential", par = NULL, center_by = "median", ...) {
  list(
    n = n,
    effect_size = effect_size,
    sd = sd,
    dist = dist,
    par = par,
    center_by = center_by
  )
}
#' -----------------------------------------------------------------
#'(ii) Extract Residuals from ANOVA Model
#' Fits one-way ANOVA and returns residuals for normality testing.
#' ----------------------------------------------------------------------------
anova_residuals <- function(data) {
  #' Validate input
  if (!all(c("group", "value") %in% names(data))) {
    stop("data must contain 'group' and 'value' columns")
  }
  #' Fit ANOVA model
  model <- aov(value ~ group, data = data)
  #' Extract residuals
  return(residuals(model))
}

#' -----------------------------------------------------------------
#'(iii) Perform One-Way ANOVA. Tests whether group means differ.
#' -----------------------------------------------------------------
one_way_anova <- function(data) {
  # Validate input
  if (!all(c("group", "value") %in% names(data))) {
    stop("data must contain 'group' and 'value' columns")
  }
  #' Check for sufficient groups
  n_groups <- length(unique(data$group))
  if (n_groups < 2) {
    stop("Need at least 2 groups for ANOVA")
  }
  #' Perform ANOVA test
  aov_model <- aov(value ~ group, data = data)
  p_value <- summary(aov_model)[[1]]$"Pr(>F)"[1]
  return(list(p.value = p_value))
}

#' ---------------------------------------------------------------------------
#' (iv) Perform Kruskal-Wallis Test
#' 
#' Nonparametric alternative to one-way ANOVA. Tests whether
#' groups have different distributions (specifically, different medians).
#' ----------------------------------------------------------------------------
kruskal_wallis_test <- function(data) {
  #' Validate input
  if (!all(c("group", "value") %in% names(data))) {
    stop("data must contain 'group' and 'value' columns")
  }
  #' Check for sufficient groups
  n_groups <- length(unique(data$group))
  if (n_groups < 2) {
    stop("Need at least 2 groups for Kruskal-Wallis test")
  }
  #' Perform test
  test_result <- kruskal.test(value ~ group, data = data)
  return(list(p.value = test_result$p.value))
}

#' -----------------------------------------------------------------
#'(v) Permutation-Based ANOVA
#' 
#' Nonparametric alternative to ANOVA using permutation distribution.
#' Uses coin package implementation.
#' -------------------------------------------------------------------
permutation_anova <- function(data) {
  #' Validate input
  if (!all(c("group", "value") %in% names(data))) {
    stop("data must contain 'group' and 'value' columns")
  }
  #' Ensure group is a factor
  data$group <- as.factor(data$group)
  #' Perform permutation test
  test_result <- coin::oneway_test(value ~ group, data = data, distribution = approximate(nresample = 1e3)
  )
  
  return(list(p.value = coin::pvalue(test_result)))
}

#' ----------------------------------------------------------------------
#' D) LINEAR REGRESSION TEST FUNCTIONS
#' ----------------------------------------------------------------------
#' (i) Generate Data for Linear Regression
#' Generates predictor X and response Y from linear model Y = beta0 + beta1*X + error
#' 
reg_data <- function(n = 10, beta0 = 0, beta1 = 0, x_dist = "exponential", error_sd = 1, par = NULL, dist = "normal", center_by = "median", ...) {
  #' Generate predictor variable X from specified distribution
  x <- generate_data(n, dist = x_dist, par = par, center_by = center_by)
  #' Generate error term: 
  error <- error_sd * generate_data(n, dist = dist, par = par, center_by = center_by)
  #' Compute response variable using linear model: Y = beta0 + beta1*X + error
  y <- beta0 + beta1 * x + error
  return(data.frame(x = x, y = y))
}

#' -----------------------------------------------------------------
#' (ii) Regression Parameters
#' Creates parameter list for regression simulations
#' 
reg_parameters <- function(n = 10, effect_size = 0, beta0 = 0, beta1 = NULL,  x_dist = "exponential", error_sd = 1, dist = "normal", par = NULL, center_by = "median", ...) {
  
  #' If beta1 not specified, use effect_size instead
  if (is.null(beta1)) {
    beta1 <- effect_size
  }
  #' Return named list of all regression parameters
  list(
    n = n,
    beta0 = beta0,
    beta1 = beta1,
    x_dist = x_dist,
    error_sd = error_sd,
    dist = dist,
    par = par,
    center_by = center_by
  )
}

#' -----------------------------------------------------------------
#' (iii) Regression Residuals
#' Extracts residuals from simple linear regression model
#' 
reg_residuals <- function(data){
  #' Fit simple linear regression model
  model <- lm(y ~ x, data = data)
  #' Extract and return residuals 
  return(residuals(model))
}

#' -----------------------------------------------------------------
#' (iv) Simple Linear Regression
#' Performs simple linear regression and extracts p-value for slope
#' 
simple_linear_reg <- function(data) {
  #' Fit linear model Y ~ X
  model <- lm(y ~ x, data = data)
  #' Extract p-value for slope coefficient (tests H0: beta1 = 0)
  p_value <- summary(model)$coefficients["x", "Pr(>|t|)"]
  return(list(p.value = p_value))
}

#' -----------------------------------------------------------------
#' (v) Permutation Regression
#' Performs permutation test for regression slope significance
#' 
perm_regression <- function(data, nresample = 10000) {
  #' Remove non-finite values  to ensure valid analysis
  data <- data[is.finite(data$x) & is.finite(data$y), , drop = FALSE]
  #' Calculate observed slope coefficient from original data
  original_slope <- coef(lm(y ~ x, data = data))[["x"]]
  #' Extract x and y vectors 
  y <- data$y
  x <- data$x
  #' Generate permutation distribution of slope under null hypothesis
  #' (randomly shuffle y values while keeping x fixed)
  perm_slopes <- replicate(nresample, {
    y_perm <- sample(y, replace = FALSE)  # Permute y values
    coef(lm(y_perm ~ x))[["x"]]           # Calculate slope for permuted data
  })
  #' Calculate two-tailed p-value using absolute values
  #' Add 1 to numerator and denominator for continuity correction
  p_value <- (sum(abs(perm_slopes) >= abs(original_slope)) + 1) / (nresample + 1)
  return(list(p.value = p_value))
}

#' -----------------------------------------------------------------
#' (vi) Rank-Based Regression
#' Performs rank-based (robust) regression using Rfit package
#' Requires Rfit package to be installed
#' 
rank_regression <- function(data) {
  # Fit rank-based regression (robust to outliers and non-normality)
  result <- Rfit::rfit(y ~ x, data = data)
  # Extract p-value for slope coefficient 
  p_value <- summary(result)$coefficients[2, 4]
  return(list(p.value = p_value))
}


#' -----------------------------------------------------------------------------
#' 6.9. HELPER FUNCTION 1: CONVERT VARIED INPUT FORMATS TO STANDARD SAMPLE LIST
#' -----------------------------------------------------------------------------
convert_to_sample_list <- function(data) {
  
  # Identify column type for data frame parsing
  col_type <- function(x) {
    if (is.numeric(x)) return("numeric")
    if (is.character(x) || is.factor(x)) return("character")
    return("other")
  }
  
  #' --- CASE 1: Single numeric vector -----------------------------------------
  if (is.numeric(data) && !is.matrix(data)) {
    return(list(
      samples = list(Sample1 = as.numeric(data)),
      sample_names = "Sample1",
      n_samples = 1
    ))
  }
  
  #' --- CASE 2: List of numeric vectors ---------------------------------------
  if (is.list(data) && !is.data.frame(data) && all(sapply(data, is.numeric))) {
    samples <- lapply(data, as.numeric)
    sample_names <- names(data)
    if (is.null(sample_names)) sample_names <- paste0("Sample_", seq_along(samples))
    
    return(list(
      samples = samples,
      sample_names = sample_names,
      n_samples = length(samples)
    ))
  }
  
  #' --- CASE 3: Numeric matrix (each column = sample) -------------------------
  if (is.matrix(data) && is.numeric(data)) {
    samples <- lapply(seq_len(ncol(data)), function(i) as.numeric(data[, i]))
    sample_names <- colnames(data)
    if (is.null(sample_names)) sample_names <- paste0("Sample_", seq_along(samples))
    
    return(list(
      samples = samples,
      sample_names = sample_names,
      n_samples = length(samples)
    ))
  }
  
  #' --- CASE 4: Data frame (multiple possible layouts) -----------------------
  if (is.data.frame(data)) {
    types <- sapply(data, col_type)
    numeric_cols <- types == "numeric"
    char_cols <- types == "character"
    
    # Case 4A: All numeric -> each column is a sample (wide format)
    if (all(numeric_cols)) {
      samples <- lapply(data, as.numeric)
      sample_names <- colnames(data)
      if (is.null(sample_names)) sample_names <- paste0("Sample_", seq_along(samples))
      
      return(list(
        samples = samples,
        sample_names = sample_names,
        n_samples = length(samples)
      ))
    }
    
    # Case 4B: 1 char + 1 numeric -> group/value format (long format)
    if (sum(char_cols) == 1 && sum(numeric_cols) == 1) {
      group_col <- names(data)[char_cols][1]
      value_col <- names(data)[numeric_cols][1]
      
      samples <- split(data[[value_col]], data[[group_col]])
      samples <- lapply(samples, as.numeric)
      
      return(list(
        samples = samples,
        sample_names = names(samples),
        n_samples = length(samples)
      ))
    }
    
    # Case 4C: Mix of numeric/other -> use numeric columns as samples
    if (any(numeric_cols)) {
      samples <- lapply(data[numeric_cols], as.numeric)
      
      return(list(
        samples = samples,
        sample_names = names(data)[numeric_cols],
        n_samples = length(samples)
      ))
    }
  }
  
  #' --- CHECK ACCEPTABLE FORMAT -----------------------------------------------
  stop(paste(
    "sample_data format not recognized. Expected formats:\n",
    "1. Numeric vector\n",
    "2. List of numeric vectors\n",
    "3. Data frame with 1 character/factor column and 1 numeric column\n",
    "4. Data frame with numeric columns (each column = sample)\n",
    "5. Matrix (each column = sample)\n",
    "Received class:", paste(class(data), collapse = ", ")
  ))
}


#' -----------------------------------------------------------------------------
#' 6.10. HELPER FUNCTION 2: RETRIEVE APPROPRIATE ML MODEL BUNDLE
#' -----------------------------------------------------------------------------

get_ml_model <- function(trained_models,
                         actual_size,
                         requested_size = NULL) {
  
  target_size <- if (!is.null(requested_size)) {
    requested_size
  } else {
    actual_size
  }
  
  available_sizes <- as.numeric(names(trained_models))
  
  if (anyNA(available_sizes)) {
    stop("trained_models must be a named list with numeric sample-size names.")
  }
  
  closest_size <- available_sizes[
    which.min(abs(available_sizes - target_size))
  ]
  
  if (closest_size != target_size) {
    warning(
      sprintf(
        "No model trained for n = %d. Using closest available model: n = %d.",
        target_size,
        closest_size
      )
    )
  }
  
  model_bundle <- trained_models[[as.character(closest_size)]]
  
  if (is.null(model_bundle)) {
    stop("No valid ML model found for n = ", closest_size)
  }
  
  model_bundle
}


#' ============================================================
#' ML ADAPTER: CLASSIFY SAMPLES USING SOURCED ML FRAMEWORK
#' ============================================================
#' This function replaces the old classify_sample_prob().
#' It uses the new ML framework stored in ml_env.
#' ============================================================

classify_sample_prob <- function(sample_data,
                                 trained_models,
                                 model_name         = "RF",
                                 use_majority_vote  = TRUE,
                                 decision_threshold = 0.50,
                                 custom_sample_size = NULL) {
  
  if (is.null(trained_models)) stop("trained_models must be provided.")
  if (is.null(sample_data))    stop("sample_data cannot be NULL.")
  
  valid_models <- c("LR", "RF", "ANN", "GBM", "SVM", "KNN")
  
  invalid_models <- setdiff(model_name, valid_models)
  
  if (length(invalid_models) > 0L) {
    stop(
      "Invalid model(s): ",
      paste(invalid_models, collapse = ", "),
      ". Valid options: ",
      paste(valid_models, collapse = ", ")
    )
  }
  
  ## Convert input to list of numeric samples.
  sample_info  <- convert_to_sample_list(sample_data)
  sample_list  <- sample_info$samples
  sample_names <- sample_info$sample_names
  n_samples    <- sample_info$n_samples
  
  if (n_samples == 0L) {
    stop("No valid samples found in sample_data.")
  }
  
  ## Select the closest trained model bundle.
  actual_size <- length(sample_list[[1]])
  
  model_bundle <- get_ml_model(
    trained_models  = trained_models,
    actual_size     = actual_size,
    requested_size  = custom_sample_size
  )
  
  available_models <- names(model_bundle$models)
  models_available <- intersect(model_name, available_models)
  
  if (length(models_available) == 0L) {
    stop(
      "None of the requested models are available.\n",
      "Requested: ",
      paste(model_name, collapse = ", "),
      "\nAvailable: ",
      paste(available_models, collapse = ", ")
    )
  }
  
  if (length(models_available) < length(model_name)) {
    warning(
      "Some requested models are not available: ",
      paste(setdiff(model_name, models_available), collapse = ", "),
      "\nUsing only: ",
      paste(models_available, collapse = ", ")
    )
  }
  
  multiple_models <- length(models_available) > 1L
  
  results <- vector("list", n_samples)
  
  for (i in seq_len(n_samples)) {
    
    x_i         <- sample_list[[i]]
    sample_name <- sample_names[i]
    sample_size <- length(x_i)
    
    sample_result <- data.frame(
      Sample_Name = sample_name,
      Sample_Size = sample_size,
      stringsAsFactors = FALSE
    )
    
    ## Compute features and apply the exact preprocessing from the ML module.
    feats_std <- tryCatch(
      ml_env$prepare_prediction_features(
        x            = x_i,
        model_bundle = model_bundle
      ),
      error = function(e) {
        warning(
          "Feature preparation failed for sample '",
          sample_name,
          "': ",
          conditionMessage(e)
        )
        NULL
      }
    )
    
    if (is.null(feats_std)) {
      
      for (m in models_available) {
        sample_result[[paste0(m, "_Prob_Normal")]]     <- NA_real_
        sample_result[[paste0(m, "_Prob_Non_Normal")]] <- NA_real_
        sample_result[[paste0(m, "_Class")]]           <- NA_character_
      }
      
      if (multiple_models && use_majority_vote) {
        sample_result$Majority_Vote_Class   <- NA_character_
        sample_result$Proportion_Non_Normal <- NA_real_
      }
      
      results[[i]] <- sample_result
      next
    }
    
    prob_nn_vec <- numeric(length(models_available))
    names(prob_nn_vec) <- models_available
    
    for (m in models_available) {
      
      pred_probs <- predict(
        model_bundle$models[[m]],
        newdata = feats_std,
        type    = "prob"
      )
      
      prob_normal     <- as.numeric(pred_probs[, "Normal"])
      prob_non_normal <- as.numeric(pred_probs[, "Non_Normal"])
      
      pred_class <- if (prob_non_normal >= decision_threshold) {
        "Non_Normal"
      } else {
        "Normal"
      }
      
      sample_result[[paste0(m, "_Prob_Normal")]]     <- prob_normal
      sample_result[[paste0(m, "_Prob_Non_Normal")]] <- prob_non_normal
      sample_result[[paste0(m, "_Class")]]           <- pred_class
      
      prob_nn_vec[m] <- prob_non_normal
    }
    
    if (multiple_models && use_majority_vote) {
      
      mean_prob_nn <- mean(prob_nn_vec, na.rm = TRUE)
      
      sample_result$Majority_Vote_Class <- if (mean_prob_nn >= decision_threshold) {
        "Non_Normal"
      } else {
        "Normal"
      }
      
      sample_result$Proportion_Non_Normal <- mean_prob_nn
    }
    
    results[[i]] <- sample_result
  }
  
  out <- do.call(rbind, results)
  rownames(out) <- NULL
  
  attr(out, "settings") <- list(
    models_requested   = model_name,
    models_used        = models_available,
    decision_threshold = decision_threshold,
    use_majority_vote  = use_majority_vote,
    model_sample_size  = model_bundle$n,
    actual_sample_size = actual_size,
    timestamp          = Sys.time()
  )
  
  class(out) <- c("ml_normality_result", "data.frame")
  
  out
}

#' ----------------------------------------------------------------------------
#' 7.1 Normality test dispatcher — unified interface for multiple data formats
#'
#' Applies generate_tests() across however many samples are present in the
#' input and returns the decision-relevant quantity for each sample.
#' What that quantity is (p-value, moment, etc.) is entirely determined by
#' generate_tests() — this function does not need to know.
#'
#' Supported input formats:
#'   Case 1 — single numeric vector
#'   Case 2 — named or unnamed list of numeric vectors
#'   Case 3 — wide-format data frame (each numeric column = one sample)
#'   Case 4 — long-format data frame / matrix (group col + value col)
#' -----------------------------------------------------------------------------
normality_test <- function(data,
                           test      = "SW",
                           mu        = 0,
                           sigma     = 1,
                           group_col = 1,
                           value_col = 2) {
  
  
  apply_test <- function(x) {
    if (!is.numeric(x)) stop("All samples must be numeric vectors.")
    generate_tests(x, test = test, mu = mu, sigma = sigma)
  }
  
  #' ---------------------------------------------------------------------------
  #' Dispatch on input type — collect one decision_threshold per sample
  #' ---------------------------------------------------------------------------
  
  #' -- Case 1: single numeric vector ------------------------------------------
  if (is.numeric(data) && is.null(dim(data))) {
    decision_threshold <- apply_test(data)
    names(decision_threshold) <- "Sample1"
    
    #' -- Case 2: list of numeric vectors ----------------------------------------
  } else if (is.list(data) && !is.data.frame(data)) {
    decision_threshold <- sapply(data, apply_test)
    names(decision_threshold) <- if (!is.null(names(data))) {
      names(data)
    } else {
      paste0("Sample", seq_along(decision_threshold))
    }
    
    #' -- Case 3: wide-format data frame (each numeric column = one sample) ------
  } else if (is.data.frame(data) && all(sapply(data, is.numeric))) {
    decision_threshold <- sapply(as.list(data), apply_test)
    names(decision_threshold) <- names(data)
    
    #' -- Case 4: long-format data frame / matrix (group col + value col) --------
  } else if ((is.data.frame(data) || is.matrix(data)) && ncol(data) >= 2) {
    values <- data[[value_col]]
    if (!is.numeric(values)) {
      stop("value_col must point to a numeric column. ", "Check group_col and value_col arguments.")
    }
    grouped_samples    <- split(as.numeric(values), data[[group_col]])
    decision_threshold <- sapply(grouped_samples, apply_test)
    
  } else {
    stop(
      "Unsupported input format. 'data' must be one of:\n",
      "  - a numeric vector\n",
      "  - a list of numeric vectors\n",
      "  - a wide-format data frame (all numeric columns)\n",
      "  - a long-format data frame with a group column and a value column"
    )
  }
  
  return(decision_threshold)
}

#' -----------------------------------------------------------------------------
#' 7.2a: FUNCTION TO COMPUTE ROC CURVES FOR CUSTOM NORMALITY TEST FUNCTIONS
#'
#' This function handles custom normality scores, including:
#'   1. A single custom score, such as Fisher-combined p-value.
#'   2. Multiple ML scores, such as RF, GBM, ANN, and MajorityVote.
#'
#' Important convention:
#'   The custom function should return a numeric value where smaller values
#'   indicate stronger evidence against normality.
#'
#' For ML, ml_fn_roc() returns P(Normal), so smaller values mean less normal.
#'
#' For multiple groups, such as two-sample data, the function applies the custom
#' normality function to each group and combines the results using the minimum.
#' This means the whole dataset is treated as non-normal if any group appears
#' non-normal.
#' -----------------------------------------------------------------------------

fn_for_custom_test_roc_curve <- function(n                  = 10,
                                         threshold_grid     = seq(0, 1, by = 0.025),
                                         H1_dist            = "exponential",
                                         Nsim               = 1e3,
                                         gen_data           = two_sample_data,
                                         get_parameters     = twosample_parameters,
                                         fn_to_get_norm_obj = raw_data,
                                         custom_fn,
                                         custom_args        = list(),
                                         ...) {
  
  #' --------------------------------------------------------------------------
  #' STEP 1: VALIDATE INPUTS
  #' --------------------------------------------------------------------------
  if (missing(custom_fn) || !is.function(custom_fn)) {
    stop("custom_fn must be supplied and must be a function.")
  }
  
  #' --------------------------------------------------------------------------
  #' STEP 2: PROBE CUSTOM FUNCTION TO DETERMINE OUTPUT LENGTH
  #' --------------------------------------------------------------------------
  #' The probe tells us whether custom_fn returns:
  #'   - one value, e.g. Fisher p-value;
  #'   - several values, e.g. one ML score per model.
  #' --------------------------------------------------------------------------
  
  probe_paras <- get_parameters(
    n    = n,
    dist = "Normal",
    par  = NULL,
    ...
  )
  
  probe_data    <- do.call(gen_data, probe_paras)
  probe_norm    <- fn_to_get_norm_obj(probe_data)
  probe_samples <- convert_to_sample_list(probe_norm)$samples
  
  probe_out <- do.call(
    custom_fn,
    c(list(probe_samples[[1]]), custom_args)
  )
  
  if (!is.numeric(probe_out)) {
    stop("custom_fn must return a numeric vector.")
  }
  
  n_curves <- length(probe_out)
  
  curve_names <- if (!is.null(names(probe_out))) {
    names(probe_out)
  } else {
    paste0("custom_", seq_len(n_curves))
  }
  
  #' --------------------------------------------------------------------------
  #' STEP 3: HELPER TO EVALUATE ONE GENERATED DATASET
  #' --------------------------------------------------------------------------
  #' This helper is called for each H0 and H1 simulation replicate.
  #'
  #' For one-sample data:
  #'   samples has length 1, so this simply returns custom_fn(sample).
  #'
  #' For two-sample data:
  #'   samples has length 2, so custom_fn is applied to each group.
  #'   We then take the minimum score per model.
  #'
  #' For ML:
  #'   custom_fn returns several values, one per model.
  #'   The output is always returned as a numeric vector of length n_curves.
  #' --------------------------------------------------------------------------
  
  evaluate_custom_score <- function(data_obj) {
    
    samples <- convert_to_sample_list(data_obj)$samples
    
    vals <- lapply(samples, function(s) {
      as.numeric(do.call(custom_fn, c(list(s), custom_args)))
    })
    
    vals_mat <- do.call(cbind, vals)
    
    if (is.null(dim(vals_mat))) {
      vals_mat <- matrix(vals_mat, nrow = n_curves)
    }
    
    ## Combine multiple samples/groups by taking the minimum score.
    ## This is appropriate because smaller values mean less normal.
    out <- apply(vals_mat, 1, min, na.rm = TRUE)
    
    if (length(out) != n_curves) {
      stop("custom_fn returned inconsistent output length.")
    }
    
    names(out) <- curve_names
    out
  }
  
  #' --------------------------------------------------------------------------
  #' STEP 4: STORAGE
  #' --------------------------------------------------------------------------
  cat("Computing ROC curve(s) for custom normality function\n")
  
  out_H0 <- matrix(
    NA_real_,
    nrow     = n_curves,
    ncol     = Nsim,
    dimnames = list(curve_names, NULL)
  )
  
  out_H1 <- matrix(
    NA_real_,
    nrow     = n_curves,
    ncol     = Nsim,
    dimnames = list(curve_names, NULL)
  )
  
  pb      <- txtProgressBar(min = 0, max = 2 * Nsim, style = 3)
  counter <- 0L
  
  #' --------------------------------------------------------------------------
  #' STEP 5: SIMULATE CUSTOM SCORES UNDER H0
  #' --------------------------------------------------------------------------
  for (i in seq_len(Nsim)) {
    
    paras_H0 <- get_parameters(
      n    = n,
      dist = "Normal",
      par  = NULL,
      ...
    )
    
    data_H0 <- do.call(gen_data, paras_H0)
    norm_H0 <- fn_to_get_norm_obj(data_H0)
    
    out_H0[, i] <- tryCatch(
      evaluate_custom_score(norm_H0),
      error = function(e) rep(NA_real_, n_curves)
    )
    
    counter <- counter + 1L
    setTxtProgressBar(pb, counter)
  }
  
  #' --------------------------------------------------------------------------
  #' STEP 6: SIMULATE CUSTOM SCORES UNDER H1
  #' --------------------------------------------------------------------------
  for (i in seq_len(Nsim)) {
    
    paras_H1 <- get_parameters(
      n    = n,
      dist = H1_dist,
      ...
    )
    
    data_H1 <- do.call(gen_data, paras_H1)
    norm_H1 <- fn_to_get_norm_obj(data_H1)
    
    out_H1[, i] <- tryCatch(
      evaluate_custom_score(norm_H1),
      error = function(e) rep(NA_real_, n_curves)
    )
    
    counter <- counter + 1L
    setTxtProgressBar(pb, counter)
  }
  
  close(pb)
  
  #' --------------------------------------------------------------------------
  #' STEP 7: COMPUTE ROC COORDINATES
  #' --------------------------------------------------------------------------
  #' Since smaller scores mean stronger evidence against normality:
  #'   reject normality when score <= threshold.
  #' --------------------------------------------------------------------------
  
  FPR_mat <- matrix(
    NA_real_,
    nrow     = n_curves,
    ncol     = length(threshold_grid),
    dimnames = list(curve_names, paste0("t_", threshold_grid))
  )
  
  TPR_mat <- matrix(
    NA_real_,
    nrow     = n_curves,
    ncol     = length(threshold_grid),
    dimnames = list(curve_names, paste0("t_", threshold_grid))
  )
  
  for (i in seq_len(n_curves)) {
    FPR_mat[i, ] <- sapply(
      threshold_grid,
      function(t) mean(out_H0[i, ] <= t, na.rm = TRUE)
    )
    
    TPR_mat[i, ] <- sapply(
      threshold_grid,
      function(t) mean(out_H1[i, ] <= t, na.rm = TRUE)
    )
  }
  
  #' --------------------------------------------------------------------------
  #' STEP 8: RETURN
  #' --------------------------------------------------------------------------
  list(
    FPR       = FPR_mat,
    TPR       = TPR_mat,
    threshold = threshold_grid,
    raw       = list(
      H0 = out_H0,
      H1 = out_H1
    )
  )
}

#' -----------------------------------------------------------------------------
#' 7.2b: FUNCTION TO COMPUTE ROC CURVES FOR CLASSICAL NORMALITY TESTS
#' -----------------------------------------------------------------------------
#' -----------------------------------------------------------------------------
fn_for_norm_test_roc_curve <- function(n = 10, 
                                       alpha_pretest = seq(from = 0, to = 1, by = 0.25), 
                                       H1_dist = "exponential", 
                                       tests = c("SW", "SF", "LF", "JB"), 
                                       Nsim = 1e3, 
                                       gen_data = two_sample_data,
                                       get_parameters = twosample_parameters,
                                       fn_to_get_norm_obj = raw_data,
                                       ...) {
  
  #' --------------------------------------------------------------------------
  #' STEP 1: INITIALIZE RESULT MATRICES
  #' --------------------------------------------------------------------------
  #' Rows = different normality tests, Columns = different alpha thresholds
  FPR <- matrix(0, nrow = length(tests), ncol = length(alpha_pretest))
  TPR <- matrix(0, nrow = length(tests), ncol = length(alpha_pretest))
  
  #' Set row and column names for readability
  rownames(FPR) <- rownames(TPR) <- tests
  colnames(FPR) <- colnames(TPR) <- paste0("alpha_", alpha_pretest)
  
  #' --------------------------------------------------------------------------
  #' STEP 2: SETUP PROGRESS BAR
  #' --------------------------------------------------------------------------
  cat("Computing ROC curves for normality tests\n")
  total_iterations <- Nsim * length(tests) * length(alpha_pretest)
  pb <- txtProgressBar(min = 0, max = total_iterations, style = 3)
  counter <- 0
  
  #' --------------------------------------------------------------------------
  #' STEP 3: MAIN SIMULATION LOOP
  #'   - Loop over each normality test
  #'   - Loop over each alpha threshold
  #'   - Run Nsim simulations for each combination
  #' --------------------------------------------------------------------------
  for (i in seq_along(tests)) {
    test_name <- tests[i]
    
    for (j in seq_along(alpha_pretest)) {
      alpha <- alpha_pretest[j]
      
      #' Storage for rejection decisions across Nsim iterations
      reject_H0 <- numeric(Nsim)  
      reject_H1 <- numeric(Nsim)  
      
      for (k in 1:Nsim) {
        
        #' ----- H0: Generate data from Normal distribution -----
        paras_H0 <- get_parameters(n, dist = "Normal", par = NULL, ...)
        normal_data <- do.call(gen_data, paras_H0)
        normal_obj <- fn_to_get_norm_obj(normal_data)
        pvals_H0 <- normality_test(normal_obj, test = test_name)
        #' Reject if ANY p-value < alpha (for multiple groups/tests)
        reject_H0[k] <- any(pvals_H0 < alpha, na.rm = TRUE)
        
        #' ----- H1: Generate data from alternative distribution -----
        paras_H1 <- get_parameters(n, dist = H1_dist, ...)
        non_normal_data <- do.call(gen_data, paras_H1)
        non_normal_obj <- fn_to_get_norm_obj(non_normal_data)
        pvals_H1 <- normality_test(non_normal_obj, test = test_name)
        #' Reject if ANY p-value < alpha
        reject_H1[k] <- any(pvals_H1 < alpha, na.rm = TRUE)
        
        #' Update progress bar
        counter <- counter + 1
        setTxtProgressBar(pb, counter)
      }
      
      #' ----------------------------------------------------------------------
      #' STEP 4: COMPUTE PERFORMANCE METRICS
      #' ----------------------------------------------------------------------
      #' FPR (False Positive Rate) 
      FPR[i, j] <- mean(reject_H0, na.rm = TRUE)
      
      #' TPR (True Positive Rate) 
      TPR[i, j] <- mean(reject_H1, na.rm = TRUE)
    }
  }
  
  #' Close progress bar
  close(pb)
  
  #' --------------------------------------------------------------------------
  #' STEP 5: RETURN RESULTS
  #' --------------------------------------------------------------------------
  return(list(
    FPR = FPR,
    TPR = TPR,
    alpha = alpha_pretest
  ))
}



#' ----------------------------------------------------------------------------
#' 7.4 assess_normality(): unified normality assessment wrapper
#'
#' A single entry point for normality assessment with two branches:
#'
#'   Branch 1 — "classical": delegates to normality_test(), which calls
#'              generate_tests() with a built-in test code (e.g. "SW", "AD").
#'
#'   Branch 2 — "custom": applies any user-supplied function fn() to each
#'              sample. This branch handles BOTH the ML approach (by wrapping
#'              classify_sample_prob() inside a closure before passing it in)
#'              and any other user-defined normality test.
#'
#' ----------------------------------------------------------------------------
assess_normality <- function(data, method = "classical", config = list()) {
  
  #' -------------------------------------------------------------------------
  #' 0. Validate method
  #' -------------------------------------------------------------------------
  method <- tolower(trimws(method))
  if (!method %in% c("classical", "custom")) {
    stop("method must be 'classical' or 'custom'.")
  }
  
  #' -------------------------------------------------------------------------
  #' BRANCH 1: Classical — delegate entirely to normality_test()
  #' -------------------------------------------------------------------------
  if (method == "classical") {
    
    norm_test <- config$norm_test %||% "SW"
    mu        <- config$mu        %||% 0
    sigma     <- config$sigma     %||% 1
    group_col <- config$group_col %||% 1
    value_col <- config$value_col %||% 2
    
    return(normality_test(data,
                          test      = norm_test,
                          mu        = mu,
                          sigma     = sigma,
                          group_col = group_col,
                          value_col = value_col))
  }
  
  #' -------------------------------------------------------------------------
  #' BRANCH 2: Custom — apply user-supplied function fn() to each sample.
  #'
  #' This branch is intentionally general: fn() can be
  #'   (a) any user-defined test that returns a single named numeric scalar, OR
  #'   (b) a closure that wraps classify_sample_prob() for the ML approach.
  #'
  #' The caller is responsible for binding any extra arguments (e.g.
  #' trained_models, decision_threshold) inside the closure before passing fn.
  #' -------------------------------------------------------------------------
  if (method == "custom") {
    
    if (is.null(config$fn) || !is.function(config$fn)) {
      stop("config$fn must be a function for method = 'custom'.")
    }
    
    fn         <- config$fn
    extra_args <- config[setdiff(names(config), "fn")]
    
    ## Convert data to a standard list of numeric vectors
    sample_info  <- convert_to_sample_list(data)
    sample_list  <- sample_info$samples
    sample_names <- sample_info$sample_names
    
    decision <- sapply(sample_list, function(x) {
      out <- do.call(fn, c(list(x), extra_args))
      if (!is.numeric(out) || length(out) != 1) {
        stop("config$fn must return a single named numeric scalar. ",
             "Got: ", paste(class(out), collapse = ", "),
             " of length ", length(out))
      }
      out
    })
    
    names(decision) <- sample_names
    return(decision)
  }
}

#' ========================================================
#' Unified plotting function for ROC curves
#' Improved legibility, original line types preserved
#'========================================================

plot_norm_roc_curve <- function(curves,
                                title = NULL,
                                dist_name = NULL) {
  
  if (is.null(title)) {
    suffix <- if (!is.null(dist_name)) paste(" | Alternative:", dist_name) else ""
    title  <- paste("ROC Curves for Normality Tests", suffix)
  }
  
  n_curves <- length(curves)
  
  # Better color palette
  colors <- 1:n_curves
  
  # Point shapes
  pchs <- seq(16, length.out = n_curves)
  
  # Save and restore graphics settings
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  # Improve spacing and text size
  par(mar = c(5, 5, 4, 2) + 0.1,
      cex.lab = 1.2,
      cex.axis = 1.1,
      cex.main = 1.4)
  
  # Empty plot
  plot(0, 0, type = "n",
       xlim = c(0, 1), ylim = c(0, 1),
       xlab = "False Positive Rate (FPR)",
       ylab = "True Positive Rate (TPR)",
       main = title)
  
  # Light grid for readability
  grid(col = "lightgray", lty = "dotted")
  
  # Reference diagonal
  abline(0, 1, lty = 2, col = "gray60", lwd = 2)
  
  auc_values <- numeric(n_curves)
  
  for (i in seq_along(curves)) {
    cv  <- curves[[i]]
    ord <- order(cv$fpr)
    fpr <- cv$fpr[ord]
    tpr <- cv$tpr[ord]
    
    # KEEP original line types from cv$lty
    lines(fpr, tpr,
          col = colors[i],
          lwd = 3,
          lty = cv$lty)
    
    # Points (visible but not cluttered)
    points(fpr, tpr,
           col = colors[i],
           pch = pchs[i],
           cex = 0.7)
    
    auc_values[i] <- compute_auc(fpr, tpr, ensure_endpoints = TRUE)
  }
  
  legend_text <- sprintf("%s (AUC = %.3f)",
                         sapply(curves, `[[`, "label"),
                         auc_values)
  
  legend("bottomright",
         legend = legend_text,
         col    = colors,
         pch    = pchs,
         lty    = sapply(curves, `[[`, "lty"),  # preserve original
         lwd    = 3,
         pt.cex = 0.9,
         title  = "Normality Tests",
         cex    = 0.9,
         bty    = "b")
  
  invisible(auc_values)
}


#' ---------------------------------------------------------------------------
#' 7.0 Dispatches to the appropriate ROC simulation function based on the
#' method in norm_config, then passes the result to plot_norm_roc_curve().
#'
#'   "classical" : calls fn_for_norm_test_roc_curve() to simulate p-values
#'                 across the alpha_pretest grid for all selected classical
#'                 tests, then plots all curves together.
#'
#'   "custom"    : calls fn_for_custom_test_roc_curve() using the function
#'                 and threshold_grid supplied in norm_config$config. Plots
#'                 one ROC curve per value returned by custom_fn. If custom_fn
#'                 wraps an ML classifier, all ML-specific logic is handled
#'                 inside the closure — this branch has no knowledge of those
#'                 details.
#' ---------------------------------------------------------------------------
plot_normality_roc_wrapper <- function(norm_config,
                                       single_n,
                                       threshold_grid,
                                       distributions,
                                       Nsim               = 1e3,
                                       norm_test          = c("SW","SF","LF","KS",
                                                              "JB","SKEW","DAP",
                                                              "AD","CVM"),
                                       selected_tests     = c("SW","SF","LF","JB",
                                                              "SKEW","DAP","AD","CVM"),
                                       gen_data,
                                       get_parameters,
                                       fn_to_get_norm_obj,
                                       center_by          = "median",
                                       ...) {
  
  method <- tolower(norm_config$method %||% "classical")
  if (!method %in% c("classical", "custom"))
    stop("norm_config$method must be 'classical' or 'custom'.")
  
  #' --------------------------------------------------------------------------
  #' BRANCH 1: Classical
  #' --------------------------------------------------------------------------
  if (method == "classical") {
    
    norm_test_method <- norm_config$config$norm_test %||% "SW"
    tests_to_run     <- if (norm_test_method %in% norm_test) norm_test else norm_test_method
    
    roc_obj <- fn_for_norm_test_roc_curve(
      n                  = single_n,
      alpha_pretest      = threshold_grid,
      H1_dist            = distributions[1],
      tests              = tests_to_run,
      Nsim               = Nsim,
      center_by          = center_by,
      gen_data           = gen_data,
      get_parameters     = get_parameters,
      fn_to_get_norm_obj = fn_to_get_norm_obj,
      ...
    )
    
    computed_tests <- rownames(roc_obj$FPR)
    tests_to_plot  <- intersect(selected_tests, computed_tests)
    if (length(tests_to_plot) == 0) tests_to_plot <- computed_tests
    
    curves <- lapply(tests_to_plot, function(t) {
      list(fpr   = as.numeric(roc_obj$FPR[t, ]),
           tpr   = as.numeric(roc_obj$TPR[t, ]),
           label = t,
           lty   = 1)
    })
    plot_norm_roc_curve(curves    = curves,
                        dist_name = distributions[1])
    
    return(invisible(roc_obj))
  }
  
  #' --------------------------------------------------------------------------
  #' BRANCH 2: Custom
  #' --------------------------------------------------------------------------
  if (method == "custom") {
    
    if (is.null(norm_config$config$fn) || !is.function(norm_config$config$fn))
      stop("norm_config$config$fn must be a function for method = 'custom'.")
    
    fn          <- norm_config$config$fn
    custom_args <- norm_config$config[setdiff(names(norm_config$config), "fn")]
    curve_label <- norm_config$config$label %||% norm_config$config$name %||% "Custom"
    
    roc_obj <- fn_for_custom_test_roc_curve(
      n                  = single_n,
      threshold_grid     = threshold_grid,
      H1_dist            = distributions[1],
      Nsim               = Nsim,
      gen_data           = gen_data,
      get_parameters     = get_parameters,
      fn_to_get_norm_obj = fn_to_get_norm_obj,
      custom_fn          = fn,
      custom_args        = custom_args,
      center_by          = center_by,
      ...
    )
    
    #' For single-curve custom tests (e.g. Fisher), apply user label.
    #' For multi-curve custom tests (e.g. ML), rownames are already set
    #' by fn_for_custom_test_roc_curve from names(probe_one) — keep them.
    if (nrow(roc_obj$FPR) == 1) {
      label <- norm_config$config$label %||% norm_config$config$name %||% "Custom"
      rownames(roc_obj$FPR) <- rownames(roc_obj$TPR) <- label
    }
    
    curve_names <- rownames(roc_obj$FPR)
    curves <- lapply(seq_along(curve_names), function(i) {
      list(fpr   = as.numeric(roc_obj$FPR[i, ]),
           tpr   = as.numeric(roc_obj$TPR[i, ]),
           label = curve_names[i],
           lty   = 1)
    })
    plot_norm_roc_curve(curves    = curves,
                        dist_name = distributions[1])
    
    return(invisible(roc_obj))
  }
}

#' --------------------------------------------------------------------------
#' Internal parallel helper (FORK-based, macOS/Linux friendly)
#'
#' Maps `fun` over seq_len(n_iter) in parallel and returns an ordered list.
#' Uses doMC (forking via mclapply): workers inherit the master process's
#' already-attached packages and defined functions, so NO .packages / .export
#' is required. Reproducibility is handled by doRNG (%dorng%) and is invariant
#' to the number of cores. A text progress bar is advanced in the master after
#' each chunk completes, giving live, meaningful feedback. On Windows (no fork)
#' it falls back to sequential execution so nothing needs to be exported.
#' --------------------------------------------------------------------------
ds_parallel_map <- function(n_iter, fun, label = NULL, n_cores = NULL,
                            max_updates = 100L) {
  if (is.null(n_cores)) {
    n_cores <- getOption("ds.parallel.cores", max(1L, parallel::detectCores() - 1L))
  }
  
  ## Register a FORK backend on Unix-like systems; sequential elsewhere.
  if (.Platform$OS.type != "windows" && requireNamespace("doMC", quietly = TRUE)) {
    doMC::registerDoMC(n_cores)
  } else {
    foreach::registerDoSEQ()
  }
  
  if (!is.null(label)) cat(label, "\n")
  
  ## Split iterations into chunks so the progress bar updates incrementally.
  n_chunks <- max(1L, min(n_iter, as.integer(max_updates)))
  chunks   <- split(seq_len(n_iter), cut(seq_len(n_iter), n_chunks, labels = FALSE))
  
  pb <- txtProgressBar(min = 0, max = n_iter, style = 3)
  on.exit(close(pb), add = TRUE)
  
  out  <- vector("list", n_iter)
  done <- 0L
  for (ch in chunks) {
    res     <- foreach::foreach(i = ch) %dorng% fun(i)
    out[ch] <- res
    done    <- done + length(ch)
    setTxtProgressBar(pb, done)
  }
  out
}

#' --------------------------------------------------------------------------
#' 8.0 Simulate normality pretest decision values and downstream 
#' test p-values under H0 and H1 across Nsim replicates
#' --------------------------------------------------------------------------

generate_pval <- function(Nsim,
                          n,
                          effect_size_H1,
                          effect_size_H0     = zero_like(effect_size_H1),
                          norm_config        = list(method = "classical",config = list(norm_test = "SW")),
                          dist               = "normal",
                          gen_data,
                          get_parameters,
                          fn_to_get_norm_obj,
                          fn_for_ds_test_1,
                          fn_for_ds_test_2,
                          ...) {
  
  #' ---------------------------------------------------------------------------
  #' Validate norm_config
  #' ---------------------------------------------------------------------------
  if (is.null(norm_config$method)) {
    stop("norm_config must have a $method slot: 'classical', or 'custom'.")
  }
  if (!tolower(norm_config$method) %in% c("classical", "custom")) {
    stop("norm_config$method must be 'classical', or 'custom'.")
  }
  
  #' ---------------------------------------------------------------------------
  #' Pre-allocate storage
  #' ---------------------------------------------------------------------------
  pval_ds_test1_H0 <- rep(NA_real_, Nsim)
  pval_ds_test1_H1 <- rep(NA_real_, Nsim)
  pval_ds_test2_H0 <- rep(NA_real_, Nsim)
  pval_ds_test2_H1 <- rep(NA_real_, Nsim)
  
  ## norm_decision stores the named numeric vector from assess_normality()
  ## per replicate — same structure regardless of which method was chosen
  norm_decision_H0 <- vector("list", Nsim)
  norm_decision_H1 <- vector("list", Nsim)
  
  #' ---------------------------------------------------------------------------
  #' Simulation loop (parallelised over Nsim replicates)
  #' ---------------------------------------------------------------------------
  ## Capture dots once so they can be forwarded safely inside the parallel body
  dots <- list(...)
  
  ## FORK-based parallel map (workers inherit packages/functions; no exports).
  ## Progress bar ticks as replicate chunks complete.
  sim_out <- ds_parallel_map(Nsim, function(i) {
    
    ## Generate data under H0 and H1
    paras_H0 <- do.call(get_parameters, c(list(n, effect_size = effect_size_H0, dist = dist), dots))
    data_H0  <- do.call(gen_data, paras_H0)
    
    paras_H1 <- do.call(get_parameters, c(list(n, effect_size = effect_size_H1, dist = dist), dots))
    data_H1  <- do.call(gen_data, paras_H1)
    
    ## Extract the object to be assessed for normality
    normality_obj_H0 <- fn_to_get_norm_obj(data_H0)
    normality_obj_H1 <- fn_to_get_norm_obj(data_H1)
    
    ## Downstream p-values + normality decisions (named numeric vectors)
    list(
      p1_H0 = fn_for_ds_test_1(data_H0)$p.value,
      p1_H1 = fn_for_ds_test_1(data_H1)$p.value,
      p2_H0 = fn_for_ds_test_2(data_H0)$p.value,
      p2_H1 = fn_for_ds_test_2(data_H1)$p.value,
      nd_H0 = assess_normality(normality_obj_H0,
                               method = norm_config$method,
                               config = norm_config$config),
      nd_H1 = assess_normality(normality_obj_H1,
                               method = norm_config$method,
                               config = norm_config$config)
    )
  })
  
  ## Reassemble flat result vectors / lists in replicate order
  pval_ds_test1_H0 <- vapply(sim_out, `[[`, numeric(1), "p1_H0")
  pval_ds_test1_H1 <- vapply(sim_out, `[[`, numeric(1), "p1_H1")
  pval_ds_test2_H0 <- vapply(sim_out, `[[`, numeric(1), "p2_H0")
  pval_ds_test2_H1 <- vapply(sim_out, `[[`, numeric(1), "p2_H1")
  norm_decision_H0 <- lapply(sim_out, `[[`, "nd_H0")
  norm_decision_H1 <- lapply(sim_out, `[[`, "nd_H1")
  
  #' ---------------------------------------------------------------------------
  #' Return
  #' ---------------------------------------------------------------------------
  list(
    pval_ds_test1_H0 = pval_ds_test1_H0,
    pval_ds_test2_H0 = pval_ds_test2_H0,
    pval_ds_test1_H1 = pval_ds_test1_H1,
    pval_ds_test2_H1 = pval_ds_test2_H1,
    norm_decision_H0 = norm_decision_H0,
    norm_decision_H1 = norm_decision_H1
  )
}



#' ----------------------------------------------------------------------------
#' 8.1 Perform adaptive statistical analysis with normality pretest
#'
#' Sweeps over a single threshold_grid and computes Type I error rates and
#' power for the adaptive strategy and each individual downstream test.
#'
#' The decision rule is uniform across all methods:
#'   use parametric test (test_1) when assess_normality() output > threshold
#'   fall back to non-parametric test (test_2) otherwise
#' -----------------------------------------------------------------------------

perform_adaptive_analysis <- function(Nsim               = 1e3,
                                      n                  = 10,
                                      effect_size_H1     = 0.5,
                                      effect_size_H0     = zero_like(effect_size_H1),
                                      distributions      = c("exponential", "normal"),
                                      norm_config        = list(method = "classical", config = list(norm_test = "SW")),
                                      threshold_grid     = seq(0, 1, by = 0.025),
                                      test_alpha         = 0.05,
                                      gen_data           = onesample_data,
                                      get_parameters     = onesample_parameters,
                                      fn_to_get_norm_obj = raw_data,
                                      fn_for_ds_test_1   = one_sample_t_test,
                                      fn_for_ds_test_2   = sign_test,
                                      ...) {
  
  #' ---------------------------------------------------------------------------
  #' Validate inputs
  #' ---------------------------------------------------------------------------
  if (is.null(norm_config$method)) {
    stop("norm_config must have a $method slot: 'classical', 'ML', or 'custom'.")
  }
  if (length(distributions) < 2) {
    stop("distributions must have at least 2 elements: ", "[1] = non-normal (H1), [2] = 'normal' (H0).")
  }
  if (tolower(distributions[2]) != "normal") {
    stop("distributions[2] must be 'normal' (the H0 reference distribution).", "Got: '", distributions[2], "'.")
  }
  
  #' ---------------------------------------------------------------------------
  #' Sweep over distributions
  #' ---------------------------------------------------------------------------
  all_pvalues <- list()
  typeI_rates <- list()
  power_rates <- list()
  
  pb_dist <- txtProgressBar(min = 0, max = length(distributions), style = 3)
  
  for (dist_idx in seq_along(distributions)) {
    
    dist <- distributions[dist_idx]
    
    sim_results <- generate_pval(
      Nsim               = Nsim,
      n                  = n,
      effect_size_H1     = effect_size_H1,
      effect_size_H0     = effect_size_H0,
      norm_config        = norm_config,
      dist               = dist,
      gen_data           = gen_data,
      get_parameters     = get_parameters,
      fn_to_get_norm_obj = fn_to_get_norm_obj,
      fn_for_ds_test_1   = fn_for_ds_test_1,
      fn_for_ds_test_2   = fn_for_ds_test_2,
      ...
    )
    all_pvalues[[dist]] <- sim_results
    
    pval_H0_test1 <- sim_results$pval_ds_test1_H0
    pval_H0_test2 <- sim_results$pval_ds_test2_H0
    pval_H1_test1 <- sim_results$pval_ds_test1_H1
    pval_H1_test2 <- sim_results$pval_ds_test2_H1
    
    ## Marginal rates for each individual test
    typeI_rates[[dist]] <- list(
      test_1   = mean(pval_H0_test1 < test_alpha, na.rm = TRUE),
      test_2   = mean(pval_H0_test2 < test_alpha, na.rm = TRUE),
      adaptive = numeric(length(threshold_grid))
    )
    power_rates[[dist]] <- list(
      test_1   = mean(pval_H1_test1 < test_alpha, na.rm = TRUE),
      test_2   = mean(pval_H1_test2 < test_alpha, na.rm = TRUE),
      adaptive = numeric(length(threshold_grid))
    )
    
    ## adaptive test
    for (j in seq_along(threshold_grid)) {
      threshold <- threshold_grid[j]
      
      use_par_H0 <- sapply(sim_results$norm_decision_H0,
                           function(d) all(d > threshold, na.rm = TRUE))
      use_par_H1 <- sapply(sim_results$norm_decision_H1,
                           function(d) all(d > threshold, na.rm = TRUE))
      
      adaptive_H0 <- ifelse(use_par_H0, pval_H0_test1, pval_H0_test2)
      adaptive_H1 <- ifelse(use_par_H1, pval_H1_test1, pval_H1_test2)
      
      typeI_rates[[dist]]$adaptive[j] <- mean(adaptive_H0 < test_alpha, na.rm = TRUE)
      power_rates[[dist]]$adaptive[j] <- mean(adaptive_H1 < test_alpha, na.rm = TRUE)
    }
    
    setTxtProgressBar(pb_dist, dist_idx)
  }
  close(pb_dist)
  
  #' ---------------------------------------------------------------------------
  #' Return
  #' ---------------------------------------------------------------------------
  list(
    typeI_rates = typeI_rates,
    power_rates = power_rates,
    all_pvalues = all_pvalues,
    param_grid  = threshold_grid,
    param_name  = "threshold",
    settings    = list(
      Nsim               = Nsim,
      n                  = n,
      effect_size_H1     = effect_size_H1,
      effect_size_H0     = effect_size_H0,
      distributions      = distributions,
      norm_config        = norm_config,
      threshold_grid     = threshold_grid,
      test_alpha         = test_alpha,
      gen_data           = deparse(substitute(gen_data)),
      get_parameters     = deparse(substitute(get_parameters)),
      fn_to_get_norm_obj = deparse(substitute(fn_to_get_norm_obj)),
      fn_for_ds_test_1   = deparse(substitute(fn_for_ds_test_1)),
      fn_for_ds_test_2   = deparse(substitute(fn_for_ds_test_2))
    )
  )
}


#' =============================================================================
#' 8.2. Compute ROC-style metrics: EPL, EPG, and ETIE
#' Evaluates the adaptive pretest strategy by comparing it against each
#' individual downstream test across the threshold parameter grid.
#' =============================================================================

compute_roc_metrics <- function(typeI_rates, power_rates, test_alpha) {
  #' ---------------------------------------------------------------------------
  #' Validate distribution ordering
  #' ---------------------------------------------------------------------------
  nms <- names(power_rates)
  if (!is.null(nms) && length(nms) >= 2) {
    if (tolower(nms[2]) != "normal") {
      stop("Distributions must be ordered: (non-normal, normal). ", "Second element is '", nms[2], "'.")
    }
  } else {
    warning("Cannot verify distribution order: power_rates is not a named list.")
  }
  
  #' ---------------------------------------------------------------------------
  #' Extract performance under each distributional scenario
  #' ---------------------------------------------------------------------------
  
  ## Non-normal scenario (H1): assesses power of the adaptive strategy
  non_normal <- list(
    power_test1    = power_rates[[1]]$test_1,
    power_test2    = power_rates[[1]]$test_2,
    power_adaptive = power_rates[[1]]$adaptive,
    error_adaptive = typeI_rates[[1]]$adaptive
  )
  
  ## Normal scenario (H0): assesses Type I error and power loss of the adaptive strategy
  normal <- list(
    power_test1    = power_rates[[2]]$test_1,
    power_test2    = power_rates[[2]]$test_2,
    power_adaptive = power_rates[[2]]$adaptive,
    error_adaptive = typeI_rates[[2]]$adaptive
  )
  
  #' ---------------------------------------------------------------------------
  #' Compute evaluation metrics
  #' ---------------------------------------------------------------------------
  
  ## Expected Power Loss (EPL): power sacrificed by the adaptive strategy
  ## relative to test_1 when data are actually normal
  EPL <- normal$power_test1 - normal$power_adaptive
  
  ## Expected Power Gain (EPG): power gained by the adaptive strategy
  ## relative to test_1 when data are non-normal
  EPG <- non_normal$power_adaptive - non_normal$power_test1
  
  ## Expected Type I Error Inflation (ETIE): how much the adaptive strategy
  ## inflates the Type I error above the nominal level, under each scenario
  ETIE_normal     <- normal$error_adaptive     - test_alpha
  ETIE_non_normal <- non_normal$error_adaptive - test_alpha
  
  ## Benchmark: direct comparison of test_2 vs test_1 (no adaptive strategy)
  benchmark_power_gain <- non_normal$power_test2 - non_normal$power_test1
  benchmark_power_loss <- normal$power_test1     - normal$power_test2
  
  #' ---------------------------------------------------------------------------
  #' Return
  #' ---------------------------------------------------------------------------
  list(
    Expected_power_loss                       = EPL,
    Expected_power_gain                       = EPG,
    Expected_type1_error_inflation_normal     = ETIE_normal,
    Expected_type1_error_inflation_non_normal = ETIE_non_normal,
    benchmark_power_gain                      = benchmark_power_gain,
    benchmark_power_loss                      = benchmark_power_loss
  )
}


#' =============================================================================#
#' 8.3 Select the optimal pretest threshold from the trade-off metrics
#' Implements a two-stage selection logic:
#' Stage 1 (preferred): among all thresholds where worst-case positive Type I
#'  inflation <= tol_pos, maximize power gain subject to power loss <= loss_tol
#'  Stage 2 (fallback):  if no threshold satisfies Stage 1, minimize inflation
#'  first, then apply the same power criterion within near-optimal candidates
#'   Stage 3 (extreme fallback): all thresholds are conservative (both
#'  distributions show negative inflation) — select on power criterion alone
#' =============================================================================#

select_optimal_parameter <- function(distributions = c("exponential", "normal"),
                                     param_grid,
                                     metrics,
                                     alpha      = 0.05,
                                     tol_pos    = 0.01,
                                     loss_tol   = 0.01,
                                     param_name = "threshold") {
  
  #' --------------------------------------------------------------------------
  #' Extract and coerce metric vectors
  #' --------------------------------------------------------------------------
  power_loss_normal      <- as.numeric(metrics$Expected_power_loss)
  power_gain_non_normal  <- as.numeric(metrics$Expected_power_gain)
  type1_infl_normal      <- as.numeric(metrics$Expected_type1_error_inflation_normal)
  type1_infl_non_normal  <- as.numeric(metrics$Expected_type1_error_inflation_non_normal)
  
  #' --------------------------------------------------------------------------
  #' Length alignment: param_grid and metrics vectors must match
  #' --------------------------------------------------------------------------
  metrics_length <- length(power_loss_normal)
  grid_length    <- length(param_grid)
  
  if (metrics_length != grid_length) {
    min_length <- min(metrics_length, grid_length)
    warning(sprintf(
      "param_grid length (%d) does not match metrics length (%d). ",
      "Truncating both to %d elements.",
      grid_length, metrics_length, min_length
    ))
    param_grid            <- param_grid[seq_len(min_length)]
    power_loss_normal     <- power_loss_normal[seq_len(min_length)]
    power_gain_non_normal <- power_gain_non_normal[seq_len(min_length)]
    type1_infl_normal     <- type1_infl_normal[seq_len(min_length)]
    type1_infl_non_normal <- type1_infl_non_normal[seq_len(min_length)]
  }
  
  #' --------------------------------------------------------------------------#
  #' Filter out non-finite values before any comparisons
  #' --------------------------------------------------------------------------#
  valid <- is.finite(param_grid) &
    is.finite(power_loss_normal) &
    is.finite(power_gain_non_normal) &
    is.finite(type1_infl_normal) &
    is.finite(type1_infl_non_normal)
  
  param_values          <- param_grid[valid]
  power_loss            <- power_loss_normal[valid]
  power_gain            <- power_gain_non_normal[valid]
  type1_infl_normal     <- type1_infl_normal[valid]
  type1_infl_non_normal <- type1_infl_non_normal[valid]
  
  #' --------------------------------------------------------------------------#
  #' Inflation classification
  #' --------------------------------------------------------------------------#
  eps <- .Machine$double.eps^0.5   # floating-point comparison tolerance
  
  ## Positive inflation: undesirable, needs to be controlled
  pos_infl_normal     <- pmax(type1_infl_normal,     0)
  pos_infl_non_normal <- pmax(type1_infl_non_normal, 0)
  worst_pos_infl      <- pmax(pos_infl_normal, pos_infl_non_normal)
  
  ## Conservative points: both distributions show negative inflation
  both_negative     <- type1_infl_normal < 0 & type1_infl_non_normal < 0
  not_both_negative <- !both_negative
  
  ## --------------------------------------------------------------------------#
  ## Inner helper: given a set of candidate indices, select the one that
  ## maximizes power gain subject to power loss <= loss_tol; ties broken by
  ## smallest parameter value.
  ## --------------------------------------------------------------------------#
  select_best <- function(candidates) {
    
    loss_vals <- power_loss[candidates]
    gain_vals <- power_gain[candidates]
    
    ## Prefer candidates where power loss is within tolerance
    good_loss <- candidates[loss_vals <= loss_tol + eps]
    
    if (length(good_loss) > 0) {
      ## Among low-loss candidates, maximize power gain
      max_gain <- max(power_gain[good_loss], na.rm = TRUE)
      best     <- good_loss[abs(power_gain[good_loss] - max_gain) <= eps * (1 + abs(max_gain))]
    } else {
      ## No candidate meets loss tolerance: minimize loss first
      min_loss  <- min(loss_vals, na.rm = TRUE)
      min_group <- candidates[abs(loss_vals - min_loss) <= eps * (1 + abs(min_loss))]
      max_gain  <- max(power_gain[min_group], na.rm = TRUE)
      best      <- min_group[abs(power_gain[min_group] - max_gain) <= eps * (1 + abs(max_gain))]
    }
    
    ## Tie-break: smallest parameter value (most conservative threshold)
    if (length(best) > 1) best[which.min(param_values[best])] else best
  }
  
  
  ## --------------------------------------------------------------------------#
  ## Ordered selection stages
  ## --------------------------------------------------------------------------#
  
  ## Candidates satisfying the Type I error inflation tolerance
  controlled_idx <- which(worst_pos_infl <= tol_pos + eps)
  
  ## Candidates that are not conservative under both scenarios
  nonconservative_idx <- which(not_both_negative)
  
  ## Stage 1: controlled and non-conservative
  preferred_idx <- intersect(controlled_idx, nonconservative_idx)
  
  if (length(preferred_idx) > 0) {
    
    optimal_index  <- select_best(preferred_idx)
    selection_note <- sprintf(
      "Stage 1: inflation controlled (<= %.3g) and not conservative under both scenarios",
      tol_pos
    )
    error_controlled <- TRUE
    feasible         <- TRUE
    
  } else if (length(controlled_idx) > 0) {
    
    ## Stage 2: error control is achievable only through conservative candidates
    optimal_index  <- select_best(controlled_idx)
    selection_note <- sprintf(
      "Stage 2: inflation controlled (<= %.3g), but all controlled candidates are conservative under both scenarios",
      tol_pos
    )
    error_controlled <- TRUE
    feasible         <- FALSE
    
  } else if (length(nonconservative_idx) > 0) {
    
    ## Stage 3: no candidate controls inflation; minimize inflation first
    min_infl <- min(worst_pos_infl[nonconservative_idx], na.rm = TRUE)
    
    near_min <- nonconservative_idx[
      abs(worst_pos_infl[nonconservative_idx] - min_infl) <=
        eps * (1 + abs(min_infl))
    ]
    
    optimal_index  <- select_best(near_min)
    selection_note <- sprintf(
      "Stage 3: no candidate controls inflation; minimized worst-case inflation (min = %.3g)",
      min_infl
    )
    error_controlled <- FALSE
    feasible         <- FALSE
    
  } else {
    
    ## Stage 4: defensive fallback
    min_infl <- min(worst_pos_infl, na.rm = TRUE)
    
    near_min <- which(
      abs(worst_pos_infl - min_infl) <= eps * (1 + abs(min_infl))
    )
    
    optimal_index  <- select_best(near_min)
    selection_note <- "Stage 4: defensive fallback over the full valid grid"
    error_controlled <- min_infl <= tol_pos + eps
    feasible         <- FALSE
  }
  
  
  ## --------------------------------------------------------------------------#
  ## Return
  ## --------------------------------------------------------------------------#
  list(
    param_star            = param_values[optimal_index],
    optimal_index         = optimal_index,
    power_loss            = power_loss[optimal_index],
    power_gain            = power_gain[optimal_index],
    inflation_normal      = type1_infl_normal[optimal_index],
    inflation_non_normal  = type1_infl_non_normal[optimal_index],
    selection_note        = selection_note,
    feasible              = feasible,
    param_name            = param_name,
    all_metrics = list(
      param_values                  = param_values,
      power_loss                    = power_loss,
      power_gain                    = power_gain,
      type1_infl_normal             = type1_infl_normal,
      type1_infl_non_normal         = type1_infl_non_normal,
      pos_infl_normal               = pos_infl_normal,
      pos_infl_non_normal           = pos_infl_non_normal,
      worst_case_positive_inflation = worst_pos_infl,
      distributions                 = distributions,
      alpha                         = alpha,
      tol_pos                       = tol_pos,
      loss_tol                      = loss_tol,
      param_name                    = param_name
    )
  )
}


## ============================================================================#
## 8.4 Plot the four trade-off curves with the selected optimal threshold
## Produces a 2×2 panel plot:
##   Panel 1 (top-left)  : Expected Power Loss under Normal (H0)
##   Panel 2 (top-right) : Expected Power Gain under non-normal (H1)
##   Panel 3 (bottom-left) : Type I Error Inflation under Normal
##   Panel 4 (bottom-right): Type I Error Inflation under non-normal
## ============================================================================#

plot_tradeoff_results_generic <- function(optimal_result,
                                          outer_title      = "Power and Type I Error Trade-off for Normality Pretesting",
                                          text_size_main   = 1.0,
                                          text_size_labels = 1.0,
                                          point_size       = 1.2) {
  
  ## --------------------------------------------------------------------------#
  ## Unpack inputs
  ## --------------------------------------------------------------------------#
  metrics               <- optimal_result$all_metrics
  param_values          <- metrics$param_values
  param_name            <- metrics$param_name
  power_loss            <- metrics$power_loss
  power_gain            <- metrics$power_gain
  type1_infl_normal     <- metrics$type1_infl_normal
  type1_infl_non_normal <- metrics$type1_infl_non_normal
  distributions         <- metrics$distributions
  tol_pos               <- metrics$tol_pos
  loss_tol              <- metrics$loss_tol
  
  optimal_index        <- optimal_result$optimal_index
  opt_param            <- optimal_result$param_star
  optimal_power_loss   <- optimal_result$power_loss
  optimal_power_gain   <- optimal_result$power_gain
  optimal_infl_normal  <- optimal_result$inflation_normal
  optimal_infl_non_normal <- optimal_result$inflation_non_normal
  selection_note       <- optimal_result$selection_note
  
  ## --------------------------------------------------------------------------#
  ## x-axis label: human-readable expression matching param_name
  ## --------------------------------------------------------------------------#
  xlab_expr <- switch(param_name,
                      threshold          = "Pretest Threshold",
                      pretest_threshold  = expression(alpha[pre]),
                      decision_threshold = "Decision Threshold",
                      param_name
  )
  
  ## --------------------------------------------------------------------------#
  ## Helper: add a highlighted optimal point clamped to the plot region
  ## Clamping uses 0.1% of the axis range to stay visibly inside the frame
  ## --------------------------------------------------------------------------#
  add_optimal_point <- function(x_val, y_val, color = "darkred", size = 1) {
    usr   <- par("usr")
    x_pad <- 0.001 * (usr[2] - usr[1])
    y_pad <- 0.001 * (usr[4] - usr[3])
    cx    <- min(max(x_val, usr[1] + x_pad), usr[2] - x_pad)
    cy    <- min(max(y_val, usr[3] + y_pad), usr[4] - y_pad)
    points(cx, cy, pch = 19, cex = size, col = color)
    text(cx, cy, labels = sprintf("threshold = %.4f", x_val), pos = 3, cex = 0.7, col = color, font = 2)
  }
  
  ## --------------------------------------------------------------------------#
  ## Core plotting function
  ## hline: optional numeric value for a horizontal reference line
  ## --------------------------------------------------------------------------#
  create_panel <- function(y, color, ylab, main, hline = NULL) {
    
    ## Proportional y-axis padding: 2% of data range
    y_range <- range(y, na.rm = TRUE)
    pad     <- 0.02 * max(diff(y_range), 1e-6)
    ylim    <- y_range + c(-pad, pad)
    
    plot(param_values, y,
         type     = "l",
         col      = color,
         lwd      = 3,
         ylim     = ylim,
         ylab     = ylab,
         xlab     = xlab_expr,
         cex.lab  = text_size_labels,
         main     = main,
         cex.main = text_size_main,
         font.lab = 2)
    
    ## Vertical line at selected threshold
    abline(v = param_values[optimal_index], lty = 2, col = "darkred", lwd = 1.5)
    
    ## Optional horizontal reference line (e.g. zero for inflation panels)
    if (!is.null(hline)) {
      abline(h = hline, lty = 3, col = "grey60")
    }
    
    ## Light scatter of all points
    #points(param_values, y, pch = 21, cex = 0.6, col = adjustcolor(color, alpha.f = 0.3))
    
    ## Highlighted optimal point with label
    #add_optimal_point(param_values[optimal_index], y[optimal_index], color = "darkred", size = point_size)
  }
  
  ## --------------------------------------------------------------------------#
  ## Layout and plot
  ## --------------------------------------------------------------------------#
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  
  par(mfrow = c(2, 2),
      oma  = c(0.8, 0.2, 3.5, 0.2),
      mar  = c(4, 3.5, 1.5, 0.8),
      mgp  = c(2.0, 0.6, 0))
  
  ## Panel 1: Expected Power Loss — Normal (H0)
  create_panel(power_loss, "blue",
               ylab  = "Expected Power Loss",
               main  = sprintf("Expected Power Loss (EPL): %s", distributions[2]))
  
  ## Panel 2: Expected Power Gain — Non-normal (H1)
  create_panel(power_gain, "red",
               ylab  = "Expected Power Gain",
               main  = sprintf("Expected Power Gain (EPG): %s", distributions[1]))
  
  ## Panel 3: Type I Error Inflation — Normal (H0), zero reference line
  create_panel(type1_infl_normal, "orange",
               ylab  = "Type I Error Inflation",
               main  = sprintf("Type I Error Inflation: %s", distributions[2]),
               hline = 0)
  
  ## Panel 4: Type I Error Inflation — Non-normal (H1), zero reference line
  create_panel(type1_infl_non_normal, "green4",
               ylab  = "Type I Error Inflation",
               main  = sprintf("Type I Error Inflation: %s", distributions[1]),
               hline = 0)
  
  ## Outer title
  mtext(outer_title, side = 3, outer = TRUE, line = 1.2, cex = 1.5, font = 2)
  
  ## Summary line beneath the title
  param_label <- param_name
  mtext(
    sprintf(
      "Selected threshold = %.4f | EPG = %.4f | EPL = %.4f | Inflation = (%.4f, %.4f) | tol_pos = %.4f | loss_tol = %.4f",
      opt_param,
      optimal_power_gain, optimal_power_loss,
      optimal_infl_normal, optimal_infl_non_normal,
      tol_pos, loss_tol
    ),
    side = 3, outer = TRUE, line = 0.0, cex = 0.95, col = "darkblue"
  )
  
  invisible(optimal_result)
}

#' =============================================================================
#' 8.5 Generate ROC data for Power vs Type I Error analysis
#'
#' Runs generate_pval() for each distribution at a fixed pretest threshold then
#' sweeps over sig_levels to build power and Type I error rate curves for
#' test_1, test_2, and the adaptive strategy.
#' ==============================================================================
power_vs_error_roc_data <- function(N                  = 1e3,
                                    n                  = 10,
                                    distributions      = c("exponential", "normal"),
                                    norm_config        = list(method = "classical", config = list(norm_test = "SW")),
                                    effect_size_H1,
                                    effect_size_H0     = zero_like(effect_size_H1),
                                    pretest_threshold  = 0.05,
                                    optimal_result     = NULL,
                                    sig_levels         = seq(0, 1, by = 0.05),
                                    sim_cache          = NULL,
                                    gen_data           = anova_gen_data,
                                    get_parameters     = anova_parameters,
                                    fn_to_get_norm_obj = anova_residuals,
                                    fn_for_ds_test_1   = one_way_anova,
                                    fn_for_ds_test_2   = kruskal_wallis_test,
                                    ...) {
  
  #' ---------------------------------------------------------------------------
  #' Validate norm_config
  #' ---------------------------------------------------------------------------
  if (is.null(norm_config$method)) {
    stop("norm_config must have a $method slot: 'classical', 'ML', or 'custom'.")
  }
  
  #' ---------------------------------------------------------------------------
  #' Override pretest_threshold from optimal_result if supplied
  #' param_name is always "threshold" in the current pipeline
  #' ---------------------------------------------------------------------------
  if (!is.null(optimal_result) && !is.null(optimal_result$param_star)) {
    pretest_threshold <- optimal_result$param_star
  }
  
  #' ---------------------------------------------------------------------------
  #' Storage and progress bar
  #' ---------------------------------------------------------------------------
  roc_results <- data.frame()
  
  total_ops <- length(distributions) * length(sig_levels)
  pb        <- txtProgressBar(0, total_ops, style = 3)
  on.exit(close(pb), add = TRUE)
  counter   <- 0L
  
  #' ---------------------------------------------------------------------------
  #' Cache reuse setup
  #'
  #' sim_cache, when supplied, is a per-distribution list of generate_pval()
  #' outputs (exactly the shape of perform_adaptive_analysis()$all_pvalues).
  #' If a distribution's entry is present and complete, it is reused verbatim
  #' instead of re-simulating. generate_pval() output depends only on
  #' (dist, n, effect sizes, norm_config, gen fns) — NOT on the pretest
  #' threshold or alpha — so reuse here is exact.
  #' ---------------------------------------------------------------------------
  reused_any <- FALSE
  needed     <- c("pval_ds_test1_H0", "pval_ds_test1_H1",
                  "pval_ds_test2_H0", "pval_ds_test2_H1",
                  "norm_decision_H0", "norm_decision_H1")
  
  #' ---------------------------------------------------------------------------
  #' Main loop: one call to generate_pval() per distribution, then sweep
  #' ---------------------------------------------------------------------------
  for (dist in distributions) {
    
    cat("Processing distribution:", dist, "\n")
    
    #' Reuse cached generate_pval() output when available and complete
    cached <- if (!is.null(sim_cache)) sim_cache[[dist]] else NULL
    if (!is.null(cached) && !all(needed %in% names(cached))) {
      warning("sim_cache for '", dist, "' is incomplete; re-simulating.")
      cached <- NULL
    }
    
    if (!is.null(cached)) {
      reused_any  <- TRUE
      message(sprintf("  Reusing cached generate_pval() for '%s' (%d replicates).",
                      dist, length(cached$pval_ds_test1_H0)))
      sim_results <- cached
    } else {
      sim_results <- generate_pval(
        Nsim               = N,
        n                  = n,
        effect_size_H1     = effect_size_H1,
        effect_size_H0     = effect_size_H0,
        norm_config        = norm_config,
        dist               = dist,
        gen_data           = gen_data,
        get_parameters     = get_parameters,
        fn_to_get_norm_obj = fn_to_get_norm_obj,
        fn_for_ds_test_1   = fn_for_ds_test_1,
        fn_for_ds_test_2   = fn_for_ds_test_2,
        ...
      )
    }
    
    #' Determine which replicates use the parametric test (test_1)
    #' Decision rule is uniform across all methods: norm_decision > threshold
    use_test1_H0 <- sapply(sim_results$norm_decision_H0,
                           function(d) all(d > pretest_threshold, na.rm = TRUE))
    use_test1_H1 <- sapply(sim_results$norm_decision_H1,
                           function(d) all(d > pretest_threshold, na.rm = TRUE))
    
    #' Adaptive p-values: test_1 when normality supported, test_2 otherwise
    adaptive_pvals_H0 <- ifelse(use_test1_H0,
                                sim_results$pval_ds_test1_H0,
                                sim_results$pval_ds_test2_H0)
    adaptive_pvals_H1 <- ifelse(use_test1_H1,
                                sim_results$pval_ds_test1_H1,
                                sim_results$pval_ds_test2_H1)
    
    #' Sweep over significance levels
    for (alpha in sig_levels) {
      
      roc_results <- rbind(
        roc_results,
        data.frame(Distribution = dist, Method = "test_1", Alpha = alpha,
                   Power       = mean(sim_results$pval_ds_test1_H1 <= alpha, na.rm = TRUE),
                   TypeI_error = mean(sim_results$pval_ds_test1_H0 <= alpha, na.rm = TRUE),
                   stringsAsFactors = FALSE),
        data.frame(Distribution = dist, Method = "test_2", Alpha = alpha,
                   Power       = mean(sim_results$pval_ds_test2_H1 <= alpha, na.rm = TRUE),
                   TypeI_error = mean(sim_results$pval_ds_test2_H0 <= alpha, na.rm = TRUE),
                   stringsAsFactors = FALSE),
        data.frame(Distribution = dist, Method = "adaptive", Alpha = alpha,
                   Power       = mean(adaptive_pvals_H1 <= alpha, na.rm = TRUE),
                   TypeI_error = mean(adaptive_pvals_H0 <= alpha, na.rm = TRUE),
                   stringsAsFactors = FALSE)
      )
      
      counter <- counter + 1L
      setTxtProgressBar(pb, counter)
    }
  }
  
  #' ---------------------------------------------------------------------------
  #' Attach settings for provenance and downstream use by the plot function
  #' ---------------------------------------------------------------------------
  #' effective_N is the actual replicate count used. When reusing a cache built
  #' at a larger Nsim (e.g. Phase 2's N_tradeoff), this exceeds the requested N.
  effective_N <- length(sim_results$pval_ds_test1_H0)
  
  attr(roc_results, "settings") <- list(
    N                 = effective_N,
    N_requested       = N,
    reused_from_cache = reused_any,
    n                 = n,
    distributions     = distributions,
    norm_config       = norm_config,
    pretest_threshold = pretest_threshold,
    sig_levels        = sig_levels,
    effect_size_H1    = effect_size_H1,
    effect_size_H0    = effect_size_H0,
    timestamp         = Sys.time()
  )
  
  roc_results
}


#' =============================================================================
#' 8.6 Plot Power vs Type I Error ROC-like curves
#' 
#' Produces a layout with:
#'   Row 1: Full view   (Type I error 0–1)
#'   Row 2: Zoomed view (Type I error 0–zoom_xlim)
#'   Row 3: Legend with performance coordinates at nominal_alpha
#' =============================================================================

power_vs_error_roc_plot <- function(roc_results,
                                    nominal_alpha  = 0.05,
                                    optimal_result = NULL,
                                    zoom_xlim      = 0.10) {
  
  #' ---------------------------------------------------------------------------
  #' Extract settings stored by power_vs_error_roc_data()
  #' ---------------------------------------------------------------------------
  settings <- attr(roc_results, "settings")
  if (is.null(settings)) {
    warning("roc_results is missing a settings attribute. Using defaults.")
    settings <- list(norm_config       = list(method = "classical"),
                     pretest_threshold = 0.05)
  }
  
  ## The threshold actually used when generating the data
  pretest_threshold <- settings$pretest_threshold
  
  ## If an optimal_result is supplied, its param_star should already 
  #' match but we read it anyway for the subtitle display
  if (!is.null(optimal_result) && !is.null(optimal_result$param_star)) {
    pretest_threshold <- optimal_result$param_star
  }
  
  ## Method label for subtitle — taken from norm_config$method if available
  method_label <- if (!is.null(settings$norm_config$method)) {
    switch(tolower(settings$norm_config$method),
           classical = "Classical",
           custom    = "Custom",
           settings$norm_config$method)
  } else "Unknown"
  
  #' ---------------------------------------------------------------------------
  #' Layout: 3 rows × n_dist columns
  #' ---------------------------------------------------------------------------
  distributions <- unique(roc_results$Distribution)
  n_dist        <- length(distributions)
  
  layout_matrix <- matrix(seq_len(3 * n_dist), nrow = 3, ncol = n_dist, byrow = TRUE)
  
  layout(layout_matrix, heights = c(1, 1, 0.28))
  
  par(mar = c(4, 4, 3.0, 1), 
      oma = c(1.2, 0, 4.2, 0),
      mgp = c(2.2, 0.7, 0),
      cex.axis = 1.0,
      cex.lab  = 1.2,
      font.lab  = 2,
      font.axis = 2
  )
  
  #' ---------------------------------------------------------------------------
  #' Colours and display names
  #' ---------------------------------------------------------------------------
  test_color   <- c(test_1 = "red", test_2 = "blue", adaptive = "green")
  method_names <- c(test_1 = "Test 1", test_2 = "Test 2", adaptive = "Adaptive")
  
  ## Store operating-point coordinates (at nominal_alpha) for the legend row
  legend_coords <- list()
  
  #' ---------------------------------------------------------------------------
  #' Helper: plot one panel (shared by full and zoomed rows)
  #' ---------------------------------------------------------------------------
  plot_roc_panel <- function(dist_data, dist, xlim, main, store_coords = FALSE) {
    plot(NA, xlim = xlim, ylim = c(0, 1), xlab = "P(Type I Error)", ylab = "Power", main = main,
         cex.axis = 1.1, cex.lab = 1.2)
    abline(0, 1,              col = "gray80", lty = 2)
    abline(v = nominal_alpha, col = "red",    lty = 2)
    
    coords <- list()
    for (method in names(test_color)) {
      if (!method %in% dist_data$Method) next
      md <- dist_data[dist_data$Method == method, ]
      md <- md[order(md$Alpha), ]
      lines(md$TypeI_error, md$Power, col = test_color[method], lwd = 3)
      
      idx <- which.min(abs(md$Alpha - nominal_alpha))
      points(md$TypeI_error[idx], md$Power[idx], col = test_color[method], pch = 16, cex = 2.0)
      
      if (store_coords) coords[[method]] <- c(md$TypeI_error[idx], md$Power[idx])
    }
    if (store_coords) coords else invisible(NULL)
  }
  
  #' ---------------------------------------------------------------------------
  #' Row 1: Full view (Type I error 0–1)
  #' ---------------------------------------------------------------------------
  for (dist in distributions) {
    dist_data <- roc_results[roc_results$Distribution == dist, ]
    coords    <- plot_roc_panel(dist_data, dist, xlim = c(0, 1), main = paste(dist, "- Full View"), store_coords = TRUE)
    legend_coords[[dist]] <- coords
  }
  
  #' ---------------------------------------------------------------------------
  #' Row 2: Zoomed view (Type I error 0–zoom_xlim)
  #' ---------------------------------------------------------------------------
  for (dist in distributions) {
    dist_data <- roc_results[roc_results$Distribution == dist, ]
    plot_roc_panel(dist_data, dist, xlim = c(0, zoom_xlim), main = paste(dist, "- Zoomed View"))
  }
  
  #' ---------------------------------------------------------------------------
  #' Row 3: Legend panels — one per distribution
  #' ---------------------------------------------------------------------------
  par(mar = c(0.2, 0.2, 0.2, 0.2))
  for (dist in distributions) {
    plot.new()
    
    legend_text   <- character(0)
    legend_colors <- character(0)
    
    for (method in names(test_color)) {
      if (!method %in% names(legend_coords[[dist]])) next
      coords <- legend_coords[[dist]][[method]]
      legend_text <- c(
        legend_text,
        sprintf("%s (alpha = %.3f): TIE = %.3f, Power = %.3f", method_names[method], nominal_alpha, coords[1], coords[2])
      )
      legend_colors <- c(legend_colors, test_color[method])
    }
    
    legend("center",
           legend    = legend_text,
           col       = legend_colors,
           lwd       = 4,
           pch       = 16,
           pt.cex    = 1.6,
           bty       = "b",
           cex       = 1.15,
           text.font = 2,
           title     = paste(dist, "Distribution"),
           title.adj = 0.5,
           title.cex = 1.2)
  }
  
  #' ---------------------------------------------------------------------------
  #' Outer titles
  #' ---------------------------------------------------------------------------
  mtext(
    sprintf("ROC-Like Curves: Power vs Type I Error Rate (%s Approach)", method_label),
    side = 3, outer = TRUE, cex = 1.3, line = 2.4, font = 2
  )
  
  mtext(
    sprintf("Nominal alpha = %.3f | Pretest threshold = %.4f",
            nominal_alpha, pretest_threshold),
    side = 3, outer = TRUE, cex = 1.0, line = 0.65
  )
  
  invisible(roc_results)
}

#' =============================================================================
#' 9.0 Apply the adaptive downstream test to a single dataset
#'
#' Runs assess_normality() on the normality object, compares the result to
#' pretest_threshold, and dispatches to fn_for_ds_test_1 (parametric) or
#' fn_for_ds_test_2 (non-parametric) accordingly.
#' =============================================================================

run_adaptive_test_dual <- function(data,
                                   norm_config = list(method = "classical", config = list(norm_test = "SW")),
                                   pretest_threshold,
                                   fn_to_get_norm_obj = raw_data,
                                   fn_for_ds_test_1   = twosample_t_test,
                                   fn_for_ds_test_2   = Mann_whitney_U_test,
                                   ...) {
  
  if (missing(pretest_threshold) || is.null(pretest_threshold)) {
    stop("pretest_threshold must be provided.")
  }
  if (is.null(norm_config$method)) {
    stop("norm_config must have a $method slot: 'classical' or 'custom'.")
  }
  
  ## Extract the object to be assessed for normality
  norm_obj <- fn_to_get_norm_obj(data)
  
  ## Assess normality — returns a named numeric vector regardless of method
  norm_decision <- assess_normality(norm_obj,
                                    method = norm_config$method,
                                    config = norm_config$config)
  
  ## Decision: use parametric test when all decision values exceed the threshold
  use_test_1 <- all(norm_decision > pretest_threshold, na.rm = TRUE)
  ds_result  <- if (use_test_1) fn_for_ds_test_1(data) else fn_for_ds_test_2(data)
  
  ds_result$p.value
}


#' =============================================================================
#' 9.1 Comprehensive downstream test simulation across sample sizes
#'
#' For each distribution and sample size, runs N replicates under H0 and H1,
#' computes power and Type I error rates for test_1, test_2, and the adaptive
#' strategy, and returns results structured for plotting and AUC comparison.
#' =============================================================================

run_ds_simulation <- function(sample_sizes       = c(10, 20, 30, 40, 50),
                              distributions      = c("exponential", "normal"),
                              N                  = 1e3,
                              alpha              = 0.05,
                              effect_size_H1     = 0.5,
                              effect_size_H0     = zero_like(effect_size_H1),
                              norm_config        = list(method = "classical",
                                                        config = list(norm_test = "SW")),
                              pretest_threshold  = 0.05,
                              optimal_result     = NULL,
                              threshold_by_n     = NULL,        # <-- NEW
                              ds_test_methods    = c("test_1", "test_2", "adaptive"),
                              gen_data           = onesample_data,
                              get_parameters     = onesample_parameters,
                              fn_to_get_norm_obj = raw_data,
                              fn_for_ds_test_1   = one_sample_t_test,
                              fn_for_ds_test_2   = sign_test,
                              ...) {
  
  #' ---------------------------------------------------------------------------
  #' Validate inputs
  #' ---------------------------------------------------------------------------
  if (is.null(norm_config$method)) {
    stop("norm_config must have a $method slot: 'classical' or 'custom'.")
  }
  
  ds_test_methods <- match.arg(ds_test_methods, choices   = c("test_1", "test_2", "adaptive"), several.ok = TRUE)
  
  #' ---------------------------------------------------------------------------
  #' pretest_threshold resolution — three-tier priority:
  #'   1. threshold_by_n (named vector, one value per n) — highest priority
  #'   2. optimal_result$param_star (single optimal from Phase 2)
  #'   3. pretest_threshold scalar — fallback default
  #' ---------------------------------------------------------------------------
  use_per_n <- !is.null(threshold_by_n) && is.numeric(threshold_by_n) && !is.null(names(threshold_by_n))
  
  if (!use_per_n) {
    if (!is.null(optimal_result) && !is.null(optimal_result$param_star)) {
      pretest_threshold <- optimal_result$param_star
      message(sprintf("Using optimal pretest threshold = %.4f", pretest_threshold))
    } else {
      warning(sprintf("No optimal_result supplied. Using pretest_threshold = %.4f.", pretest_threshold))
    }
  } else {
    message("Using per-n optimal pretest thresholds.")
  }
  
  #' ---------------------------------------------------------------------------
  #' Capture dots once for safe forwarding inside the parallel body. Parallelism
  #' itself is handled per (distribution, sample size) cell by ds_parallel_map().
  #' ---------------------------------------------------------------------------
  dots <- list(...)
  
  #' ---------------------------------------------------------------------------
  #' Main simulation loop
  #' ---------------------------------------------------------------------------
  sim_results <- list()
  plot_data   <- list()
  
  for (dist in distributions) {
    cat("\n=== Distribution:", dist, "===\n")
    
    dist_results <- list(
      power  = setNames(vector("list", length(ds_test_methods)), ds_test_methods),
      type1  = setNames(vector("list", length(ds_test_methods)), ds_test_methods),
      pvals  = list(
        H0 = setNames(vector("list", length(ds_test_methods)), ds_test_methods),
        H1 = setNames(vector("list", length(ds_test_methods)), ds_test_methods)
      ),
      timing = list()
    )
    for (method in ds_test_methods) {
      dist_results$power[[method]]      <- numeric(length(sample_sizes))
      dist_results$type1[[method]]      <- numeric(length(sample_sizes))
      dist_results$pvals$H0[[method]]   <- vector("list", length(sample_sizes))
      dist_results$pvals$H1[[method]]   <- vector("list", length(sample_sizes))
    }
    
    for (i in seq_along(sample_sizes)) {
      n <- sample_sizes[i]
      
      #' Resolve the threshold for this specific n
      n_threshold <- if (use_per_n) {
        tv <- threshold_by_n[as.character(n)]
        if (is.na(tv)) {
          warning(sprintf("No threshold_by_n entry for n = %d. Falling back to pretest_threshold = %.4f.", n, pretest_threshold))
          pretest_threshold
        } else {
          as.numeric(tv)
        }
      } else {
        pretest_threshold
      }
      pval_H0    <- matrix(NA_real_, nrow = N, ncol = length(ds_test_methods))
      pval_H1    <- matrix(NA_real_, nrow = N, ncol = length(ds_test_methods))
      colnames(pval_H0) <- colnames(pval_H1) <- ds_test_methods
      
      start_time <- Sys.time()
      
      ## Parallelise the N replicates for this (distribution, n) cell.
      ## The progress bar (inside ds_parallel_map) ticks as sim chunks complete.
      sim_rows <- ds_parallel_map(
        N,
        function(sim) {
          paras_H0 <- do.call(get_parameters, c(list(n, effect_size = effect_size_H0, dist = dist), dots))
          data_H0  <- do.call(gen_data, paras_H0)
          paras_H1 <- do.call(get_parameters, c(list(n, effect_size = effect_size_H1, dist = dist), dots))
          data_H1  <- do.call(gen_data, paras_H1)
          
          h0 <- setNames(rep(NA_real_, length(ds_test_methods)), ds_test_methods)
          h1 <- setNames(rep(NA_real_, length(ds_test_methods)), ds_test_methods)
          
          for (method in ds_test_methods) {
            if (method == "test_1") { # test 1
              h0[method] <- fn_for_ds_test_1(data_H0)$p.value
              h1[method] <- fn_for_ds_test_1(data_H1)$p.value
            } else if (method == "test_2") { # test 2
              h0[method] <- fn_for_ds_test_2(data_H0)$p.value
              h1[method] <- fn_for_ds_test_2(data_H1)$p.value
            } else { ## Adaptive test: assess normality then dispatch
              h0[method] <- do.call(run_adaptive_test_dual, c(list(
                data               = data_H0,
                norm_config        = norm_config,
                pretest_threshold  = n_threshold,
                fn_to_get_norm_obj = fn_to_get_norm_obj,
                fn_for_ds_test_1   = fn_for_ds_test_1,
                fn_for_ds_test_2   = fn_for_ds_test_2
              ), dots))
              h1[method] <- do.call(run_adaptive_test_dual, c(list(
                data               = data_H1,
                norm_config        = norm_config,
                pretest_threshold  = n_threshold,
                fn_to_get_norm_obj = fn_to_get_norm_obj,
                fn_for_ds_test_1   = fn_for_ds_test_1,
                fn_for_ds_test_2   = fn_for_ds_test_2
              ), dots))
            }
          }
          list(H0 = h0, H1 = h1)
        },
        label = sprintf("  n = %d (%d sims)", n, N)
      )
      
      ## Reassemble per-replicate rows into the N x methods matrices
      pval_H0[] <- do.call(rbind, lapply(sim_rows, `[[`, "H0"))
      pval_H1[] <- do.call(rbind, lapply(sim_rows, `[[`, "H1"))
      
      for (method in ds_test_methods) {
        dist_results$type1[[method]][i]    <- mean(pval_H0[, method] < alpha, na.rm = TRUE)
        dist_results$power[[method]][i]    <- mean(pval_H1[, method] < alpha, na.rm = TRUE)
        dist_results$pvals$H0[[method]][[i]] <- pval_H0[, method]
        dist_results$pvals$H1[[method]][[i]] <- pval_H1[, method]
      }
      dist_results$timing[[paste0("n=", n)]] <- difftime(Sys.time(), start_time, units = "secs")
    }
    
    dist_results$sample_sizes <- sample_sizes
    sim_results[[dist]]       <- dist_results
    
    plot_data[[dist]] <- list(
      power = data.frame(
        n            = sample_sizes,
        test_1       = dist_results$power$test_1,
        test_2       = dist_results$power$test_2,
        adaptive     = dist_results$power$adaptive,
        Distribution = dist
      ),
      type1 = data.frame(
        n            = sample_sizes,
        test_1       = dist_results$type1$test_1,
        test_2       = dist_results$type1$test_2,
        adaptive     = dist_results$type1$adaptive,
        Distribution = dist
      )
    )
    cat("Completed:", dist, "\n")
  }
  
  #' ---------------------------------------------------------------------------
  #' Combine and compute AUC tables
  #' ---------------------------------------------------------------------------
  combined_power <- do.call(rbind, lapply(plot_data, function(x) x$power))
  combined_type1 <- do.call(rbind, lapply(plot_data, function(x) x$type1))
  rownames(combined_power) <- rownames(combined_type1) <- NULL
  
  auc_tables <- list()
  for (dist in distributions) {
    dp <- plot_data[[dist]]
    auc_tables[[dist]] <- data.frame(
      Method    = ds_test_methods,
      AUC_Power = c(
        compute_auc(sample_sizes, dp$power$test_1,   ensure_endpoints = FALSE, normalize = TRUE),
        compute_auc(sample_sizes, dp$power$test_2,   ensure_endpoints = FALSE, normalize = TRUE),
        compute_auc(sample_sizes, dp$power$adaptive, ensure_endpoints = FALSE, normalize = TRUE)
      ),
      AUC_TypeI = c(
        compute_auc(sample_sizes, dp$type1$test_1,   ensure_endpoints = FALSE, normalize = TRUE),
        compute_auc(sample_sizes, dp$type1$test_2,   ensure_endpoints = FALSE, normalize = TRUE),
        compute_auc(sample_sizes, dp$type1$adaptive, ensure_endpoints = FALSE, normalize = TRUE)
      ),
      row.names = NULL
    )
  }
  
  #' ---------------------------------------------------------------------------
  #' Return
  #' ---------------------------------------------------------------------------
  results <- list(
    raw       = sim_results,
    plot_data = list(combined_power = combined_power, combined_type1 = combined_type1),
    auc_tables = auc_tables,
    summary   = list(
      sample_sizes      = sample_sizes,
      distributions     = distributions,
      N                 = N,
      alpha             = alpha,
      effect_size_H1    = effect_size_H1,
      effect_size_H0    = effect_size_H0,
      norm_config       = norm_config,
      pretest_threshold = pretest_threshold,
      ds_test_methods   = ds_test_methods,
      timestamp         = Sys.time()
    )
  )
  class(results) <- "ds_simulation_results"
  results
}


#' =============================================================================
#' Helper: extract data structured for plotting
#' =============================================================================
extract_for_plotting <- function(results) {
  list(
    combined_power = results$plot_data$combined_power,
    combined_type1 = results$plot_data$combined_type1,
    sample_sizes   = results$summary$sample_sizes
  )
}


#' =============================================================================
#' 9.2 Plot power and Type I error results across sample sizes
#'   Produces a layout with:
#'   Row 1: Power curves per distribution
#'   Row 2: Type I error curves per distribution
#'   Row 3: Legend + title summary
#' =============================================================================

plot_power_type1_results <- function(combined_power,
                                     combined_type1,
                                     ds_test_methods = c("test_1", "test_2", "adaptive"),
                                     distributions   = c("exponential", "normal"),
                                     sample_sizes,
                                     test_alpha      = 0.05,
                                     optimal_result  = NULL,
                                     norm_config     = NULL) {
  
  #' ---------------------------------------------------------------------------
  #' Determine the pretest threshold and approach label for the subtitle
  #' ---------------------------------------------------------------------------
  optimal_parameter <- if (!is.null(optimal_result$param_star))
    optimal_result$param_star else NULL
  
  param_display <- if (!is.null(optimal_parameter)) {
    sprintf("threshold = %.4f", optimal_parameter)
  } else {
    "default threshold"
  }
  
  #' Approach label: read from norm_config$method if available, otherwise
  #' use Custom Approach
  approach_display <- if (!is.null(norm_config$method)) {
    switch(tolower(norm_config$method),
           classical = "Classical Approach",
           custom    = "Custom Approach",
           norm_config$method)
  } else {
    "Custom Approach"
  }
  
  #' ---------------------------------------------------------------------------
  #' Plot aesthetics
  #' ---------------------------------------------------------------------------
  colors <- c("red", "blue", "green")
  shapes <- c(19, 17, 15)
  
  #' ---------------------------------------------------------------------------
  #' Layout: 2 plot rows + 1 legend row
  #' ---------------------------------------------------------------------------
  n_dist       <- length(distributions)
  n_panels     <- 2 * n_dist          
  legend_panel <- n_panels + 1
  layout(
    matrix(c(seq_len(n_panels), rep(legend_panel, n_dist)), nrow = 3, byrow = TRUE),
    heights = c(1, 1, 0.25)
  )
  
  op <- par(mar  = c(3.4, 3.4, 2.7, 0.8),
            oma  = c(0, 0, 2.5, 0),
            mgp  = c(2.0, 0.50, 0),
            tcl  = -0.2,
            xaxs = "r", yaxs = "r",
            cex.axis = 1.0, cex.lab = 1.2,
            font.lab  = 2,   # axis titles
            font.axis = 2    # axis tick labels 
  )
  
  on.exit(par(op))
  
  x_range <- range(sample_sizes, finite = TRUE)
  x_pad   <- 0.02 * diff(x_range)
  x_lim   <- c(x_range[1] - x_pad, x_range[2] + x_pad)
  
  #' ---------------------------------------------------------------------------
  #' Row 1: Power curves
  #' ---------------------------------------------------------------------------
  for (dist in distributions) {
    dist_power <- subset(combined_power, Distribution == dist)
    
    plot(NA, xlim = x_lim, ylim = c(0, 1.01), xlab = "Sample Size", ylab = "Power", main = paste("Power -", dist), cex.lab = 1.5, cex.main = 1.5, font.main = 2)
    
    for (j in seq_along(ds_test_methods)) { 
      method <- ds_test_methods[j]
      lines(dist_power$n, dist_power[[method]], type = "b", col = colors[j], pch = shapes[j], lwd = 3)
    }
  }
  
  #' ---------------------------------------------------------------------------
  #' Row 2: Type I error curves
  #' ---------------------------------------------------------------------------
  for (dist in distributions) {
    dist_type1 <- subset(combined_type1, Distribution == dist)
    
    y_max <- max(dist_type1[, ds_test_methods], na.rm = TRUE)
    y_pad <- 0.05 * y_max
    
    plot(NA, xlim = x_lim, ylim = c(0, y_max + max(y_pad, 1e-6)),
         xlab = "Sample Size", ylab = "P(Type I Error)",
         main = paste("P(Type I Error) -", dist),
         cex.lab = 1.5, cex.main = 1.5, font.main = 2)
    
    for (j in seq_along(ds_test_methods)) {
      method <- ds_test_methods[j]
      lines(dist_type1$n, dist_type1[[method]], type = "b", col = colors[j], pch = shapes[j], lwd = 3)
    }
    
    abline(h = test_alpha, col = "gray50", lty = 3, lwd = 1.2)
  }
  
  #' ---------------------------------------------------------------------------
  #' Row 3: Legend
  #' ---------------------------------------------------------------------------
  par(mar = c(0.1, 0.1, 0.1, 0.1))
  plot.new()
  
  legend("center",
         legend = ds_test_methods,
         col    = colors,
         title  = "Test Methods",
         pch    = shapes,
         lwd    = 3,
         horiz  = TRUE,
         bty    = "n",
         cex    = 1.5)
  
  mtext(
    sprintf("Test Methods Comparison (%s) | %s", approach_display, param_display),
    side = 3, outer = TRUE, line = -0.25, cex = 1.5
  )
}


#' ==============================================================================
#' 10.0 Power vs effect size analysis
#' For each distribution and effect size, runs Nsim replicates under H1 and
#' computes power for test_1, test_2, and the adaptive strategy at a fixed sample size.
#' ==============================================================================

perform_ds_power_by_effect <- function(fixed_n           = 10,
                                       effect_sizes      = seq(0.1, 1.0, by = 0.1),
                                       distributions     = c("exponential", "normal"),
                                       Nsim              = 1e3,
                                       norm_config       = list(method = "classical",
                                                                config = list(norm_test = "SW")),
                                       pretest_threshold = 0.05,
                                       optimal_result    = NULL,
                                       gen_data,
                                       get_parameters,
                                       fn_to_get_norm_obj,
                                       fn_for_ds_test_1,
                                       fn_for_ds_test_2,
                                       ds_test_methods   = c("test_1", "test_2", "adaptive"),
                                       test_alpha        = 0.05,
                                       effect_template   = NULL,
                                       verbose           = TRUE,
                                       ...) {
  
  #' ---------------------------------------------------------------------------
  #' Validate
  #' ---------------------------------------------------------------------------
  if (is.null(norm_config$method)) {
    stop("norm_config must have a $method slot: 'classical', or 'custom'.")
  }
  
  ds_test_methods <- match.arg(ds_test_methods,
                               choices    = c("test_1", "test_2", "adaptive"),
                               several.ok = TRUE)
  
  #' ---------------------------------------------------------------------------
  #' Resolve pretest_threshold
  #' ---------------------------------------------------------------------------
  if (!is.null(optimal_result) && !is.null(optimal_result$param_star)) {
    pretest_threshold <- optimal_result$param_star
    message(sprintf("Using optimal pretest threshold = %.4f", pretest_threshold))
  } else {
    warning(sprintf(
      "No optimal_result supplied. Using pretest_threshold = %.4f.", pretest_threshold
    ))
  }
  
  #' ---------------------------------------------------------------------------
  #' Infer effect_template if not supplied
  #' ---------------------------------------------------------------------------
  if (is.null(effect_template)) {
    base_paras      <- get_parameters(fixed_n, dist = distributions[1], ...)
    effect_template <- base_paras$effect_size
  }
  eff_len <- length(effect_template)
  
  #' ---------------------------------------------------------------------------
  #' Parallel backend (FORK on macOS/Linux via doMC: workers inherit attached
  #' packages and functions, so no .packages / .export is required). Registered
  #' once; the inner Nsim replicate loop is parallelised with %dorng%. The
  #' existing effect-size progress bar is advanced once per effect size.
  #' ---------------------------------------------------------------------------
  dots    <- list(...)
  n_cores <- getOption("ds.parallel.cores", max(1L, parallel::detectCores() - 1L))
  if (.Platform$OS.type != "windows" && requireNamespace("doMC", quietly = TRUE)) {
    doMC::registerDoMC(n_cores)
  } else {
    foreach::registerDoSEQ()
  }
  
  #' ---------------------------------------------------------------------------
  #' Main loop
  #' ---------------------------------------------------------------------------
  results_by_dist <- list()
  
  for (dist in distributions) {
    if (verbose) cat("\n=== Distribution:", dist, "===\n")
    
    results_power <- setNames(vector("list", length(ds_test_methods)), ds_test_methods)
    timing        <- setNames(vector("list", length(ds_test_methods)), ds_test_methods)
    pval_storage  <- list(
      H1 = setNames(vector("list", length(ds_test_methods)), ds_test_methods)
    )
    
    for (method in ds_test_methods) {
      if (verbose) cat("  Method:", method, "\n")
      
      pow_vec <- numeric(length(effect_sizes))
      names(pow_vec) <- paste0("d=", effect_sizes)
      pval_storage$H1[[method]] <- vector("list", length(effect_sizes))
      
      pb <- txtProgressBar(min = 0, max = length(effect_sizes), style = 3)
      t0 <- Sys.time()
      
      for (i in seq_along(effect_sizes)) {
        d <- effect_sizes[i]
        
        ## Build effect_size argument with correct structure
        current_effect <- if (eff_len > 1) {
          e         <- rep(0, eff_len)
          e[eff_len] <- d
          e
        } else d
        
        ## Parallelise the Nsim replicates for this effect size.
        ## The effect-size progress bar (pb) is advanced once per effect size,
        ## exactly as before.
        pval_H1 <- foreach::foreach(sim = seq_len(Nsim), .combine = c) %dorng% {
          paras_H1 <- do.call(get_parameters, c(list(fixed_n, effect_size = current_effect, dist = dist), dots))
          data_H1  <- do.call(gen_data, paras_H1)
          
          if (method == "test_1") {
            fn_for_ds_test_1(data_H1)$p.value
          } else if (method == "test_2") {
            fn_for_ds_test_2(data_H1)$p.value
          } else {
            do.call(run_adaptive_test_dual, c(list(
              data               = data_H1,
              norm_config        = norm_config,
              pretest_threshold  = pretest_threshold,
              fn_to_get_norm_obj = fn_to_get_norm_obj,
              fn_for_ds_test_1   = fn_for_ds_test_1,
              fn_for_ds_test_2   = fn_for_ds_test_2
            ), dots))
          }
        }
        
        pow_vec[i]                     <- mean(pval_H1 < test_alpha, na.rm = TRUE)
        pval_storage$H1[[method]][[i]] <- pval_H1
        setTxtProgressBar(pb, i)
      }
      close(pb)
      
      timing[[method]]        <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      results_power[[method]] <- pow_vec
      if (verbose) cat("  Completed in", round(timing[[method]], 2), "seconds\n")
    }
    
    results_by_dist[[dist]] <- list(
      power             = results_power,
      pvalues           = pval_storage,
      timing            = timing,
      effect_sizes      = effect_sizes,
      pretest_threshold = pretest_threshold,
      norm_config       = norm_config,
      fixed_n           = fixed_n,
      distribution      = dist,
      test_alpha        = test_alpha,
      eff_len           = eff_len
    )
  }
  
  format_power_effect_results(results_by_dist, distributions, ds_test_methods)
}


#' =============================================================================
#' 10.1 Format power-effect results for plotting
#' (internal helper called by perform_ds_power_by_effect)
#' =============================================================================

format_power_effect_results <- function(results_by_dist, distributions, methods) {
  
  plot_data <- list()
  
  for (dist in distributions) {
    dist_results <- results_by_dist[[dist]]
    
    power_df <- data.frame(
      effect_value = dist_results$effect_sizes,
      Distribution = dist,
      n            = dist_results$fixed_n,
      stringsAsFactors = FALSE
    )
    for (method in methods) {
      power_df[[method]] <- dist_results$power[[method]]
    }
    plot_data[[dist]] <- power_df
  }
  
  combined_power <- do.call(rbind, plot_data)
  rownames(combined_power) <- NULL
  
  ## Read shared settings from the first distribution's stored results
  first <- results_by_dist[[distributions[1]]]
  
  list(
    raw_results = results_by_dist,
    plot_data   = list(power_by_effect = combined_power),
    parameters  = list(
      distributions     = distributions,
      methods           = methods,
      norm_config       = first$norm_config,
      pretest_threshold = first$pretest_threshold,
      fixed_n           = first$fixed_n,
      effect_sizes      = first$effect_sizes,
      test_alpha        = first$test_alpha,
      eff_len           = first$eff_len,
      timestamp         = Sys.time()
    )
  )
}


#' =============================================================================
#' 10.2 Plot power vs effect size
#' =============================================================================
plot_power_by_effect_size <- function(combined_power,
                                      fixed_n,
                                      optimal_result = NULL,
                                      norm_config    = NULL,
                                      eff_len        = 1) {
  
  #' ---------------------------------------------------------------------------
  #' Threshold and approach display — same pattern as plot_power_type1_results
  #' ---------------------------------------------------------------------------
  optimal_parameter <- if (!is.null(optimal_result$param_star))
    optimal_result$param_star else NULL
  
  param_display <- if (!is.null(optimal_parameter)) {
    sprintf("threshold = %.4f", optimal_parameter)
  } else {
    "default threshold"
  }
  
  #' Approach label: read from norm_config$method if available, otherwise
  #' use Custom Approach
  approach_display <- if (!is.null(norm_config$method)) {
    switch(tolower(norm_config$method),
           classical = "Classical Approach",
           custom    = "Custom Approach",
           norm_config$method)
  } else {
    "Custom Approach"
  }
  
  #' ---------------------------------------------------------------------------
  #' Plot aesthetics
  #' ---------------------------------------------------------------------------
  colors <- c("red", "blue", "green")
  shapes <- c(16, 17, 15)
  
  distributions <- unique(combined_power$Distribution)
  
  x_label <- if (eff_len > 1) {
    paste("Effect Size of Group", eff_len)
  } else {
    "Effect Size"
  }
  
  #' ---------------------------------------------------------------------------
  #' Layout: one panel per distribution + legend row
  #' ---------------------------------------------------------------------------
  n_dist <- length(distributions)
  layout(
    matrix(c(seq_len(n_dist), rep(n_dist + 1, n_dist)), nrow = 2, byrow = TRUE),
    heights = c(4, 0.45)
  )
  
  par(mar = c(4, 4, 3, 1), 
      oma = c(1, 0, 3, 0),
      cex.axis = 1.0,
      cex.lab  = 1.2,
      font.lab  = 2,   # axis titles
      font.axis = 2    # axis tick labels (numbers)
  )
  
  #' ---------------------------------------------------------------------------
  #' Power panels
  #' ---------------------------------------------------------------------------
  for (dist in distributions) {
    dist_data <- combined_power[combined_power$Distribution == dist, ]
    
    plot(NA,
         xlim = range(dist_data$effect_value),
         ylim = c(0, 1),
         xlab = x_label,
         ylab = "Power",
         main = paste(dist, "(n =", fixed_n, "per group)"))
    
    lines(dist_data$effect_value, dist_data$test_1,
          type = "b", col = colors[1], lwd = 3, pch = shapes[1])
    lines(dist_data$effect_value, dist_data$test_2,
          type = "b", col = colors[2], lwd = 3, pch = shapes[2])
    lines(dist_data$effect_value, dist_data$adaptive,
          type = "b", col = colors[3], lwd = 3, pch = shapes[3])
    
    abline(v = 0.5, col = "gray", lty = 2)
  }
  
  #' ---------------------------------------------------------------------------
  #' Outer title
  #' ---------------------------------------------------------------------------
  mtext(
    sprintf("Power vs %s (%s) | %s", x_label, approach_display, param_display),
    side = 3, outer = TRUE, line = 1, cex = 1.2, font = 2
  )
  
  #' ---------------------------------------------------------------------------
  #' Legend panel
  #' ---------------------------------------------------------------------------
  par(mar = c(0.2, 0.2, 0.2, 0.2))
  plot.new()
  legend("center",
         title  = "Test Methods",
         legend = c("test_1", "test_2", "adaptive"),
         col    = colors,
         lwd    = 3,
         pch    = shapes,
         horiz  = TRUE,
         bty    = "n",
         cex    = 1)
}


#' =============================================================================
#' run_simulation()
#'
#' Orchestrates the full pipeline in six phases:
#'   Phase 1 : Normality pretest ROC curves
#'   Phase 2 : Trade-off analysis (EPL, EPG, ETIE) → selects optimal threshold
#'   Phase 3 : Power vs Type I error ROC at the optimal threshold
#'   Phase 4 : Power and Type I error across sample sizes
#'   Phase 5 : Power vs effect size at a fixed sample size
#'   Phase 6 : Save results and workspace to disk
#'
#' Phase dependencies
#' ------------------
#'   Phase 1  — independent; can always be skipped
#'   Phase 2  — independent; produces param_star used by Phases 3–5
#'   Phase 3  — needs param_star (from Phase 2 or pretest_threshold)
#'   Phase 4  — needs param_star (from Phase 2 or pretest_threshold)
#'   Phase 5  — needs param_star (from Phase 2 or pretest_threshold)
#'   Phase 6  — always runs when any phase ran; saves whatever was produced
#'
#' Skipping rules
#' --------------
#'    Any phase can be excluded via the `phases` argument.
#'    If Phases 3, 4, or 5 are requested and Phase 2 is skipped, the user
#'     MUST supply `pretest_threshold`. An error is raised otherwise.
#'    Phase 6 (save) cannot be skipped independently — it always executes
#'     whenever at least one other phase ran.
#'
#' Key parameters
#' ------------------
#'   phases            : integer vector of phases to run, e.g. c(1, 2, 4).
#'                       Default c(1,2,3,4,5,6) runs everything.
#'   pretest_threshold : numeric scalar. Required when Phase 2 is skipped and
#'                       any of Phases 3–5 are requested. Ignored otherwise.
#' =============================================================================

run_simulation <- function(Nsim           = 1e3,
                           N_tradeoff     = 1e3,
                           test_type      = "one_sample_t_vs_sign_test",
                           distributions  = c("exponential", "normal"),
                           
                           norm_config    = list(method = "classical",
                                                 config = list(norm_test = "SW")),
                           threshold_grid = seq(0.009, 1, by = 0.025),
                           
                           tol_pos    = 0.01,
                           loss_tol   = 0.01,
                           test_alpha = 0.05,
                           center_by  = "median",
                           effect_size = 0.5,
                           
                           sample_sizes = c(10, 20, 30, 40, 50),
                           single_n     = NULL,
                           
                           effect_sizes_plot = c(0.0, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0),
                           sig_levels        = seq(0.005, 1, by = 0.005),
                           
                           ds_test_methods = c("test_1", "test_2", "adaptive"),
                           
                           norm_test      = c("SW","SF","LF","KS","JB","SKEW", "KURT",
                                              "DAP","AD","CVM"),
                           selected_tests = c("SW","SF","LF","JB","DAP","AD"),
                           
                           phases            = c(1, 2, 3, 4, 5, 6),
                           pretest_threshold = NULL,
                           per_n_thresholds  = FALSE,    
                           
                           gen_data           = onesample_data,
                           get_parameters     = onesample_parameters,
                           fn_to_get_norm_obj = raw_data,
                           fn_for_ds_test_1   = one_sample_t_test,
                           fn_for_ds_test_2   = sign_test,
                           ...) {
  
  #' --------------------------------------------------------------------------
  #' STEP 1: VALIDATE
  #' --------------------------------------------------------------------------
  test_type <- trimws(test_type)
  method    <- tolower(trimws(norm_config$method %||% ""))
  
  if (!method %in% c("classical", "custom")) {
    stop("norm_config$method must be 'classical' or 'custom'.")
  }
  
  #' --------------------------------------------------------------------------
  #' STEP 2: VALIDATE PHASES AND DEPENDENCIES
  #' --------------------------------------------------------------------------
  valid_phases <- 1:6
  bad_phases   <- setdiff(phases, valid_phases)
  if (length(bad_phases) > 0) {
    stop("Invalid phase(s): ", paste(bad_phases, collapse = ", "),
         ". Valid phases are 1 through 6.")
  }
  
  run_phase <- function(p) p %in% phases
  
  #' Phases 3, 4, 5 all need param_star.
  #' It comes from Phase 2 unless the user bypasses Phase 2 with pretest_threshold.
  needs_threshold <- any(c(3, 4, 5) %in% phases)
  skip_phase2     <- !run_phase(2)
  
  if (needs_threshold && skip_phase2) {
    if (is.null(pretest_threshold) || !is.numeric(pretest_threshold)) {
      stop(
        "Phase 2 (trade-off analysis) is skipped but Phases ",
        paste(intersect(c(3, 4, 5), phases), collapse = ", "),
        " are requested.\n",
        "These phases require a pretest threshold.\n",
        "Please supply: pretest_threshold = <scalar> or ",
        "pretest_threshold = c('10' = 0.05, '20' = 0.04, ...)."
      )
    }
    if (!is.null(names(pretest_threshold))) {
      message("Phase 2 skipped. Using user-supplied per-n pretest thresholds.")
    } else {
      message(sprintf(
        "Phase 2 skipped. Using user-supplied pretest_threshold = %.4f.",
        pretest_threshold
      ))
    }
  }
  
  #' --------------------------------------------------------------------------
  #' STEP 3: TIMING / INIT
  #' --------------------------------------------------------------------------
  start_time     <- Sys.time()
  method_display <- switch(method,
                           classical = "Classical",
                           custom    = "Custom")
  
  cat("\n", strrep("=", 60), "\n", sep = "")
  cat("STARTING SIMULATION WORKFLOW (", method_display, " Approach)\n", sep = "")
  cat("Test type    :", test_type, "\n")
  cat("Distributions:", paste(distributions, collapse = " vs "), "\n")
  cat("Phases       :", paste(sort(phases), collapse = ", "), "\n")
  cat(strrep("=", 60), "\n\n")
  
  eff       <- split_effect_size(effect_size)
  effect_H0 <- eff$H0
  effect_H1 <- eff$H1
  
  #' --------------------------------------------------------------------------
  #' Results container
  #' --------------------------------------------------------------------------
  results <- list(
    params = list(
      Nsim              = Nsim,
      N_tradeoff        = N_tradeoff,
      sample_sizes      = sample_sizes,
      single_n          = single_n,
      test_alpha        = test_alpha,
      norm_config       = norm_config,
      threshold_grid    = threshold_grid,
      effect_sizes_plot = effect_sizes_plot,
      sig_levels        = sig_levels,
      ds_test_methods   = ds_test_methods,
      test_type         = test_type,
      distributions     = distributions,
      tol_pos           = tol_pos,
      loss_tol          = loss_tol,
      effect_size_input = effect_size,
      effect_size_H0    = effect_H0,
      effect_size_H1    = effect_H1,
      center_by         = center_by,
      phases_run        = sort(phases),
      pretest_threshold = pretest_threshold,
      timestamp_start   = start_time
    ),
    objects = list(),
    files   = list()
  )
  
  results$root_dir <- file.path(getwd(), "results", test_type)
  if (!dir.exists(results$root_dir)) {
    dir.create(results$root_dir, recursive = TRUE, showWarnings = FALSE)
  }
  message("Output directory: ", results$root_dir)
  
  #' --------------------------------------------------------------------------
  #' Helper: resolve param_star for Phases 3–5
  #' Either from Phase 2 output or from user-supplied pretest_threshold.
  #' Called only when Phase 2 was skipped.
  #' --------------------------------------------------------------------------
  make_bypass_optimal <- function(threshold) {
    list(
      param_star   = threshold,
      param_name   = "threshold",
      selection_note = sprintf(
        "User-supplied pretest_threshold = %.4f (Phase 2 skipped)", threshold
      ),
      feasible     = NA,
      power_loss   = NA,
      power_gain   = NA,
      inflation_normal     = NA,
      inflation_non_normal = NA,
      all_metrics  = list(
        param_values  = threshold,
        distributions = distributions,
        param_name    = "threshold"
      )
    )
  }
  
  #' ==========================================================================
  #' PHASE 1: NORMALITY PRETEST ROC ANALYSIS
  #' ==========================================================================
  if (run_phase(1)) {
    cat("\n1. GENERATING NORMALITY TEST ROC CURVES\n")
    
    results$files$norm_roc <- file.path(
      results$root_dir, paste0(test_type, "_normality_roc.pdf")
    )
    pdf(results$files$norm_roc, width = 9, height = 7)
    
    phase1_out <- tryCatch({
      plot_normality_roc_wrapper(
        norm_config        = norm_config,
        single_n           = single_n,
        threshold_grid     = threshold_grid,
        distributions      = distributions,
        Nsim               = Nsim,
        norm_test          = norm_test,
        selected_tests     = selected_tests,
        gen_data           = gen_data,
        get_parameters     = get_parameters,
        fn_to_get_norm_obj = fn_to_get_norm_obj,
        center_by          = center_by,
        ...
      )
    }, error = function(e) {
      message("  Phase 1 plotting error: ", conditionMessage(e))
      NULL
    }, finally = {
      dev.off()
    })
    
    results$objects$roc_pval_ds_test <- phase1_out
    if (!is.null(phase1_out)) {
      cat("  Saved:", results$files$norm_roc, "\n")
    } else {
      cat("  Phase 1 ROC plot could not be saved — check the error above.\n")
    }
  } else {
    cat("\n1. NORMALITY ROC CURVES  [skipped]\n")
  }
  
  #' ==========================================================================
  #' PHASE 2: TRADE-OFF ANALYSIS
  #' ==========================================================================
  if (run_phase(2)) {
    
    #' ------------------------------------------------------------------
    #' CASE A: per-n thresholds — run trade-off analysis for every n
    #'
    #' Each sample size requires an independent N_tradeoff-replicate
    #' simulation to find its optimal threshold. These runs share no
    #' state. A sequential foreach loop is used here so that:
    #'   (a) a progress bar is visible in the parent process — %dopar%
    #'       runs workers in separate R sessions where txtProgressBar
    #'       output cannot reach the console;
    #'   (b) the pdf() device can be opened and written inside the loop
    #'       without process-safety concerns.
    #' The inner perform_adaptive_analysis() loop is already parallelised
    #' via pbapply, so the ML prediction load is distributed there.
    #' ------------------------------------------------------------------
    if (per_n_thresholds) {
      cat("\n2. PERFORMING TRADE-OFF ANALYSIS (per sample size)\n")
      
      ## Capture dots (...) as a list so they can be passed inside foreach.
      ## foreach workers cannot access '...' from the enclosing scope directly.
      extra_args <- list(...)
      
      ## Progress bar: one tick per completed sample size.
      pb_tradeoff <- txtProgressBar(
        min   = 0,
        max   = length(sample_sizes),
        style = 3
      )
      on.exit(close(pb_tradeoff), add = TRUE)
      
      ## Open the trade-off PDF before the loop so each plot is written
      ## as its sample size completes (sequential, device-safe).
      results$files$tradeoff <- file.path(
        results$root_dir, paste0(test_type, "_tradeoff_analysis.pdf")
      )
      pdf(results$files$tradeoff, width = 10, height = 8)
      
      ## Run one trade-off analysis per sample size sequentially.
      ## %do% (sequential foreach) keeps the same result structure as
      ## the parallel version so the collection code below is unchanged.
      per_n_results <- foreach(
        n_i            = sample_sizes,
        .errorhandling = "stop"
      ) %do% {
        
        ## Tradeoff simulation for this specific sample size.
        analysis_i <- do.call(
          perform_adaptive_analysis,
          c(
            list(
              Nsim               = N_tradeoff,
              n                  = n_i,
              effect_size_H1     = effect_H1,
              effect_size_H0     = effect_H0,
              distributions      = distributions,
              norm_config        = norm_config,
              threshold_grid     = threshold_grid,
              test_alpha         = test_alpha,
              gen_data           = gen_data,
              get_parameters     = get_parameters,
              fn_to_get_norm_obj = fn_to_get_norm_obj,
              fn_for_ds_test_1   = fn_for_ds_test_1,
              fn_for_ds_test_2   = fn_for_ds_test_2,
              center_by          = center_by
            ),
            extra_args
          )
        )
        
        ## Compute ROC metrics from the simulation output.
        metrics_i <- compute_roc_metrics(
          typeI_rates = analysis_i$typeI_rates,
          power_rates = analysis_i$power_rates,
          test_alpha  = test_alpha
        )
        
        ## Select the optimal pretest threshold for this sample size.
        optimal_i <- select_optimal_parameter(
          distributions = distributions,
          param_grid    = analysis_i$param_grid,
          metrics       = metrics_i,
          alpha         = test_alpha,
          tol_pos       = tol_pos,
          loss_tol      = loss_tol,
          param_name    = analysis_i$param_name
        )
        
        ## Write the plot for this sample size immediately.
        plot_tradeoff_results_generic(
          optimal_result = optimal_i,
          outer_title    = sprintf(
            "Trade-off (%s Approach) | n = %d | threshold = %.4f",
            method_display, n_i, optimal_i$param_star
          )
        )
        
        ## Advance the progress bar after this sample size completes.
        setTxtProgressBar(pb_tradeoff, which(sample_sizes == n_i))
        
        ## Return all three objects for collection below.
        list(
          n_i      = n_i,
          analysis = analysis_i,
          metrics  = metrics_i,
          optimal  = optimal_i
        )
      }
      
      dev.off()
      
      ## %do% returns results in order, but sort defensively for consistency.
      per_n_results <- per_n_results[order(sapply(per_n_results, `[[`, "n_i"))]
      
      ## Collect optimal thresholds and result objects into named vectors/lists.
      threshold_by_n      <- setNames(
        sapply(per_n_results, function(r) r$optimal$param_star),
        as.character(sample_sizes)
      )
      optimal_result_by_n <- setNames(
        lapply(per_n_results, `[[`, "optimal"),
        as.character(sample_sizes)
      )
      
      cat("\n")
      cat("  Saved:", results$files$tradeoff, "\n")
      cat("  Per-n thresholds:\n")
      for (i in seq_along(sample_sizes)) {
        cat(sprintf("    n = %d  ->  threshold = %.4f\n",
                    sample_sizes[i], threshold_by_n[i]))
      }
      
      ## For Phases 3 and 5, which operate at single_n, extract that entry.
      optimal_result <- optimal_result_by_n[[as.character(single_n)]]
      param_star     <- threshold_by_n[as.character(single_n)]
      
      ## Expose the single_n trade-off analysis (which holds its all_pvalues)
      ## so Phase 3 can reuse that simulation instead of regenerating it.
      single_n_idx <- which(sapply(per_n_results, `[[`, "n_i") == single_n)
      if (length(single_n_idx) == 1L) {
        results$objects$analysis_ds_tests <- per_n_results[[single_n_idx]]$analysis
      }
      
      results$objects$threshold_by_n     <- threshold_by_n
      results$objects$optimal_result_by_n <- optimal_result_by_n
      results$objects$optimal_result     <- optimal_result
      results$objects$param_star         <- param_star
      
      #' ------------------------------------------------------------------
      #' CASE B: single threshold — existing behaviour unchanged
      #' ------------------------------------------------------------------
    } else {
      cat("\n2. PERFORMING TRADE-OFF ANALYSIS\n")
      
      analysis_ds_tests <- perform_adaptive_analysis(
        Nsim               = N_tradeoff,
        n                  = single_n,
        effect_size_H1     = effect_H1,
        effect_size_H0     = effect_H0,
        distributions      = distributions,
        norm_config        = norm_config,
        threshold_grid     = threshold_grid,
        test_alpha         = test_alpha,
        gen_data           = gen_data,
        get_parameters     = get_parameters,
        fn_to_get_norm_obj = fn_to_get_norm_obj,
        fn_for_ds_test_1   = fn_for_ds_test_1,
        fn_for_ds_test_2   = fn_for_ds_test_2,
        center_by          = center_by,
        ...
      )
      results$objects$analysis_ds_tests <- analysis_ds_tests
      
      metrics <- compute_roc_metrics(
        typeI_rates = analysis_ds_tests$typeI_rates,
        power_rates = analysis_ds_tests$power_rates,
        test_alpha  = test_alpha
      )
      results$objects$metrics <- metrics
      
      optimal_result <- select_optimal_parameter(
        distributions = distributions,
        param_grid    = analysis_ds_tests$param_grid,
        metrics       = metrics,
        alpha         = test_alpha,
        tol_pos       = tol_pos,
        loss_tol      = loss_tol,
        param_name    = analysis_ds_tests$param_name
      )
      results$objects$optimal_result <- optimal_result
      param_star <- optimal_result$param_star
      results$objects$param_star <- param_star
      threshold_by_n <- NULL
      
      results$files$tradeoff <- file.path(
        results$root_dir, paste0(test_type, "_tradeoff_analysis.pdf")
      )
      pdf(results$files$tradeoff, width = 10, height = 8)
      plot_tradeoff_results_generic(
        optimal_result = optimal_result,
        outer_title    = sprintf("Power and Type I Error Trade-off (%s Approach)", method_display)
      )
      dev.off()
      cat("  Saved:", results$files$tradeoff, "\n")
      cat("  Optimal threshold =", round(param_star, 4), "\n")
    }
    
  } else {
    cat("\n2. TRADE-OFF ANALYSIS  [skipped]\n")
    threshold_by_n <- NULL
    
    if (needs_threshold) {
      ## Accept either a named vector or a scalar for pretest_threshold
      if (is.numeric(pretest_threshold) && !is.null(names(pretest_threshold))) {
        threshold_by_n    <- pretest_threshold
        param_star        <- pretest_threshold[as.character(single_n)] %||%
          pretest_threshold[[1]]
        optimal_result    <- make_bypass_optimal(as.numeric(param_star))
        cat("  Using user-supplied per-n thresholds.\n")
        for (nm in names(threshold_by_n)) {
          cat(sprintf(" n = %s  ->  threshold = %.4f\n", nm, threshold_by_n[nm]))
        }
      } else {
        optimal_result <- make_bypass_optimal(pretest_threshold)
        param_star     <- pretest_threshold
        cat("  Using pretest_threshold =", round(param_star, 4), "\n")
      }
      results$objects$optimal_result <- optimal_result
      results$objects$param_star     <- param_star
      results$objects$threshold_by_n <- threshold_by_n
    }
  }
  
  #' ==========================================================================
  #' PHASE 3: POWER VS TYPE I ERROR ROC
  #' ==========================================================================
  if (run_phase(3)) {
    cat("\n3. GENERATING POWER VS TYPE I ERROR ROC DATA\n")
    
    #' Reuse Phase 2's simulation when available. analysis_ds_tests$all_pvalues
    #' is a per-distribution list of generate_pval() outputs produced at the
    #' same single_n, effect sizes, and norm_config that Phase 3 needs — so it
    #' can be consumed verbatim. NULL when Phase 2 was skipped, in which case
    #' power_vs_error_roc_data() falls back to simulating.
    sim_cache <- results$objects$analysis_ds_tests$all_pvalues
    if (!is.null(sim_cache)) {
      cat("  Reusing Phase 2 generate_pval() cache for distributions:",
          paste(names(sim_cache), collapse = ", "), "\n")
    }
    
    roc_data <- power_vs_error_roc_data(
      N                  = Nsim,
      n                  = single_n,
      distributions      = distributions,
      norm_config        = norm_config,
      effect_size_H1     = effect_H1,
      effect_size_H0     = effect_H0,
      optimal_result     = optimal_result,
      sig_levels         = sig_levels,
      sim_cache          = sim_cache,
      gen_data           = gen_data,
      get_parameters     = get_parameters,
      fn_to_get_norm_obj = fn_to_get_norm_obj,
      fn_for_ds_test_1   = fn_for_ds_test_1,
      fn_for_ds_test_2   = fn_for_ds_test_2,
      center_by          = center_by,
      ...
    )
    results$objects$roc_data <- roc_data
    
    results$files$power_error_roc <- file.path(
      results$root_dir, paste0(test_type, "_power_error_roc.pdf")
    )
    pdf(results$files$power_error_roc, width = 10, height = 10)
    power_vs_error_roc_plot(roc_results   = roc_data,
                            nominal_alpha  = test_alpha,
                            optimal_result = optimal_result)
    dev.off()
    cat("  Saved:", results$files$power_error_roc, "\n")
    
  } else {
    cat("\n3. POWER VS TYPE I ERROR ROC  [skipped]\n")
  }
  
  #' ==========================================================================
  #' PHASE 4: MAIN SIMULATION ACROSS SAMPLE SIZES
  #' ==========================================================================
  if (run_phase(4)) {
    cat("\n4. RUNNING MAIN SIMULATION (ACROSS SAMPLE SIZES)\n")
    
    sim_output <- run_ds_simulation(
      sample_sizes       = sample_sizes,
      distributions      = distributions,
      N                  = Nsim,
      alpha              = test_alpha,
      effect_size_H1     = effect_H1,
      effect_size_H0     = effect_H0,
      norm_config        = norm_config,
      optimal_result     = optimal_result,
      threshold_by_n     = results$objects$threshold_by_n,
      ds_test_methods    = ds_test_methods,
      gen_data           = gen_data,
      get_parameters     = get_parameters,
      fn_to_get_norm_obj = fn_to_get_norm_obj,
      fn_for_ds_test_1   = fn_for_ds_test_1,
      fn_for_ds_test_2   = fn_for_ds_test_2,
      center_by          = center_by,
      ...
    )
    results$objects$sim_output       <- sim_output
    results$objects$combined_power_n <- sim_output$plot_data$combined_power
    results$objects$combined_type1_n <- sim_output$plot_data$combined_type1
    results$objects$auc_tables       <- sim_output$auc_tables
    
    results$files$power_type1 <- file.path(
      results$root_dir, paste0(test_type, "_power_type1_comparison.pdf")
    )
    pdf(results$files$power_type1, width = 10, height = 8)
    plot_power_type1_results(
      combined_power  = sim_output$plot_data$combined_power,
      combined_type1  = sim_output$plot_data$combined_type1,
      ds_test_methods = ds_test_methods,
      distributions   = distributions,
      sample_sizes    = sample_sizes,
      test_alpha      = test_alpha,
      optimal_result  = optimal_result,
      norm_config     = norm_config
    )
    dev.off()
    cat("  Saved:", results$files$power_type1, "\n")
    
    ## Print AUC summary to console
    cat("  AUC Summary:\n")
    for (dist in names(sim_output$auc_tables)) {
      cat("  Distribution:", dist, "\n")
      print(sim_output$auc_tables[[dist]], row.names = FALSE)
      cat("\n")
    }
    
  } else {
    cat("\n4. MAIN SIMULATION (SAMPLE SIZES)  [skipped]\n")
  }
  
  #' ==========================================================================
  #' PHASE 5: POWER VS EFFECT SIZE
  #' ==========================================================================
  if (run_phase(5)) {
    cat("\n5. RUNNING POWER VS EFFECT SIZE ANALYSIS\n")
    
    base_paras      <- get_parameters(single_n, dist = distributions[1],
                                      center_by = center_by, ...)
    eff_len         <- length(base_paras$effect_size)
    effect_template <- if (eff_len > 1) rep(0, eff_len) else 0
    
    sim_power_list <- perform_ds_power_by_effect(
      fixed_n            = single_n,
      effect_sizes       = effect_sizes_plot,
      distributions      = distributions,
      Nsim               = Nsim,
      norm_config        = norm_config,
      optimal_result     = optimal_result,
      gen_data           = gen_data,
      get_parameters     = get_parameters,
      fn_to_get_norm_obj = fn_to_get_norm_obj,
      fn_for_ds_test_1   = fn_for_ds_test_1,
      fn_for_ds_test_2   = fn_for_ds_test_2,
      ds_test_methods    = ds_test_methods,
      test_alpha         = test_alpha,
      effect_template    = effect_template,
      center_by          = center_by,
      ...
    )
    results$objects$sim_power_list        <- sim_power_list
    results$objects$combined_power_effect <- sim_power_list$plot_data$power_by_effect
    
    results$files$power_by_effect <- file.path(
      results$root_dir, paste0(test_type, "_power_by_effect.pdf")
    )
    pdf(results$files$power_by_effect, width = 10, height = 7)
    plot_power_by_effect_size(
      combined_power = sim_power_list$plot_data$power_by_effect,
      fixed_n        = single_n,
      optimal_result = optimal_result,
      norm_config    = norm_config,
      eff_len        = eff_len
    )
    dev.off()
    cat("  Saved:", results$files$power_by_effect, "\n")
    
  } else {
    cat("\n5. POWER VS EFFECT SIZE  [skipped]\n")
  }
  
  #' ==========================================================================
  #' PHASE 6: SAVE
  #' Always runs when any other phase ran.
  #' ==========================================================================
  cat("\n6. SAVING RESULTS\n")
  
  rdata_results <- file.path(results$root_dir,paste0(test_type, "_results.RData"))
  save(results, file = rdata_results)
  
  rdata_ws <- file.path(results$root_dir, paste0(test_type, "_workspace.RData"))
  save.image(rdata_ws)
  
  cat("  Saved:", rdata_results, "\n")
  cat("  Saved:", rdata_ws, "\n")
  
  #' ==========================================================================
  #' FINALIZE
  #' ==========================================================================
  end_time   <- Sys.time()
  total_time <- difftime(end_time, start_time, units = "mins")
  results$params$timestamp_end         <- end_time
  results$params$total_runtime_minutes <- as.numeric(total_time)
  
  cat("\n", strrep("=", 60), "\n", sep = "")
  cat("SIMULATION COMPLETE (", method_display, " Approach)\n", sep = "")
  cat("  Phases run    :", paste(sort(phases), collapse = ", "), "\n")
  if (!is.null(results$objects$param_star)) {
    cat("  Threshold used:", round(results$objects$param_star, 4), "\n")
  }
  cat("  Total runtime :", round(total_time, 1), "minutes\n")
  cat(strrep("=", 60), "\n\n")
  
  invisible(results)
}


#' =============================================================================
#' SECTION 0 — SHARED SIMULATION PARAMETERS
#' =============================================================================
#' -- Simulation scale ---------------------------------------------------------
Nsim       <- 1e4
N_tradeoff <- 1e5

#' -- Distributions ------------------------------------------------------------
distributions <- c("exponential", "normal")

#' -- Threshold grids ----------------------------------------------------------
threshold_grid_classical <- seq(from = 0.00, to = 1, by = 0.01)
threshold_grid_ml        <- seq(from = 0.00, to = 1, by = 0.01)

#' -- Optimality tolerances ----------------------------------------------------
tol_pos  <- 0.005
loss_tol <- 0.01

#' -- Test settings ------------------------------------------------------------
test_alpha  <- 0.05
center_by   <- "median"
effect_size <- 0.5

#' -- Sample sizes -------------------------------------------------------------
sample_sizes      <- c(10, 20, 30, 40, 50)
single_n          <- 10
effect_sizes_plot <- c(0.0, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0)
sig_levels        <- seq(0.00, 1, by = 0.005)

#' -- Downstream test method labels --------------------------------------------
ds_test_methods <- c("test_1", "test_2", "adaptive")

#' -- Classical normality tests (Phase 1 ROC) ----------------------------------
norm_test      <- c("SW", "SF", "LF", "KS", "JB", "SKEW", "DAP", "AD", "CVM")
selected_tests <- c("SW", "SF", "LF", "JB", "SKEW", "DAP", "AD", "CVM")

#' =============================================================================
#' SECTION 1 — NORMALITY PRETEST CONFIGURATIONS
#' =============================================================================

#' -- 1a. Classical: Shapiro-Wilk ----------------------------------------------
norm_config_sw <- list(
  method = "classical",
  config = list(norm_test = "SW")
)

#' -- 1b. Classical: Anderson-Darling ------------------------------------------
norm_config_ad <- list(
  method = "classical",
  config = list(norm_test = "AD")
)

#' -- 1c. Classical: Shapiro-Francia -------------------------------------------
norm_config_sf <- list(
  method = "classical",
  config = list(norm_test = "SF")
)

#' -- 1d. Classical: Skewness --------------------------------------------------
norm_config_skewness <- list(
  method = "classical",
  config = list(norm_test = "SKEW")
)

#' -- 1e. Custom: Fisher combined SW + AD --------------------------------------
fisher_combined <- function(x) {
  p_sw <- shapiro.test(x)$p.value
  p_ad <- nortest::ad.test(x)$p.value
  
  stat <- -2 * sum(log(c(p_sw, p_ad)))
  
  c(p.value = pchisq(stat, df = 4, lower.tail = FALSE))
}

norm_config_fisher <- list(
  method = "custom",
  config = list(fn = fisher_combined)
)


#' =============================================================================
#' SECTION 1F — LOAD TRAINED SELECTED-FEATURE ML MODELS
#' =============================================================================
#' The user-framework does not train ML models.
#' It only loads models that were already trained in the ML framework.
#'
#' This file should come from the selected-feature ML run, not the full-feature
#' VIP-discovery run.
#' =============================================================================

load("results_selected_feat/trained_models.RData")


#' =============================================================================
#' SECTION 1G — ML SCORE FUNCTIONS USED BY THE USER-FRAMEWORK
#' =============================================================================

#' -----------------------------------------------------------------------------
#' ML score function for ROC curves.
#'
#' Returns P(Normal) for each requested ML model.
#' Smaller values mean stronger evidence against normality.
#'
#' If use_majority_vote = TRUE, an additional MajorityVote score is returned:
#'   MajorityVote = 1 - average P(Non_Normal)
#' -----------------------------------------------------------------------------

ml_fn_roc <- function(x,
                      trained_models,
                      model_name         = c("RF"),
                      use_majority_vote  = FALSE,
                      decision_threshold = 0.50,
                      custom_sample_size = NULL) {
  
  result <- classify_sample_prob(
    sample_data        = x,
    trained_models     = trained_models,
    model_name         = model_name,
    use_majority_vote  = use_majority_vote,
    decision_threshold = decision_threshold,
    custom_sample_size = custom_sample_size
  )
  
  per_model <- sapply(model_name, function(m) {
    as.numeric(result[[paste0(m, "_Prob_Normal")]])
  })
  
  if (isTRUE(use_majority_vote) &&
      length(model_name) > 1L &&
      "Proportion_Non_Normal" %in% names(result)) {
    
    per_model <- c(
      per_model,
      MajorityVote = as.numeric(1 - result$Proportion_Non_Normal)
    )
  }
  
  per_model
}


#' -----------------------------------------------------------------------------
#' ML score function for adaptive pretest decisions.
#'
#' Returns one scalar: P(Normal).
#'
#' The adaptive rule later uses:
#'   use parametric test if P(Normal) > threshold
#'   otherwise use the nonparametric test
#' -----------------------------------------------------------------------------

ml_fn_decision <- function(x,
                           trained_models,
                           model_name         = "RF",
                           use_majority_vote  = FALSE,
                           decision_threshold = 0.50,
                           custom_sample_size = NULL) {
  
  ## Fast path: single model, no majority vote — the normal Phases 2-6 case.
  ## Bypasses the classify_sample_prob() data.frame machinery entirely.
  if (!isTRUE(use_majority_vote) && length(model_name) == 1L) {
    
    ## Retrieve the model bundle for this sample's size (or nearest available).
    model_bundle <- get_ml_model(
      trained_models = trained_models,
      actual_size    = length(x),
      requested_size = custom_sample_size
    )
    
    ## Compute features and predict P(Normal) in one lean call.
    prob_normal <- ml_env$predict_prob_normal_fast(
      x            = x,
      model_bundle = model_bundle,
      model_name   = model_name
    )
    
    ## Return as a named scalar consistent with the full-path output format.
    return(setNames(prob_normal, "Sample1"))
  }
  
  ## Full path: majority vote or multi-model call — use classify_sample_prob().
  result <- classify_sample_prob(
    sample_data        = x,
    trained_models     = trained_models,
    model_name         = model_name,
    use_majority_vote  = use_majority_vote,
    decision_threshold = decision_threshold,
    custom_sample_size = custom_sample_size
  )
  
  if (isTRUE(use_majority_vote) &&
      length(model_name) > 1L &&
      "Proportion_Non_Normal" %in% names(result)) {
    
    out <- 1 - result$Proportion_Non_Normal
    
  } else {
    
    out <- result[[paste0(model_name[1], "_Prob_Normal")]]
  }
  
  setNames(as.numeric(out), result$Sample_Name)
}


#' =============================================================================
#' SECTION 1H — ML NORMALITY CONFIGURATIONS
#' =============================================================================

#' -- ML ROC configuration ------------------------------------------------------
#' Use this only for Phase 1, where you want to compare several ML models.
#' It returns multiple curves: RF, GBM, ANN, SVM, LR, and MajorityVote.
norm_config_ml_roc <- list(
  method = "custom",
  config = list(
    fn                 = ml_fn_roc,
    trained_models     = trained_models,
    model_name         = c("RF", "GBM", "ANN", "SVM"),#, "LR", "KNN"),
    use_majority_vote  = TRUE,
    decision_threshold = 0.50
  )
)

#' -- ML adaptive-decision configuration ----------------------------------------
#' Use this for Phases 2-6 after choosing the model you want for the adaptive rule.
#' Here RF is used as the normality pretest.
norm_config_ml <- list(
  method = "custom",
  config = list(
    fn                 = ml_fn_decision,
    trained_models     = trained_models,
    model_name         = "SVM",
    use_majority_vote  = FALSE,
    decision_threshold = 0.50
  )
)

#' =============================================================================
#' SECTION 2 — ONE-SAMPLE TESTS
#'
#' test_1 : one-sample t-test  (parametric)
#' test_2 : sign test          (non-parametric)
#' normality object: raw data vector
#' =============================================================================

onesample_fns <- list(
  gen_data           = onesample_data,
  get_parameters     = onesample_parameters,
  fn_to_get_norm_obj = raw_data,
  fn_for_ds_test_1   = one_sample_t_test,
  fn_for_ds_test_2   = sign_test
)

onesample_params <- list(
  Nsim              = Nsim,
  N_tradeoff        = N_tradeoff,
  distributions     = distributions,
  threshold_grid    = threshold_grid_classical,
  tol_pos           = tol_pos,
  loss_tol          = loss_tol,
  test_alpha        = test_alpha,
  center_by         = center_by,
  effect_size       = effect_size,
  sample_sizes      = sample_sizes,
  single_n          = single_n,
  effect_sizes_plot = effect_sizes_plot,
  sig_levels        = sig_levels,
  ds_test_methods   = ds_test_methods,
  norm_test         = norm_test,
  selected_tests    = selected_tests
)

#' ----------------------------------------------------------------------------
#' 2A. One-sample t-test vs sign test
#' ----------------------------------------------------------------------------

#' -- Classical SW -------------------------------------------------------------
res_os_ttest_sign_sw <- do.call(run_simulation, c(
  onesample_params,
  list(
    test_type = "onesample_ttest_vs_sign_sw_new_tradeoff_fn",
    norm_config = norm_config_sw,
    phases      = c(1, 2, 3, 4, 5, 6)),
  onesample_fns
))


#' -- Classical AD -------------------------------------------------------------
res_os_ttest_sign_ad <- do.call(run_simulation, c(
  onesample_params,
  list(
    test_type   = "onesample_ttest_vs_sign_ad",
    norm_config = norm_config_ad,
    phases      = c(2, 3, 4, 5, 6)),
  onesample_fns
))

#' -- Classical SW (per-n optimal thresholds) ----------------------------------
res_os_ttest_sign_sw_pern <- do.call(run_simulation, c(
  onesample_params,
  list(
    test_type        = "onesample_ttest_vs_sign_sw_pern_tradeoff",
    norm_config      = norm_config_sw,
    per_n_thresholds = TRUE,
    phases           = c(2, 3, 4, 5, 6)),
  onesample_fns
))

#' -- Classical SF -------------------------------------------------------------
res_os_ttest_sign_sf <- do.call(run_simulation, c(
  onesample_params,
  list(
    test_type   = "onesample_ttest_vs_sign_sf",
    norm_config = norm_config_sf,
    phases      = c(2, 3, 4, 5, 6)),
  onesample_fns
))


#' -- Custom Fisher ------------------------------------------------------------
res_os_ttest_sign_fisher <- do.call(run_simulation, c(
  onesample_params,
  list(
    test_type   = "onesample_ttest_vs_sign_fisher",
    norm_config = norm_config_fisher,
    phases      = c(1, 2, 3, 4, 5, 6)),
  onesample_fns
))


#' -- ML: Phase 1 (all models compared) ---------------------------------------
res_os_ttest_sign_ml_roc <- do.call(run_simulation, c(
  modifyList(onesample_params, list(threshold_grid = threshold_grid_ml)),
  list(
    test_type   = "onesample_ttest_vs_sign_t1",
    norm_config = norm_config_ml_roc,
    phases      = c(1)),
  onesample_fns
))

#' -- ML: Phases 2-6 (chosen model after reviewing Phase 1) -------------------
res_os_ttest_sign_ml <- do.call(run_simulation, c(
  modifyList(onesample_params, list(threshold_grid = threshold_grid_ml)),
  list(
    test_type   = "onesample_ttest_vs_sign_t1",
    norm_config = norm_config_ml,
    per_n_thresholds = TRUE,
    phases      = c(2, 3, 4, 5, 6)),
  onesample_fns
))



#' ----------------------------------------------------------------------------
#' 2B. One-sample t-test vs permutation test
#' ----------------------------------------------------------------------------

onesample_perm_fns <- modifyList(
  onesample_fns,
  list(fn_for_ds_test_2 = one_sample_perm_test)
)

onesample_perm_params <- modifyList(onesample_params, list(
  distributions = c("laplace", "normal"),
  center_by     = "mean"
))

#' -- Classical SW -------------------------------------------------------------
res_os_ttest_perm_sw <- do.call(run_simulation, c(
  onesample_perm_params,
  list(
    test_type = "onesample_ttest_vs_perm_sw",
    norm_config = norm_config_sw,
    phases      = c(1, 2, 3, 4, 5, 6)),
  onesample_perm_fns
))

#' -- Classical SF -------------------------------------------------------------
res_os_ttest_perm_sf <- do.call(run_simulation, c(
  onesample_perm_params,
  list(
    test_type   = "onesample_ttest_vs_perm_sf",
    norm_config = norm_config_sf,
    phases      = c(2, 3, 4, 5, 6)),
  onesample_perm_fns
))

#' -- Classical AD -------------------------------------------------------------
res_os_ttest_perm_ad <- do.call(run_simulation, c(
  onesample_perm_params,
  list(
    test_type   = "onesample_ttest_vs_perm_ad",
    norm_config = norm_config_ad,
    phases      = c(2, 3, 4, 5, 6)),
  onesample_perm_fns
))

#' -- Custom Fisher ------------------------------------------------------------
res_os_ttest_perm_fisher <- do.call(run_simulation, c(
  onesample_perm_params,
  list(
    test_type   = "onesample_ttest_vs_perm_fisher",
    norm_config = norm_config_fisher,
    phases      = c(1, 2, 3, 4, 5, 6)),
  onesample_perm_fns
))

#' -- ML: Phase 1 (all models compared) ---------------------------------------
res_os_ttest_perm_ml_roc <- do.call(run_simulation, c(
  modifyList(onesample_perm_params, list(threshold_grid = threshold_grid_ml)),
  list(
    test_type   = "onesample_ttest_vs_perm_ml",
    norm_config = norm_config_ml_roc,
    phases      = c(1)),
  onesample_perm_fns
))


#' -- ML: Phases 2-6 (chosen model after reviewing Phase 1) -------------------
res_os_ttest_perm_ml <- do.call(run_simulation, c(
  modifyList(onesample_perm_params, list(threshold_grid = threshold_grid_ml)),
  list(
    test_type   = "onesample_ttest_vs_perm_ml",
    norm_config = norm_config_ml,
    phases      = c(2, 3, 4, 5, 6)),
  onesample_perm_fns
))

#' -- Classical SW (per-n optimal thresholds) ----------------------------------
res_os_ttest_perm_sw_pern <- do.call(run_simulation, c(
  onesample_perm_params,
  list(
    test_type        = "onesample_ttest_vs_perm_sw_pern",
    norm_config      = norm_config_sw,
    per_n_thresholds = TRUE,
    phases           = c(2, 3, 4, 5, 6)),
  onesample_perm_fns
))

#' =============================================================================
#' SECTION 3 — TWO-SAMPLE TESTS
#'
#' test_1 : two-sample t-test   (parametric)
#' test_2 : Mann-Whitney U test (non-parametric)
#' normality object: raw data frame (group + value columns)
#' =============================================================================

twosample_fns <- list(
  gen_data           = two_sample_data,
  get_parameters     = twosample_parameters,
  fn_to_get_norm_obj = raw_data,
  fn_for_ds_test_1   = twosample_t_test,
  fn_for_ds_test_2   = Mann_whitney_U_test
)

twosample_params <- modifyList(onesample_params, list(
  center_by = "median"
))

#' ----------------------------------------------------------------------------
#' 3A. Two-sample t-test vs Mann-Whitney U test
#' ----------------------------------------------------------------------------
#' 
#' -- Classical SW -------------------------------------------------------------
res_ts_ttest_mw_sw <- do.call(run_simulation, c(
  twosample_params,
  list(
    test_type   = "twosample_ttest_vs_mw_sw",
    norm_config = norm_config_sw,
    phases      = c(1, 2, 3, 4, 5, 6)),
  twosample_fns
))


#' -- Classical SF -------------------------------------------------------------
res_ts_ttest_mw_sf <- do.call(run_simulation, c(
  twosample_params,
  list(
    test_type   = "twosample_ttest_vs_mw_sf",
    norm_config = norm_config_sf,
    phases      = c(2, 3, 4, 5, 6)),
  twosample_fns
))

#' -- Classical AD -------------------------------------------------------------
res_ts_ttest_mw_ad <- do.call(run_simulation, c(
  twosample_params,
  list(
    test_type   = "twosample_ttest_vs_mw_ad",
    norm_config = norm_config_ad,
    phases      = c(2, 3, 4, 5, 6)),
  twosample_fns
))


#' -- Custom Fisher ------------------------------------------------------------
res_ts_ttest_mw_fisher <- do.call(run_simulation, c(
  twosample_params,
  list(
    test_type   = "twosample_ttest_vs_mw_fisher",
    norm_config = norm_config_fisher,
    phases      = c(1, 2, 3, 4, 5, 6)),
  twosample_fns
))

#' -- ML: Phase 1 (all models compared) ---------------------------------------
res_ts_ttest_mw_ml_roc <- do.call(run_simulation, c(
  modifyList(twosample_params, list(threshold_grid = threshold_grid_ml)),
  list(
    test_type   = "twosample_ttest_vs_mw_ml",
    norm_config = norm_config_ml_roc,
    phases      = c(1)),
  twosample_fns
))

#' -- ML: Phases 2-6 (chosen model after reviewing Phase 1) -------------------
res_ts_ttest_mw_ml <- do.call(run_simulation, c(
  modifyList(twosample_params, list(threshold_grid = threshold_grid_ml)),
  list(
    test_type   = "twosample_ttest_vs_mw_ml",
    norm_config = norm_config_ml,
    phases      = c(2, 3, 4, 5, 6)),
  twosample_fns
))


#' -- Classical SW (per-n optimal thresholds) ----------------------------------
res_ts_ttest_mw_sw <- do.call(run_simulation, c(
  twosample_params,
  list(
    test_type   = "twosample_ttest_vs_mw_sw_pern",
    norm_config = norm_config_sw,
    per_n_thresholds = TRUE,
    phases      = c(2, 3, 4, 5, 6)),
  twosample_fns
))


#' ----------------------------------------------------------------------------
#' 3B. Two-sample t-test vs permutation test
#' ----------------------------------------------------------------------------

twosample_perm_fns <- modifyList(
  twosample_fns,
  list(fn_for_ds_test_2 = two_sample_perm_test)
)

twosample_perm_params <- modifyList(onesample_params, list(
  center_by = "mean"
))

#' -- Classical SW -------------------------------------------------------------
res_ts_ttest_perm_sw <- do.call(run_simulation, c(
  twosample_perm_params,
  list(
    test_type   = "twosample_ttest_vs_perm_sw",
    norm_config = norm_config_sw,
    phases      = c(2, 3, 4, 5, 6)),
  twosample_perm_fns
))


#' -- Classical SF -------------------------------------------------------------
res_ts_ttest_perm_sf <- do.call(run_simulation, c(
  twosample_perm_params,
  list(
    test_type   = "twosample_ttest_vs_perm_sf",
    norm_config = norm_config_sf,
    phases      = c(2, 3, 4, 5, 6)),
  twosample_perm_fns
))

#' -- Classical AD -------------------------------------------------------------
res_ts_ttest_mw_ad <- do.call(run_simulation, c(
  twosample_perm_params,
  list(
    test_type   = "twosample_ttest_vs_perm_ad",
    norm_config = norm_config_ad,
    phases      = c(2, 3, 4, 5, 6)),
  twosample_perm_fns
))

#' -- Custom Fisher ------------------------------------------------------------
res_ts_ttest_perm_fisher <- do.call(run_simulation, c(
  twosample_perm_params,
  list(
    test_type   = "twosample_ttest_vs_perm_fisher",
    norm_config = norm_config_fisher,
    phases      = c(1, 2, 3, 4, 5, 6)),
  twosample_perm_fns
))

#' -- ML: Phase 1 (all models compared) ---------------------------------------
res_ts_ttest_perm_ml_roc <- do.call(run_simulation, c(
  modifyList(twosample_perm_params, list(threshold_grid = threshold_grid_ml)),
  list(
    test_type   = "twosample_ttest_vs_perm_ml",
    norm_config = norm_config_ml_roc,
    phases      = c(1)),
  twosample_perm_fns
))

#' -- ML: Phases 2-6 (chosen model after reviewing Phase 1) -------------------
res_ts_ttest_perm_ml <- do.call(run_simulation, c(
  modifyList(twosample_perm_params, list(threshold_grid = threshold_grid_ml)),
  list(
    test_type   = "twosample_ttest_vs_perm_ml",
    norm_config = norm_config_ml,
    phases      = c(2, 3, 4, 5, 6)),
  twosample_perm_fns
))


#' =============================================================================
#' SECTION 4 — ANOVA
#'
#' test_1 : one-way ANOVA       (parametric)
#' test_2 : Kruskal-Wallis test (non-parametric)
#' normality object: ANOVA residuals
#'
#' effect_size: list(H0 = c(0,0,0), H1 = c(0,0,0.5))
#'   H0 = all zeros (no group differences)
#'   H1 = last group shifted
#' =============================================================================

anova_fns <- list(
  gen_data           = anova_gen_data,
  get_parameters     = anova_parameters,
  fn_to_get_norm_obj = anova_residuals,
  fn_for_ds_test_1   = one_way_anova,
  fn_for_ds_test_2   = kruskal_wallis_test
)

anova_params <- modifyList(onesample_params, list(
  effect_size = list(H0 = c(0, 0, 0), H1 = c(0, 0, 0.5))
))

#' ----------------------------------------------------------------------------
#' 4A. One-way ANOVA vs Kruskal-Wallis test
#' ----------------------------------------------------------------------------

#' -- Classical SW -------------------------------------------------------
res_anova_kw_sw <- do.call(run_simulation, c(
  anova_params,
  list(
    test_type   = "anova_vs_kw_sw",
    norm_config = norm_config_sw,
    phases      = c(1, 2, 3, 4, 5, 6)),
  anova_fns
))

#' -- Classical SF -------------------------------------------------------
res_anova_kw_sf <- do.call(run_simulation, c(
  anova_params,
  list(
    test_type   = "anova_vs_kw_sf",
    norm_config = norm_config_sf,
    phases      = c(1, 2, 3, 4, 5, 6)),
  anova_fns
))


#' -- Classical Skewness -------------------------------------------------------
res_anova_kw_skew <- do.call(run_simulation, c(
  anova_params,
  list(
    test_type   = "anova_vs_kw_skewness",
    norm_config = norm_config_skewness,
    phases      = c(1, 2, 3, 4, 5, 6)),
  anova_fns
))

#' -- Custom Fisher ------------------------------------------------------------
res_anova_kw_fisher <- do.call(run_simulation, c(
  anova_params,
  list(
    test_type   = "anova_vs_kw_fisher",
    norm_config = norm_config_fisher,
    phases      = c(1, 2, 3, 4, 5, 6)),
  anova_fns
))

#' -- ML: Phase 1 (all models compared) ---------------------------------------
res_anova_kw_ml_roc <- do.call(run_simulation, c(
  modifyList(anova_params, list(threshold_grid = threshold_grid_ml)),
  list(
    test_type   = "anova_vs_kw_ml",
    norm_config = norm_config_ml_roc,
    phases      = c(1)),
  anova_fns
))

#' -- ML: Phases 2-6 (chosen model after reviewing Phase 1) -------------------
res_anova_kw_ml <- do.call(run_simulation, c(
  modifyList(anova_params, list(threshold_grid = threshold_grid_ml)),
  list(
    test_type   = "anova_vs_kw_ml",
    norm_config = norm_config_ml,
    phases      = c(2, 3, 4, 5, 6)),
  anova_fns
))

#' -- Classical SW (per-n optimal thresholds) ----------------------------------
res_anova_kw_sw_pern <- do.call(run_simulation, c(
  anova_params,
  list(
    test_type   = "anova_vs_kw_sw_pern",
    norm_config = norm_config_sw,
    per_n_thresholds = TRUE,
    phases      = c(1, 2, 3, 4, 5, 6)),
  anova_fns
))

#' ----------------------------------------------------------------------------
#' 4B. One-way ANOVA vs permutation ANOVA
#' ----------------------------------------------------------------------------

anova_perm_fns <- modifyList(
  anova_fns,
  list(fn_for_ds_test_2 = permutation_anova)
)

#' -- Classical SW -------------------------------------------------------
res_anova_perm_sw <- do.call(run_simulation, c(
  anova_params,
  list(
    test_type   = "anova_vs_perm_sw",
    norm_config = norm_config_sw,
    phases      = c(1, 2, 3, 4, 5, 6)),
  anova_perm_fns
))

#' -- Classical SF -------------------------------------------------------
res_anova_perm_sf <- do.call(run_simulation, c(
  anova_params,
  list(
    test_type   = "anova_vs_perm_sf",
    norm_config = norm_config_sf,
    phases      = c(1, 2, 3, 4, 5, 6)),
  anova_perm_fns
))


#' -- Classical Skewness -------------------------------------------------------
res_anova_perm_skew <- do.call(run_simulation, c(
  anova_params,
  list(
    test_type   = "anova_vs_perm_skewness",
    norm_config = norm_config_skewness,
    phases      = c(1, 2, 3, 4, 5, 6)),
  anova_perm_fns
))

#' -- Custom Fisher ------------------------------------------------------------
res_anova_perm_fisher <- do.call(run_simulation, c(
  anova_params,
  list(
    test_type   = "anova_vs_perm_fisher",
    norm_config = norm_config_fisher,
    phases      = c(1, 2, 3, 4, 5, 6)),
  anova_perm_fns
))

#' -- ML: Phase 1 (all models compared) ---------------------------------------
res_anova_perm_ml_roc <- do.call(run_simulation, c(
  modifyList(anova_params, list(threshold_grid = threshold_grid_ml)),
  list(
    test_type   = "anova_vs_perm_ml",
    norm_config = norm_config_ml_roc,
    phases      = c(1)),
  anova_perm_fns
))

#' -- ML: Phases 2-6 (chosen model after reviewing Phase 1) -------------------
res_anova_perm_ml <- do.call(run_simulation, c(
  modifyList(anova_params, list(threshold_grid = threshold_grid_ml)),
  list(
    test_type   = "anova_vs_perm_ml",
    norm_config = norm_config_ml,
    phases      = c(2, 3, 4, 5, 6)),
  anova_perm_fns
))


#' =============================================================================
#' SECTION 5 — REGRESSION
#'
#' test_1 : OLS simple linear regression (parametric)
#' test_2 : rank-based regression        (non-parametric)
#' normality object: OLS residuals
#'
#' effect_size maps to beta1 (slope).
#'   H0: beta1 = 0, H1: beta1 = effect_size
#' =============================================================================

regression_fns <- list(
  gen_data           = reg_data,
  get_parameters     = reg_parameters,
  fn_to_get_norm_obj = reg_residuals,
  fn_for_ds_test_1   = simple_linear_reg,
  fn_for_ds_test_2   = rank_regression
)

regression_params <- modifyList(onesample_params, list(
  center_by = "mean"
))

#' ----------------------------------------------------------------------------
#' 5A. OLS regression vs rank-based regression
#' ----------------------------------------------------------------------------

#' -- Classical SW -------------------------------------------------------------
res_reg_rank_sw <- do.call(run_simulation, c(
  regression_params,
  list(
    test_type   = "regression_ols_vs_rank_sw",
    norm_config = norm_config_sw,
    phases      = c(1, 2, 3, 4, 5, 6)),
  regression_fns
))

#' -- Classical AD -------------------------------------------------------------
res_reg_rank_ad <- do.call(run_simulation, c(
  regression_params,
  list(
    test_type   = "regression_ols_vs_rank_ad",
    norm_config = norm_config_ad,
    phases      = c(1, 2, 3, 4, 5, 6)),
  regression_fns
))

#' -- Custom Fisher ------------------------------------------------------------
res_reg_rank_fisher <- do.call(run_simulation, c(
  regression_params,
  list(
    test_type   = "regression_ols_vs_rank_fisher",
    norm_config = norm_config_fisher,
    phases      = c(1, 2, 3, 4, 5, 6)),
  regression_fns
))

#' -- ML: Phase 1 (all models compared) ---------------------------------------
res_reg_rank_ml_roc <- do.call(run_simulation, c(
  modifyList(regression_params, list(threshold_grid = threshold_grid_ml)),
  list(
    test_type   = "regression_ols_vs_rank_ml",
    norm_config = norm_config_ml_roc,
    phases      = c(1)),
  regression_fns
))

#' -- ML: Phases 2-6 (chosen model after reviewing Phase 1) -------------------
res_reg_rank_ml <- do.call(run_simulation, c(
  modifyList(regression_params, list(threshold_grid = threshold_grid_ml)),
  list(
    test_type   = "regression_ols_vs_rank_ml",
    norm_config = norm_config_ml,
    phases      = c(2, 3, 4, 5, 6)),
  regression_fns
))

#' -- Classical SW (per-n optimal thresholds) ----------------------------------
res_reg_rank_sw_pern <- do.call(run_simulation, c(
  regression_params,
  list(
    test_type   = "regression_ols_vs_rank_sw_pern",
    norm_config = norm_config_sw,
    per_n_thresholds = TRUE,
    phases      = c(1, 2, 3, 4, 5, 6)),
  regression_fns
))


#' ----------------------------------------------------------------------------
#' 5B. OLS regression vs permutation regression
#' ----------------------------------------------------------------------------

regression_perm_fns <- modifyList(
  regression_fns,
  list(fn_for_ds_test_2 = perm_regression)
)

#' -- Classical SW -------------------------------------------------------------
res_reg_perm_sw <- do.call(run_simulation, c(
  regression_params,
  list(
    test_type   = "regression_ols_vs_perm_sw",
    norm_config = norm_config_sw,
    phases      = c(1, 2, 3, 4, 5, 6)),
  regression_perm_fns
))

#' -- Custom Fisher ------------------------------------------------------------
res_reg_perm_fisher <- do.call(run_simulation, c(
  regression_params,
  list(
    test_type   = "regression_ols_vs_perm_fisher",
    norm_config = norm_config_fisher,
    phases      = c(1, 2, 3, 4, 5, 6)),
  regression_perm_fns
))

#' -- ML: Phase 1 (all models compared) ---------------------------------------
res_reg_perm_ml_roc <- do.call(run_simulation, c(
  modifyList(regression_params, list(threshold_grid = threshold_grid_ml)),
  list(
    test_type   = "regression_ols_vs_perm_ml",
    norm_config = norm_config_ml_roc,
    phases      = c(1)),
  regression_perm_fns
))

#' -- ML: Phases 2-6 (chosen model after reviewing Phase 1) -------------------
res_reg_perm_ml <- do.call(run_simulation, c(
  modifyList(regression_params, list(threshold_grid = threshold_grid_ml)),
  list(
    test_type   = "regression_ols_vs_perm_ml",
    norm_config = norm_config_ml,
    phases      = c(2, 3, 4, 5, 6)),
  regression_perm_fns
))


#' =============================================================================
#' EXPLORATION OF THE SELECTED PRETEST THRESHOLDS ON THE NORMALITY ROC CURVES
#'  ( Shapiro-Wilk per-n   vs   SVM per-n )
#' -----------------------------------------------------------------------------

#' =============================================================================
#' SECTION 0 — LOAD SAVED RESULTS AND EXTRACT THE SELECTED THRESHOLDS
#' =============================================================================

results_dir <- file.path(
  "/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation",
  "User_framework_Rpkg/ML_application/SW_vs_ML_Application_eploration/"
)

sw_env  <- new.env()
svm_env <- new.env()

load(file.path(results_dir,
               "onesample_ttest_vs_sign_sw_pern_tradeoff",
               "onesample_ttest_vs_sign_sw_pern_tradeoff_results.RData"),
     envir = sw_env)

load(file.path(results_dir,
               "onesample_ttest_vs_sign_t1",
               "onesample_ttest_vs_sign_t1_results.RData"),
     envir = svm_env)

#' Per-n optimal thresholds 
sw_thr  <- sw_env$results$objects$threshold_by_n
svm_thr <- svm_env$results$objects$threshold_by_n

#' The threshold 
grid_used <- svm_env$results$params$threshold_grid

#' checks (print to confirm before running the analysis)
sw_env$results$params$norm_config   
sw_thr                              # SW p-value cutoffs
svm_thr                             # SVM P(Normal) cutoffs
grid_used                           


#' =============================================================================
#' SECTION 1 — HELPER FUNCTIONS
#' =============================================================================

#' -----------------------------------------------------------------------------
#' Map ONE threshold to its (FPR, TPR) operating point on a ROC grid.
#' Returns the nearest grid node (an exact match when the threshold is a
#' multiple of the grid step, which it is by construction).
#' -----------------------------------------------------------------------------
roc_operating_point <- function(grid, fpr, tpr, threshold) {
  if (is.na(threshold)) {
    return(c(threshold = NA, grid_used = NA, FPR = NA, TPR = NA))
  }
  j <- which.min(abs(grid - threshold))
  c(threshold = threshold,
    grid_used = grid[j],
    FPR       = fpr[j],
    TPR       = tpr[j])
}


#' -----------------------------------------------------------------------------
#' Compute one method's normality-detection ROC at a single sample size.
#' Returns a long data frame: grid, FPR, TPR, curve.
#' -----------------------------------------------------------------------------
norm_roc_frame <- function(norm_config, n, threshold_grid, H1_dist, Nsim,
                           gen_data, get_parameters, fn_to_get_norm_obj,
                           center_by = "median", ...) {
  
  method <- tolower(norm_config$method)
  
  if (method == "classical") {
    
    test_code <- norm_config$config$norm_test
    roc <- fn_for_norm_test_roc_curve(
      n                  = n,
      alpha_pretest      = threshold_grid,
      H1_dist            = H1_dist,
      tests              = test_code,
      Nsim               = Nsim,
      gen_data           = gen_data,
      get_parameters     = get_parameters,
      fn_to_get_norm_obj = fn_to_get_norm_obj,
      center_by          = center_by,
      ...
    )
    data.frame(
      grid  = roc$alpha,
      FPR   = as.numeric(roc$FPR[test_code, ]),
      TPR   = as.numeric(roc$TPR[test_code, ]),
      curve = test_code,
      stringsAsFactors = FALSE
    )
    
  } else {  # custom (ML)
    
    fn          <- norm_config$config$fn
    custom_args <- norm_config$config[setdiff(names(norm_config$config), "fn")]
    roc <- fn_for_custom_test_roc_curve(
      n                  = n,
      threshold_grid     = threshold_grid,
      H1_dist            = H1_dist,
      Nsim               = Nsim,
      gen_data           = gen_data,
      get_parameters     = get_parameters,
      fn_to_get_norm_obj = fn_to_get_norm_obj,
      custom_fn          = fn,
      custom_args        = custom_args,
      center_by          = center_by,
      ...
    )
    cn <- rownames(roc$FPR)[1]   # single curve (e.g. "SVM")
    data.frame(
      grid  = roc$threshold,
      FPR   = as.numeric(roc$FPR[1, ]),
      TPR   = as.numeric(roc$TPR[1, ]),
      curve = cn,
      stringsAsFactors = FALSE
    )
  }
}


#' =============================================================================
#' SECTION 2 — MAIN ANALYSIS
#' Build overlaid SW-vs-SVM detection ROCs across sample sizes and locate each
#' method's selected per-n threshold on its own curve.
#' =============================================================================
explore_norm_roc_overlay <- function(
    sample_sizes,
    sw_threshold_by_n,
    svm_threshold_by_n,
    trained_models,
    svm_model_name     = "SVM",
    threshold_grid     = seq(0, 1, by = 0.025),
    H1_dist            = "exponential",
    Nsim               = 10000,
    gen_data           = onesample_data,
    get_parameters     = onesample_parameters,
    fn_to_get_norm_obj = raw_data,
    center_by          = "median",
    ...) {
  
  #' Pretest configurations -----------------------------------------------------
  sw_config  <- list(
    method = "classical",
    config = list(norm_test = "SW")
  )
  
  svm_config <- list(
    method = "custom",
    config = list(
      fn                 = ml_fn_roc,
      trained_models     = trained_models,
      model_name         = svm_model_name,
      use_majority_vote  = FALSE,
      decision_threshold = 0.50
    )
  )
  
  sw_label  <- "SW"
  svm_label <- sprintf("%s", svm_model_name)
  
  roc_all <- list()
  pts_all <- list()
  auc_all <- list()
  
  #' Loop over sample sizes -----------------------------------------------------
  for (n in sample_sizes) {
    key <- as.character(n)
    cat(sprintf("\n--- Sample size n = %d -------------------------------\n", n))
    
    cat("  SW  detection ROC...\n")
    sw_roc  <- norm_roc_frame(sw_config,  n, threshold_grid, H1_dist, Nsim,
                              gen_data, get_parameters, fn_to_get_norm_obj,
                              center_by, ...)
    
    cat("  SVM detection ROC...\n")
    svm_roc <- norm_roc_frame(svm_config, n, threshold_grid, H1_dist, Nsim,
                              gen_data, get_parameters, fn_to_get_norm_obj,
                              center_by, ...)
    
    sw_roc$method  <- sw_label
    sw_roc$n       <- n
    svm_roc$method <- svm_label
    svm_roc$n      <- n
    roc_all[[key]] <- rbind(sw_roc, svm_roc)
    
    #' Locate each selected threshold on its own ROC curve
    sw_t  <- suppressWarnings(as.numeric(sw_threshold_by_n[key]))
    svm_t <- suppressWarnings(as.numeric(svm_threshold_by_n[key]))
    
    sw_pt  <- roc_operating_point(sw_roc$grid,  sw_roc$FPR,  sw_roc$TPR,  sw_t)
    svm_pt <- roc_operating_point(svm_roc$grid, svm_roc$FPR, svm_roc$TPR, svm_t)
    
    pts_all[[key]] <- data.frame(
      n         = n,
      method    = c(sw_label, svm_label),
      threshold = c(sw_pt["threshold"], svm_pt["threshold"]),
      grid_used = c(sw_pt["grid_used"], svm_pt["grid_used"]),
      FPR       = c(sw_pt["FPR"], svm_pt["FPR"]),
      TPR       = c(sw_pt["TPR"], svm_pt["TPR"]),
      stringsAsFactors = FALSE, row.names = NULL
    )
    
    #' AUC for each full curve
    auc_all[[key]] <- data.frame(
      n      = n,
      method = c(sw_label, svm_label),
      AUC    = c(compute_auc(sw_roc$FPR,  sw_roc$TPR,  ensure_endpoints = TRUE),
                 compute_auc(svm_roc$FPR, svm_roc$TPR, ensure_endpoints = TRUE)),
      stringsAsFactors = FALSE, 
      row.names = NULL
    )
  }
  
  list(
    roc    = do.call(rbind, roc_all),
    points = do.call(rbind, pts_all),
    auc    = do.call(rbind, auc_all),
    meta   = list(
      sw_label       = sw_label,
      svm_label      = svm_label,
      H1_dist        = H1_dist,
      Nsim           = Nsim,
      threshold_grid = threshold_grid
    )
  )
}


#' =============================================================================
#' SECTION 3 — ROC Curves
#' =============================================================================
plot_norm_roc_overlay_base <- function(explore_out,
                                       file  = NULL,
                                       width = 11,
                                       height = 8,
                                       title = "Normality-detection ROC by sample size: SW vs SVM") {
  
  roc     <- explore_out$roc
  pts     <- explore_out$points
  auc     <- explore_out$auc
  ns      <- sort(unique(roc$n))
  methods <- c(explore_out$meta$sw_label, explore_out$meta$svm_label)
  cols    <- setNames(c("blue", "red"), methods)   
  
  np   <- length(ns)
  ncol <- ceiling(sqrt(np))
  nrow <- ceiling(np / ncol)
  
  if (!is.null(file)) pdf(file, width = width, height = height)
  op <- par(mfrow    = c(nrow, ncol),
            mar      = c(4, 4, 2.6, 1),
            oma      = c(0, 0, 3, 0),
            font.lab = 2)
  on.exit({ par(op); if (!is.null(file)) dev.off() }, add = TRUE)
  
  for (n in ns) {
    r <- roc[roc$n == n, ]
    p <- pts[pts$n == n, ]
    a <- auc[auc$n == n, ]
    
    plot(NA, xlim = c(0, 1), ylim = c(0, 1),
         xlab = "FPR  P(reject | normal)",
         ylab = "TPR  P(reject | non-normal)",
         main = sprintf("n = %d", n))
    abline(0, 1, lty = 2, col = "grey70")
    
    #' ROC curves
    for (m in methods) {
      rm  <- r[r$method == m, ]
      ord <- order(rm$FPR, rm$TPR)
      lines(rm$FPR[ord], rm$TPR[ord], col = cols[m], lwd = 2.6)
    }
    
    #' Selected-threshold operating points
    for (i in seq_len(nrow(p))) {
      if (is.na(p$FPR[i])) next
      points(p$FPR[i], p$TPR[i], pch = 23, bg = cols[p$method[i]], col = "black", cex = 1.7, lwd = 1.1)
      
      text(p$FPR[i], p$TPR[i], labels = sprintf("t*=%.3f", p$threshold[i]), 
           pos = 4, offset = 0.5, cex = 0.72, col = cols[p$method[i]], font = 2)
    }
    
    #' Per-panel legend: AUC + operating point
    leg <- vapply(methods, function(m) {
      am <- a[a$method == m, ]
      pm <- p[p$method == m, ]
      sprintf("%s: AUC=%.3f | t*=%.3f (FPR=%.2f, TPR=%.2f)",
              m, am$AUC, pm$threshold, pm$FPR, pm$TPR)
    }, character(1))
    legend("bottomright", legend = leg, col = cols[methods],
           lwd = 2.6, pch = 23, pt.bg = cols[methods], bty = "n", cex = 0.68)
  }
  
  mtext(title, side = 3, outer = TRUE, line = 0.8, cex = 1.25, font = 2)
  invisible(explore_out)
}


#' =============================================================================
#' SECTION 4 — SUMMARY TABLE
#' =============================================================================
build_wide_table <- function(explore_out, digits = 4) {
  pts     <- explore_out$points
  auc     <- explore_out$auc
  sw_lab  <- explore_out$meta$sw_label
  svm_lab <- explore_out$meta$svm_label
  
  sw_p  <- pts[pts$method == sw_lab,  c("n", "threshold", "FPR", "TPR")]
  svm_p <- pts[pts$method == svm_lab, c("n", "threshold", "FPR", "TPR")]
  names(sw_p)[-1]  <- paste0(c("threshold", "FPR", "TPR"), "_SW")
  names(svm_p)[-1] <- paste0(c("threshold", "FPR", "TPR"), "_SVM")
  
  sw_a  <- auc[auc$method == sw_lab,  c("n", "AUC")]; names(sw_a)[2]  <- "AUC_SW"
  svm_a <- auc[auc$method == svm_lab, c("n", "AUC")]; names(svm_a)[2] <- "AUC_SVM"
  
  out <- Reduce(function(a, b) merge(a, b, by = "n"),list(sw_p, svm_p, sw_a, svm_a))
  out <- out[order(out$n), ]
  out$better_detector <- ifelse(out$AUC_SVM > out$AUC_SW, "SVM", "SW")
  
  out <- out[, c("n",
                 "threshold_SW", "threshold_SVM",
                 "FPR_SW",        "FPR_SVM",
                 "TPR_SW",        "TPR_SVM",
                 "AUC_SW",        "AUC_SVM",
                 "better_detector")]
  
  num <- setdiff(names(out)[sapply(out, is.numeric)], "n")
  out[num] <- lapply(out[num], round, digits)
  rownames(out) <- NULL
  out
}


#' =============================================================================
#' SECTION 5 — RUN
#' =============================================================================
expl <- explore_norm_roc_overlay(
  sample_sizes       = c(10, 20, 30, 40, 50),
  sw_threshold_by_n  = sw_thr,
  svm_threshold_by_n = svm_thr,
  trained_models     = trained_models,
  svm_model_name     = "SVM",
  threshold_grid     = grid_used,        
  H1_dist            = "exponential",    
  Nsim               = 10000,             
  center_by          = "median"
)



#' print summary table
wide <- build_wide_table(expl)
print(wide, row.names = FALSE)
write.csv(wide, file.path(results_dir, "sw_vs_svm_roc_summary.csv"), row.names = FALSE)

#' ROC figure 
plot_norm_roc_overlay_base(
  expl,
  file = file.path(results_dir, "sw_vs_svm_roc_overlay.pdf")
)

# save results
save(expl, sw_thr, svm_thr, grid_used, file = "sw_vs_svm_results.RData")


#' =============================================================================
#' SECTION 6 — FIXED-THRESHOLD VARIANT (SW = 0.05, SVM = 0.50, for every n)
#' Instead of the per-n "optimal" thresholds used above, mark each curve at
#' the conventional fixed cutoffs: alpha = 0.05 for SW, decision boundary
#' = 0.50 for SVM. Same sample sizes / ROC curves, different operating points.
#' =============================================================================
fixed_sample_sizes <- c(10, 20, 30, 40, 50)

sw_thr_fixed  <- setNames(rep(0.05, length(fixed_sample_sizes)), as.character(fixed_sample_sizes))
svm_thr_fixed <- setNames(rep(0.50, length(fixed_sample_sizes)), as.character(fixed_sample_sizes))

expl_fixed <- explore_norm_roc_overlay(
  sample_sizes       = fixed_sample_sizes,
  sw_threshold_by_n  = sw_thr_fixed,
  svm_threshold_by_n = svm_thr_fixed,
  trained_models     = trained_models,
  svm_model_name     = "SVM",
  threshold_grid     = grid_used,
  H1_dist            = "exponential",
  Nsim               = 10000,
  center_by          = "median"
)

#' print summary table (fixed-threshold variant)
wide_fixed <- build_wide_table(expl_fixed)
print(wide_fixed, row.names = FALSE)
write.csv(wide_fixed, file.path(results_dir, "sw_vs_svm_roc_summary_fixed.csv"), row.names = FALSE)

#' ROC figure (fixed-threshold variant)
plot_norm_roc_overlay_base(
  expl_fixed,
  file  = file.path(results_dir, "sw_vs_svm_roc_overlay_fixed.pdf"),
  title = "Normality-detection ROC by sample size: SW (alpha=0.05) vs SVM (t=0.50)"
)

# save results (fixed-threshold variant)
save(expl_fixed, sw_thr_fixed, svm_thr_fixed, grid_used, file = "sw_vs_svm_results_fixed.RData")
