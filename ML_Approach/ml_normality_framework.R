#' ============================================================
#' ML-BASED NORMALITY TESTING FRAMEWORK
#' Part 1: Setup, Helpers, Distribution Matching, Data Generation
#' ============================================================

## Optional: set this path manually before sourcing this file, if needed.
## Do not hard-code setwd() in a reusable script.
## Example:
project_dir <- "~/Desktop/OSU/Exploration"
if (dir.exists(project_dir)) setwd(project_dir)

#' ============================================================
#' SECTION 0: SETUP
#' ============================================================

RNGversion("4.2.0")

set.seed(
  12345,
  kind        = "Mersenne-Twister",
  normal.kind = "Inversion",
  sample.kind = "Rejection"
)

options(repos = c(CRAN = "https://cloud.r-project.org"))
options(mc_doScale_quiet = TRUE)

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}

pkgs <- c(
  ## Distributions
  "LaplacesDemon", "VGAM", "evd",
  
  ## Normality tests
  "nortest", "DescTools", "moments", "tseries",
  
  ## Feature engineering
  "Lmoments", "robustbase", "pracma", "ineq",
  
  ## Machine learning
  "caret", "glmnet", "randomForest", "gbm",
  "nnet", "kernlab", "e1071",
  
  ## Evaluation
  "pROC"
)

pacman::p_load(char = pkgs)


#' ============================================================
#' SECTION 1: SMALL HELPER FUNCTIONS
#' ============================================================

safe_calc <- function(expr, default = NA_real_) {
  out <- tryCatch(expr, error = function(e) default)
  
  if (is.null(out) || length(out) == 0L) {
    return(default)
  }
  
  out
}


make_output_dir <- function(output_dir) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  invisible(output_dir)
}


impute_na <- function(df, method = c("median", "mean")) {
  
  method <- match.arg(method)
  
  for (nm in names(df)) {
    
    if (!is.numeric(df[[nm]]) || !anyNA(df[[nm]])) {
      next
    }
    
    fill_value <- switch(
      method,
      median = median(df[[nm]], na.rm = TRUE),
      mean   = mean(df[[nm]], na.rm = TRUE)
    )
    
    if (!is.finite(fill_value)) {
      fill_value <- 0
    }
    
    df[[nm]][is.na(df[[nm]])] <- fill_value
  }
  
  df
}

#' ============================================================
#' SECTION 2: DEFAULT DISTRIBUTION SETS
#' ============================================================
#'
#' Each alternative distribution will be paired with its own
#' moment-matched Normal distribution.
#'
#' A user can add a new distribution in two ways:
#'   1. use a supported distribution name and parameters;
#'   2. provide normal_par = c(mean, sd) manually.
#'
#' This avoids keeping two separate lists:
#'   normal_configs
#'   nonnormal_specs
#'
#' which can easily become mismatched.
#' ============================================================


default_training_specs <- function() {
  list(
    list(dist = "gumbel",      par = c(0, 1)),
    list(dist = "chi_square",  par = 7),
    list(dist = "exponential", par = 3),
    list(dist = "weibull",     par = c(2, 1)),
    list(dist = "laplace",     par = c(0, 4)),
    list(dist = "beta",        par = c(8, 2)),
    list(dist = "uniform",     par = c(0, 10)),
    list(dist = "t",           par = 3),
    list(dist = "f",           par = c(6, 15)),
    list(dist = "cauchy",      par = c(0, 1))
  )
}


default_eval_specs <- function() {
  list(
    list(dist = "gumbel",       par = c(0, 1)),
    list(dist = "beta",         par = c(8, 2)),
    list(dist = "chi_square",   par = 3),
    list(dist = "laplace",      par = c(0, 4)),
    list(dist = "uniform",      par = c(0, 1)),
    list(dist = "cauchy",       par = c(0, 1)),
    list(dist = "contaminated", par = c(0.75, 0, 1, 5)),
    list(dist = "t",            par = 3)
  )
}


#' ============================================================
#' SECTION 3: MOMENT-MATCHED NORMAL PARAMETERS
#' ============================================================
#'
#' moment_matched_normal()
#'   Returns list(dist = "normal", par = c(mean, sd)).
#'
#' Priority:
#'   1. If the user supplies normal_par, use it.
#'   2. If the distribution is supported, use exact moments.
#'   3. If the user supplies rfun, estimate mean and sd by simulation.
#'   4. Otherwise, stop and ask for normal_par.
#'
#' For Cauchy:
#'   Mean and variance do not exist. We match the median and IQR.
#' ============================================================

moment_matched_normal <- function(spec,
                                  n_match = 1e6,
                                  seed    = 12345) {
  
  ## User-supplied matched Normal parameters
  if (!is.null(spec$normal_par)) {
    return(list(
      dist = "normal",
      par  = spec$normal_par
    ))
  }
  
  ## Custom random generator: estimate matching mean and SD
  if (!is.null(spec$rfun) && is.function(spec$rfun)) {
    set.seed(seed)
    x <- spec$rfun(n_match)
    
    return(list(
      dist = "normal",
      par  = c(mean(x), sd(x))
    ))
  }
  
  dist <- tolower(trimws(spec$dist))
  par  <- spec$par
  
  matched_par <- switch(
    dist,
    
    normal = {
      c(par[1], par[2])
    },
    
    gumbel = {
      loc   <- par[1]
      scale <- par[2]
      
      c(
        loc + scale * 0.5772156649,
        pi * scale / sqrt(6)
      )
    },
    
    chi_square = {
      df <- par[1]
      
      c(df, sqrt(2 * df))
    },
    
    gamma = {
      shape <- par[1]
      rate  <- par[2]
      
      c(shape / rate, sqrt(shape) / rate)
    },
    
    exponential = {
      rate <- par[1]
      
      c(1 / rate, 1 / rate)
    },
    
    weibull = {
      shape <- par[1]
      scale <- par[2]
      
      mu <- scale * gamma(1 + 1 / shape)
      sd <- scale * sqrt(
        gamma(1 + 2 / shape) - gamma(1 + 1 / shape)^2
      )
      
      c(mu, sd)
    },
    
    laplace = {
      location <- par[1]
      scale    <- par[2]
      
      c(location, sqrt(2) * scale)
    },
    
    beta = {
      a <- par[1]
      b <- par[2]
      
      mu <- a / (a + b)
      sd <- sqrt(a * b / ((a + b)^2 * (a + b + 1)))
      
      c(mu, sd)
    },
    
    uniform = {
      a <- par[1]
      b <- par[2]
      
      c((a + b) / 2, (b - a) / sqrt(12))
    },
    
    logistic = {
      location <- par[1]
      scale    <- par[2]
      
      c(location, scale * pi / sqrt(3))
    },
    
    lognormal = {
      meanlog <- par[1]
      sdlog   <- par[2]
      
      mu <- exp(meanlog + sdlog^2 / 2)
      sd <- sqrt((exp(sdlog^2) - 1) * exp(2 * meanlog + sdlog^2))
      
      c(mu, sd)
    },
    
    t = {
      df <- par[1]
      
      if (df <= 2) {
        stop("t distribution needs df > 2 for finite variance.")
      }
      
      c(0, sqrt(df / (df - 2)))
    },
    
    f = {
      df1 <- par[1]
      df2 <- par[2]
      
      if (df2 <= 4) {
        stop("F distribution needs df2 > 4 for finite variance.")
      }
      
      mu <- df2 / (df2 - 2)
      sd <- sqrt(
        2 * df2^2 * (df1 + df2 - 2) /
          (df1 * (df2 - 2)^2 * (df2 - 4))
      )
      
      c(mu, sd)
    },
    
    pareto = {
      alpha <- par[1]
      xm    <- par[2]
      
      if (alpha <= 2) {
        stop("Pareto distribution needs alpha > 2 for finite variance.")
      }
      
      mu <- xm * alpha / (alpha - 1)
      sd <- xm * sqrt(alpha / ((alpha - 1)^2 * (alpha - 2)))
      
      c(mu, sd)
    },
    
    cauchy = {
      location <- par[1]
      scale    <- par[2]
      
      ## Match Cauchy median and IQR to Normal scale.
      ## Cauchy IQR = 2 * scale.
      ## Normal IQR = 2 * qnorm(0.75) * sd.
      c(location, scale / qnorm(0.75))
    },
    
    contaminated = {
      p   <- par[1]
      mu  <- par[2]
      sd1 <- par[3]
      sd2 <- par[4]
      
      sd_mix <- sqrt(p * sd1^2 + (1 - p) * sd2^2)
      
      c(mu, sd_mix)
    },
    
    stop(
      "No automatic moment-matching rule for distribution: ", spec$dist, "\n",
      "Add normal_par = c(mean, sd) to the specification."
    )
  )
  
  list(
    dist = "normal",
    par  = matched_par
  )
}


make_paired_specs <- function(alt_specs) {
  
  lapply(alt_specs, function(alt) {
    list(
      normal = moment_matched_normal(alt),
      alt    = alt
    )
  })
}


#' ============================================================
#' SECTION 4: DATA GENERATION
#' ============================================================
#'
#' generate_data()
#'   Generates one sample from a named distribution.
#'
#' standardize = FALSE:
#'   Returns the raw sample.
#'
#' standardize = TRUE:
#'   Centers and scales using theoretical values where available.
#'   If theoretical SD is unavailable, IQR is used.
#'
#' Supported distributions:
#'   normal, gumbel, chi_square, gamma, exponential, weibull,
#'   laplace, beta, uniform, logistic, lognormal, t, f,
#'   pareto, cauchy, contaminated
#' ============================================================


get_distribution_moments <- function(dist, par) {
  
  dist <- tolower(trimws(dist))
  
  switch(
    dist,
    
    normal = {
      list(
        mean   = par[1],
        median = par[1],
        sd     = par[2],
        iqr    = qnorm(0.75, par[1], par[2]) - qnorm(0.25, par[1], par[2])
      )
    },
    
    gumbel = {
      loc   <- par[1]
      scale <- par[2]
      
      list(
        mean   = loc + scale * 0.5772156649,
        median = loc - scale * log(log(2)),
        sd     = pi * scale / sqrt(6),
        iqr    = qevd_gumbel(0.75, loc, scale) - qevd_gumbel(0.25, loc, scale)
      )
    },
    
    chi_square = {
      df <- par[1]
      
      list(
        mean   = df,
        median = qchisq(0.50, df),
        sd     = sqrt(2 * df),
        iqr    = qchisq(0.75, df) - qchisq(0.25, df)
      )
    },
    
    gamma = {
      shape <- par[1]
      rate  <- par[2]
      
      list(
        mean   = shape / rate,
        median = qgamma(0.50, shape = shape, rate = rate),
        sd     = sqrt(shape) / rate,
        iqr    = qgamma(0.75, shape = shape, rate = rate) - qgamma(0.25, shape = shape, rate = rate)
      )
    },
    
    exponential = {
      rate <- par[1]
      
      list(
        mean   = 1 / rate,
        median = log(2) / rate,
        sd     = 1 / rate,
        iqr    = qexp(0.75, rate = rate) - qexp(0.25, rate = rate)
      )
    },
    
    weibull = {
      shape <- par[1]
      scale <- par[2]
      
      list(
        mean   = scale * gamma(1 + 1 / shape),
        median = scale * log(2)^(1 / shape),
        sd     = scale * sqrt(
          gamma(1 + 2 / shape) - gamma(1 + 1 / shape)^2
        ),
        iqr    = qweibull(0.75, shape = shape, scale = scale) - qweibull(0.25, shape = shape, scale = scale)
      )
    },
    
    laplace = {
      location <- par[1]
      scale    <- par[2]
      
      list(
        mean   = location,
        median = location,
        sd     = sqrt(2) * scale,
        iqr    = 2 * scale * log(2)
      )
    },
    
    beta = {
      a <- par[1]
      b <- par[2]
      
      list(
        mean   = a / (a + b),
        median = qbeta(0.50, a, b),
        sd     = sqrt(a * b / ((a + b)^2 * (a + b + 1))),
        iqr    = qbeta(0.75, a, b) - qbeta(0.25, a, b)
      )
    },
    
    uniform = {
      a <- par[1]
      b <- par[2]
      
      list(
        mean   = (a + b) / 2,
        median = (a + b) / 2,
        sd     = (b - a) / sqrt(12),
        iqr    = (b - a) / 2
      )
    },
    
    logistic = {
      location <- par[1]
      scale    <- par[2]
      
      list(
        mean   = location,
        median = location,
        sd     = scale * pi / sqrt(3),
        iqr    = qlogis(0.75, location, scale) - qlogis(0.25, location, scale)
      )
    },
    
    lognormal = {
      meanlog <- par[1]
      sdlog   <- par[2]
      
      list(
        mean   = exp(meanlog + sdlog^2 / 2),
        median = exp(meanlog),
        sd     = sqrt((exp(sdlog^2) - 1) *
                        exp(2 * meanlog + sdlog^2)),
        iqr    = qlnorm(0.75, meanlog, sdlog) - qlnorm(0.25, meanlog, sdlog)
      )
    },
    
    t = {
      df <- par[1]
      
      list(
        mean   = if (df > 1) 0 else NA_real_,
        median = 0,
        sd     = if (df > 2) sqrt(df / (df - 2)) else NA_real_,
        iqr    = qt(0.75, df) - qt(0.25, df)
      )
    },
    
    f = {
      df1 <- par[1]
      df2 <- par[2]
      
      list(
        mean   = if (df2 > 2) df2 / (df2 - 2) else NA_real_,
        median = qf(0.50, df1, df2),
        sd     = if (df2 > 4) {
          sqrt( 2 * df2^2 * (df1 + df2 - 2) / (df1 * (df2 - 2)^2 * (df2 - 4)))
        } else {
          NA_real_
        },
        iqr    = qf(0.75, df1, df2) - qf(0.25, df1, df2)
      )
    },
    
    pareto = {
      alpha <- par[1]
      xm    <- par[2]
      
      list(
        mean   = if (alpha > 1) xm * alpha / (alpha - 1) else NA_real_,
        median = xm * 2^(1 / alpha),
        sd     = if (alpha > 2) {
          xm * sqrt(alpha / ((alpha - 1)^2 * (alpha - 2)))
        } else {
          NA_real_
        },
        iqr    = xm * (0.25)^(-1 / alpha) - xm * (0.75)^(-1 / alpha)
      )
    },
    
    cauchy = {
      location <- par[1]
      scale    <- par[2]
      
      list(
        mean   = NA_real_,
        median = location,
        sd     = NA_real_,
        iqr    = 2 * scale
      )
    },
    
    contaminated = {
      p   <- par[1]
      mu  <- par[2]
      sd1 <- par[3]
      sd2 <- par[4]
      
      sd_mix <- sqrt(p * sd1^2 + (1 - p) * sd2^2)
      
      list(
        mean   = mu,
        median = mu,
        sd     = sd_mix,
        iqr    = qnorm(0.75, mu, sd_mix) - qnorm(0.25, mu, sd_mix)
      )
    },
    
    stop("Unsupported distribution in get_distribution_moments(): ", dist)
  )
}


qevd_gumbel <- function(p, loc, scale) {
  loc - scale * log(-log(p))
}

#' Generate data
generate_data <- function(n,
                          dist,
                          par         = NULL,
                          standardize = FALSE,
                          center_by   = c("mean", "median", NULL)) {
  
  dist      <- tolower(trimws(dist))
  center_by <- match.arg(center_by)
  
  if (is.null(par)) {
    par <- switch(
      dist,
      normal       = c(0, 1),
      gumbel       = c(0, 1),
      chi_square   = 3,
      gamma        = c(3, 1),
      exponential  = 1,
      weibull      = c(2, 1),
      laplace      = c(0, 1),
      beta         = c(2, 5),
      uniform      = c(0, 1),
      logistic     = c(0, 1),
      lognormal    = c(0, 1),
      t            = 3,
      f            = c(6, 15),
      pareto       = c(3, 1),
      cauchy       = c(0, 1),
      contaminated = c(0.75, 0, 1, 5),
      stop("Unsupported distribution: ", dist)
    )
  }
  
  x <- switch(
    dist,
    
    normal = {
      rnorm(n, mean = par[1], sd = par[2])
    },
    
    gumbel = {
      evd::rgumbel(n, loc = par[1], scale = par[2])
    },
    
    chi_square = {
      rchisq(n, df = par[1])
    },
    
    gamma = {
      rgamma(n, shape = par[1], rate = par[2])
    },
    
    exponential = {
      rexp(n, rate = par[1])
    },
    
    weibull = {
      rweibull(n, shape = par[1], scale = par[2])
    },
    
    laplace = {
      LaplacesDemon::rlaplace(
        n,
        location = par[1],
        scale    = par[2]
      )
    },
    
    beta = {
      rbeta(n, shape1 = par[1], shape2 = par[2])
    },
    
    uniform = {
      runif(n, min = par[1], max = par[2])
    },
    
    logistic = {
      rlogis(n, location = par[1], scale = par[2])
    },
    
    lognormal = {
      rlnorm(n, meanlog = par[1], sdlog = par[2])
    },
    
    t = {
      rt(n, df = par[1])
    },
    
    f = {
      rf(n, df1 = par[1], df2 = par[2])
    },
    
    pareto = {
      VGAM::rpareto(n, shape = par[1], scale = par[2])
    },
    
    cauchy = {
      rcauchy(n, location = par[1], scale = par[2])
    },
    
    contaminated = {
      p   <- par[1]
      mu  <- par[2]
      sd1 <- par[3]
      sd2 <- par[4]
      
      group <- rbinom(n, size = 1, prob = p)
      sds   <- ifelse(group == 1, sd1, sd2)
      
      rnorm(n, mean = mu, sd = sds)
    },
    
    stop("Unsupported distribution: ", dist)
  )
  
  if (!standardize) {
    return(x)
  }
  
  moments <- get_distribution_moments(dist, par)
  
  center_value <- switch(
    center_by,
    mean   = moments$mean,
    median = moments$median
  )
  
  if (!is.finite(center_value)) {
    center_value <- moments$median
  }
  
  scale_value <- moments$sd
  
  if (!is.finite(scale_value) || scale_value <= 0) {
    scale_value <- moments$iqr
  }
  
  if (!is.finite(scale_value) || scale_value <= 0) {
    scale_value <- sd(x)
  }
  
  if (!is.finite(scale_value) || scale_value <= 0) {
    scale_value <- 1
  }
  
  x_std <- (x - center_value) / scale_value
  
  attr(x_std, "distribution") <- dist
  attr(x_std, "parameters")   <- par
  attr(x_std, "center_by")    <- center_by
  attr(x_std, "center_value") <- center_value
  attr(x_std, "scale_value")  <- scale_value
  attr(x_std, "standardized") <- TRUE
  
  x_std
}

#' ============================================================
#' SECTION 5: FEATURE ROW GENERATION
#' ============================================================
#'
#' This function draws repeated samples from one distribution,
#' computes features, adds the class label, and returns a data
#' frame ready for model training or evaluation.
#'
#' calculate_features() will be defined in Part 2.
#' ============================================================

generate_feature_rows <- function(n,
                                  n_rep,
                                  spec,
                                  label,
                                  feature_set = NULL,
                                  standardize_sample = FALSE,
                                  center_by = NULL) {
  
  rows <- vector("list", n_rep)
  
  for (i in seq_len(n_rep)) {
    
    x <- generate_data(
      n           = n,
      dist        = spec$dist,
      par         = spec$par,
      standardize = standardize_sample,
      center_by   = center_by
    )
    
    feats <- calculate_features(x)
    
    if (!is.null(feature_set)) {
      feats <- feats[, feature_set, drop = FALSE]
    }
    
    feats$Label <- label
    rows[[i]] <- feats
  }
  
  do.call(rbind, rows)
}


#' ============================================================
#' PART 2: FEATURE ENGINEERING
#' ============================================================
#'
#' This section computes numerical summaries of a sample.
#' Each row returned by calculate_features() represents one
#' simulated sample.
#'
#' The goal is not to decide normality from one statistic.
#' The goal is to give the ML models many complementary
#' signals about shape, tails, spacing, entropy, and Q-Q fit.
#' ============================================================

#' ============================================================
#' SECTION 6: BASIC FEATURE HELPERS
#' ============================================================

#' Standardize a numeric sample using its sample mean and sample SD.
#' Returns NA_real_ if the sample has no variance
standardize_sample <- function(x) {
  x <- as.numeric(na.omit(x))
  n <- length(x)
  
  if (n < 2L) {
    return(rep(NA_real_, n))
  }
  
  s <- sd(x)
  if (!is.finite(s) || s <= 0) {
    # Return explicit NAs if standardization is impossible
    return(rep(NA_real_, n))
  }
  
  as.numeric((x - mean(x)) / s)
}



#' Safely compute a ratio.
safe_ratio <- function(num, den) {
  
  if (!is.finite(num) || !is.finite(den)) {
    return(NA_real_)
  }
  
  if (abs(den) < .Machine$double.eps) {
    return(NA_real_)
  }
  
  num / den
}


#' Mean log spacing between consecutive order statistics.
mean_log_spacing <- function(x) {
  
  x <- sort(as.numeric(na.omit(x)))
  
  # Consecutive order-statistic spacings.
  d <- diff(x)
  
  # Log spacing requires strictly positive spacings.
  d <- d[d > 0]
  
  if (length(d) == 0L) {
    return(NA_real_)
  }
  
  mean(log(d))
}


#' ============================================================
#' SECTION 7: NORMALITY TEST STATISTICS
#' ============================================================

normality_test_stats <- function(x) {
 
  data.frame(
    SW_stat  = safe_calc(as.numeric(shapiro.test(x)$statistic)),
    SF_stat  = safe_calc(as.numeric(nortest::sf.test(x)$statistic)),
    AD_stat  = safe_calc(as.numeric(nortest::ad.test(x)$statistic)),
    LF_stat  = safe_calc(as.numeric(nortest::lillie.test(x)$statistic)),
    KS_stat  = safe_calc(as.numeric(ks.test(x, "pnorm", mean(x), sd(x))$statistic)),
    JB_stat  = safe_calc(as.numeric(tseries::jarque.bera.test(x)$statistic))
  )
}


#' ============================================================
#' SECTION 8: ENTROPY AND SPACING FEATURES
#' ============================================================
#'
#' 1. Lin-Mudholkar Zp statistic.
#'
#' The statistic is based on the correlation between the ordered
#' sample and leave-one-out variance terms. It is designed to
#' capture departures from Normality through the relationship
#' between order statistics and local dispersion.

calculate_zp_statistic <- function(x) {
  x <- sort(as.numeric(na.omit(x)))
  n <- length(x)

  if (n < 4L || length(unique(x)) < 4L) {
    return(NA_real_)
  }

  ## Leave-one-out variance terms.
  sum_x  <- sum(x)
  sum_x2 <- sum(x^2)
  n1     <- n - 1L

  u <- ((sum_x2 - x^2) - ((sum_x - x)^2 / n1))^(1 / 3)

  safe_calc(
    atanh(cor(x, u)),
    default = NA_real_
  )
}

#' ---------------------------------------------
#' 2. Vasicek entropy-based statistic.
#'
#' This statistic estimates entropy using order-statistic spacings.
#' For smooth unimodal distributions, the spacing structure is
#' informative about shape. Normal samples tend to produce a stable
#' spacing pattern after scaling by the sample SD.
calculate_vasicek_kmn <- function(x) {
  # Fast cleaning and sorting
  x <- sort(as.numeric(na.omit(x)))
  n <- length(x)
  
  # Fallback 1: Not enough data points to compute shape/entropy
  if (n < 4L) {
    return(NA_real_)
  }
  
  m <- floor(sqrt(n))
  s <- sd(x)
  
  # Fallback 2: Constant data has zero entropy variance
  if (!is.finite(s) || s <= 0.0) {
    return(NA_real_)
  }
  
  # Fast vectorized index bounding
  i <- seq_len(n)
  upper_idx <- pmin(n, i + m)
  lower_idx <- pmax(1L, i - m)
  
  spacings <- x[upper_idx] - x[lower_idx]
  
  # Fallback 3: Handle heavy ties or duplicates gracefully
  if (any(spacings <= 0.0)) {
    # Add a tiny epsilon to prevent log(0) while preserving feature signal
    spacings[spacings <= 0.0] <- 1e-7
  }
  
  # Compute entropy value
  entropy_val <- (n / (2 * m * s)) * exp(mean(log(spacings)))
  
  # Final safety check against Inf or NaN values
  if (!is.finite(entropy_val)) {
    return(NA_real_)
  }
  
  return(entropy_val)
}


#' --------------------------------------------------------
#' 3.  Rényi entropy using a histogram approximation.
#'
#' This summarizes how concentrated or dispersed the empirical
#' distribution is across bins. It is included as a broad shape
#' feature rather than as a formal Normality test.
#' alpha = 0 Hartley Entropy
#' alpha = 1 Shannon Entropy
#' alpha = 2 Collision Entropy
#' alpha = infinity Min- Entropy
renyi_entropy <- function(x, alpha = 2, bins = 10, global_range = NULL) {
  # Fast cleaning
  x <- as.numeric(na.omit(x))
  n <- length(x)
  
  # Fallback 1: Too few data points to build a histogram
  if (n < 3L) {
    return(NA_real_)
  }
  
  # Determine bin boundaries
  # Pass a fixed c(min, max) range vector to compare samples fairly
  if (!is.null(global_range)) {
    breaks <- seq(global_range[1], global_range[2], length.out = bins + 1)
  } else {
    breaks <- bins
  }
  
  # Compute histogram probabilities safely
  # explicit include.lowest handles edge cases
  p <- table(cut(x, breaks = breaks, include.lowest = TRUE))
  p <- as.numeric(p) / n
  
  # Filter out empty bins (0 * log(0) is 0 in entropy)
  p <- p[p > 0]
  if (length(p) == 0L) {
    return(NA_real_)
  }
  
  # Handle special case alpha = 1 (Shannon Entropy limit)
  if (abs(alpha - 1) < 1e-9) {
    entropy_val <- -sum(p * log(p))
  } else {
    entropy_val <- (1 / (1 - alpha)) * log(sum(p^alpha))
  }
  
  # Final safety check for numerical instability
  if (!is.finite(entropy_val)) {
    return(NA_real_)
  }
  
  return(entropy_val)
}


#' --------------------------------------------------------
#' 4. Approximate negentropy
#' Pre-computed reference constant for negentropy_approx().
negentropy_normal_ref <- local({
  set.seed(42L)
  # Pre-compute with a large draw .
  mean(log(cosh(rnorm(1e6L))))
})

#' Approximate negentropy.
#' Measures departure from Gaussianity 
negentropy_approx <- function(x) {
  z <- standardize_sample(x) 
  
  if (length(z) < 3L) {
    return(NA_real_)
  }
  
  # Log-Cosh Safe Approximation:
  # log(cosh(z)) is mathematically equal to |z| - log(2) + log(1 + exp(-2*|z|))
  abs_z <- abs(z)
  log_cosh_z <- abs_z - log(2) + log(1 + exp(-2 * abs_z))
  
  # Execute calculation using  safe_calc utility
  raw_negentropy <- safe_calc(
    (mean(log_cosh_z) - negentropy_normal_ref)^2,
    default = NA_real_
  )
  
  # Catch any boundary failures
  if (!is.finite(raw_negentropy)) {
    return(NA_real_)
  }
  
  return(raw_negentropy)
}

#' --------------------------------------------------------------
#' 5. Variance and skewness of log spacings.
#'
#' Consecutive order-statistic spacings describe how observations
#' are distributed across the support. The variance and skewness
#' of log spacings can detect clustering, gaps, and tail behavior.
calc_spacing_ratio_stats <- function(x) {
  x <- sort(as.numeric(na.omit(x)))
  
  # Compute consecutive differences
  d <- diff(x)
  d <- d[d > 0]
  
  # Fallback: Need at least 3 points to compute skewness of differences safely
  if (length(d) < 3L) {
    return(list(
      var_log_spacings = NA_real_,
      skew_log_spacings = NA_real_
    ))
  }
  
  log_d <- log(d)
  
  # Calculate variance safely
  v_val <- var(log_d)
  if (!is.finite(v_val)) {
    v_val <- NA_real_
  }
  
  # Calculate skewness safely using safe_calc utility
  s_val <- safe_calc(
    as.numeric(moments::skewness(log_d)),
    default = NA_real_
  )
  
  # Intercept NaN/Inf returns from moments package
  if (!is.finite(s_val)) {
    s_val <- NA_real_
  }
  
  list(
    var_log_spacings = v_val,
    skew_log_spacings = s_val
  )
}


#' --------------------------------------------------------------
#' 6. Greenwood Spacing Statistic
#'
#' Measures clustering and spacing unevenness.
#' Natively handles internal ties and applies precise sample size corrections.
greenwood_stat <- function(x) {
  x <- sort(as.numeric(na.omit(x)))
  
  # Remove duplicate ties to prevent artificial deflation of the statistic
  d <- diff(x)
  d <- d[d > 0]
  
  #  we need at least 3 valid intervals (4 points)
  # to compute a meaningful standardized variance structure
  n_spacings <- length(d)
  if (n_spacings < 3L) {
    return(NA_real_)
  }
  
  total_d <- sum(d)
  if (!is.finite(total_d) || total_d <= 0) {
    return(NA_real_)
  }
  
  # Raw Greenwood statistic (sum of squared normalized spacings)
  g <- sum(d^2) / (total_d^2)
  
  # EXACT uniform moments for m = n_spacings internal intervals sum to 1:
  # E[G] = 2 / (m + 1)
  # Var(G) = 4 * (m - 1) / ((m + 1)^2 * (m + 2) * (m + 3))
  m <- n_spacings
  expected_g <- 2 / (m + 1)
  var_g <- (4 * (m - 1)) / (((m + 1)^2) * (m + 2) * (m + 3))
  
  if (var_g <= 0) {
    return(NA_real_)
  }
  
  # Standardize using your safe_ratio helper
  standardized_g <- safe_ratio(g - expected_g, sqrt(var_g))
  
  # Catch any numerical boundary failures
  if (!is.finite(standardized_g)) {
    return(NA_real_)
  }
  
  return(standardized_g)
}

#' --------------------------------------------------------------
#' 7. Rao Spacing Statistic
#'
#' Measures the total absolute deviation of spacings from equal intervals.
rao_spacing_stat <- function(x) {
  x <- sort(as.numeric(na.omit(x)))
  n <- length(x)
  
  # Need at least 4 observations to compute an informative spacing profile
  if (n < 4L) {
    return(NA_real_)
  }
  
  # 1. Capture ALL n spacings, including the final wrap-around boundary gap
  internal_gaps <- diff(x)
  boundary_gap <- x[n] - x[1L] # Raw range used to scale the wrap
  
  # Filter out duplicate ties to prevent artificial inflation/division by zero
  internal_gaps <- internal_gaps[internal_gaps > 0]
  if (length(internal_gaps) == 0L || boundary_gap <= 0) {
    return(NA_real_)
  }
  
  # Combine internal intervals with the wrap-around interval
  d <- c(internal_gaps, boundary_gap)
  m <- length(d) # Total number of valid spacing intervals
  total_d <- sum(d)
  
  if (!is.finite(total_d) || total_d <= 0) {
    return(NA_real_)
  }
  
  # 2. Compute true normalized spacings (must sum to 1.0)
  s <- d / total_d
  
  # 3. Compute Rao's raw statistic: sum of absolute deviations from equal spacing
  raw_rao <- sum(abs(s - (1 / m)))
  
  # 4. Apply exact asymptotic limits for the uniform distribution:
  # Expected Mean -> 2 / e (~0.7357589)
  # Expected Variance -> (2*e - 5) / (e^2 * m)
  e_const <- exp(1)
  mean0 <- 2 / e_const
  var0 <- (2 * e_const - 5) / ((e_const^2) * m)
  
  if (var0 <= 0) {
    return(NA_real_)
  }
  
  # Standardize into a clean Z-score feature using safe_ratio helper
  standardized_rao <- safe_ratio(raw_rao - mean0, sqrt(var0))
  
  # Catch any unmanaged numerical instabilities or boundary failures
  if (!is.finite(standardized_rao)) {
    return(NA_real_)
  }
  
  return(standardized_rao)
}


#' ============================================================
#' SECTION 9: TAIL, OUTLIER, AND ROBUSTNESS FEATURES
#' ============================================================
#' -----------------------------------------------------------
#' 8. Tail Asymmetry Ratio Feature
#'
#' Compares the absolute structural span of the upper 10% tail 
#' to the lower 10% tail relative to the median.
calculate_tail_asymmetry <- function(x) {
  x <- sort(as.numeric(na.omit(x)))
  n <- length(x)
  
  # Need at least 10 observations to extract a valid 10% tail slice (1 element each)
  if (n < 10L) {
    return(NA_real_)
  }
  
  # Use the Median as the anchor instead of the Mean
  md <- median(x)
  
  # Calculate the exact number of elements that fit into a 10% tail slice
  k <- floor(n * 0.10)
  if (k < 1L) {
    return(NA_real_)
  }
  
  # Extract structurally identical slice lengths directly using sorted index positions
  lower_tail <- x[1:k]
  upper_tail <- x[(n - k + 1L):n]
  
  # Calculate absolute structural distance from the center anchor
  upper_excess <- mean(upper_tail) - md
  lower_deficit <- md - mean(lower_tail)
  
  # If the data is completely flat or has zero tail variance, return NA_real_
  if (!is.finite(upper_excess) || !is.finite(lower_deficit) || 
      upper_excess <= 0 || lower_deficit <= 0) {
    return(NA_real_)
  }
  
  # Compute the ratio via safe_ratio helper
  ratio_val <- safe_ratio(upper_excess, lower_deficit)
  
  # Catch any unmanaged numerical instabilities or boundary failures
  if (!is.finite(ratio_val)) {
    return(NA_real_)
  }
  
  return(ratio_val)
}

#' ---------------------------------------------------------------
#' 9. Moors Robust Kurtosis Feature
#'
#' Evaluates peak flatness and tail weight using octile points.
#' Completely independent of mean/SD standardization biases.
#' For a standard Normal distribution, this metric equals exactly 0.0.
calculate_moors_kurtosis <- function(x) {
  x <- as.numeric(na.omit(x))
  if (length(x) < 8L) {
    return(NA_real_)
  }
  
  # Extract the 7 precise octile boundaries (12.5% up to 87.5%)
  probs <- seq(0.125, 0.875, by = 0.125)
  q <- quantile(x, probs = probs, names = FALSE)
  
  # Moors Formula calculation
  num <- (q[7L] - q[5L]) + (q[3L] - q[1L])
  den <- q[6L] - q[2L]
  
  if (!is.finite(num) || !is.finite(den) || den <= 0) {
    return(NA_real_)
  }
  
  # Standardize so that a normal baseline shape maps perfectly to 0.0
  # Raw normal value for the Moors layout is ~1.2331
  moors_stat <- (num / den) - 1.233093
  
  if (!is.finite(moors_stat)) {
    # Handles extreme tie breakdowns safely
    return(NA_real_) 
  }
  
  return(moors_stat)
}


#' ----------------------------------------------------------------
#' 10. Excess Tail Proportion Feature
#'
#' Calculates the proportion of standardized observations beyond +/- 2 SD
#' relative to the exact mathematical standard Normal reference value.
calculate_excess_tail_prop <- function(x) {
  # Standardize the vector to mean 0, SD 1
  z <- standardize_sample(x)
  
  # Strip out NAs to get the true count of valid standardized elements
  z <- as.numeric(na.omit(z))
  n <- length(z)
  
  # Need an adequate sample size to calculate a meaningful proportion
  if (n < 5L) {
    return(NA_real_)
  }
  
  # Calculate the exact theoretical normal reference: P(|Z| > 2) = 0.04550026...
  normal_reference <- 2 * pnorm(-2)
  
  # Compute the empirical proportion of exceedances
  observed_proportion <- mean(abs(z) > 2.0)
  
  # Execute the ratio using the global safe_ratio helper
  ratio_val <- safe_ratio(observed_proportion, normal_reference)
  
  # Catch any numerical instabilities or boundary failures
  if (!is.finite(ratio_val)) {
    return(NA_real_)
  }
  
  return(ratio_val)
}


#' -------------------------------------------------------------------------
#' 11. Tukey Outlier Proportion Feature
#'
#' Calculates the proportion of observations falling outside the Tukey boxplot fences.
calculate_outlier_proportion <- function(x) {
  x <- as.numeric(na.omit(x))
  n <- length(x)
  
  # Need at least 4 elements to establish distinct quartile bounds
  if (n < 4L) {
    return(NA_real_)
  }
  
  # Extract quartiles exactly once 
  q <- quantile(x, probs = c(0.25, 0.75), names = FALSE)
  
  # Derive the robust scale interval directly
  interquartile_range <- q[2L] - q[1L]
  
  # If the data has heavy ties or is a flat line, IQR is 0. 
  if (!is.finite(interquartile_range) || interquartile_range <= 0) {
    return(NA_real_)
  }
  
  # Compute the standard 1.5 Tukey fence multiplier
  h <- 1.5 * interquartile_range
  lower_fence <- q[1L] - h
  upper_fence <- q[2L] + h
  
  # Calculate the empirical proportion of points breaching the fences
  outlier_prop <- mean(x < lower_fence | x > upper_fence)
  
  if (!is.finite(outlier_prop)) {
    return(NA_real_)
  }
  
  return(outlier_prop)
}


#' ----------------------------------------------------------------------------
#' 12. Qn Robust Scale Ratio Feature
#'
#' Compares a robust scale estimator (Qn) to an unbiased standard deviation.
calc_qn_robust_spread <- function(x) {
  x <- as.numeric(na.omit(x))
  n <- length(x)
  
  # Qn requires at least 2 points
  if (n < 5L) {
    return(NA_real_)
  }
  
  s <- sd(x)
  if (!is.finite(s) || s <= 0) {
    return(NA_real_)
  }
  
  # Compute the exact finite-sample correction factor c4(n) for the normal SD
  # This ensures E[s / c4(n)] = sigma, matching the calibration of Qn
  c4 <- sqrt(2 / (n - 1L)) * exp(lgamma(n / 2) - lgamma((n - 1L) / 2))
  unbiased_sd <- s / c4
  
  # Extract the robust scale estimator safely
  qn_val <- safe_calc(
    as.numeric(robustbase::Qn(x, finite.corr = TRUE)),
    default = NA_real_
  )
  
  if (is.na(qn_val) || !is.finite(qn_val) || qn_val <= 0) {
    return(NA_real_)
  }
  
  # Compute the scale ratio via the global safe_ratio helper
  ratio_val <- safe_ratio(qn_val, unbiased_sd)
  
  if (!is.finite(ratio_val)) {
    return(NA_real_)
  }
  
  return(ratio_val)
}


#' -----------------------------------------------------------------
#' 13. Medcouple Robust Skewness Feature
#'
#' Computes the medcouple estimator of skewness. Robust to outliers 
calc_medcouple <- function(x) {
  x <- as.numeric(na.omit(x))
  n <- length(x)
  
  # Medcouple requires at least 3 distinct observations
  if (n < 5L) {
    return(NA_real_)
  }
  
  # If the most common number takes up more than 50% of the row, 
  # the median collapses into a tie, causing robustbase::mc to freeze or crash.
  max_tie_count <- max(table(x))
  if (max_tie_count > (n / 2.0)) {
    return(NA_real_)
  }
  
  # Compute the medcouple metric safely
  mc_val <- safe_calc(
    as.numeric(robustbase::mc(x)),
    default = NA_real_
  )
  
  # Catch any unmanaged numerical instabilities or boundary failures
  if (!is.finite(mc_val)) {
    return(NA_real_)
  }
  
  return(mc_val)
}

#' -----------------------------------------------------------------------
#' 14. Robust Left and Right Tail Weights (Brys Methodology)
#'
#' Computes independent left and right tail heaviness metrics using octile ratios.
calc_tail_weights_brys <- function(x) {
  x <- as.numeric(na.omit(x))
  n <- length(x)
  
  # Brys octile tail metrics require at least 8 elements 
  if (n < 8L) {
    return(list(
      LTW = NA_real_,
      RTW = NA_real_
    ))
  }
  
  # Extract the specific quantiles needed for the Brys tail formulas:
  probs <- c(0.03125, 0.0625, 0.25, 0.50, 0.75, 0.9375, 0.96875)
  q <- quantile(x, probs = probs, names = FALSE)
  
  # Denominators: Inner body spread from the median on each side
  left_body  <- q[4L] - q[3L]  # Median - Q(0.25)
  right_body <- q[5L] - q[4L]  # Q(0.75) - Median
  
  # Guard against completely flat core bodies or heavy ties
  if (left_body <= 0 || right_body <= 0) {
    return(list(
      LTW = NA_real_,
      RTW = NA_real_
    ))
  }
  
  # TRUE BRYS FORMULAS:
  # Left Tail Weight (LTW): Ratio of extreme left spread to inner left spread
  ltw_val <- (q[2L] - q[1L]) / left_body
  
  # Right Tail Weight (RTW): Ratio of extreme right spread to inner right spread
  rtw_val <- (q[7L] - q[6L]) / right_body
  
  # Catch any floating-point or boundary anomalies
  if (!is.finite(ltw_val) || ltw_val < 0) ltw_val <- NA_real_
  if (!is.finite(rtw_val) || rtw_val < 0) rtw_val <- NA_real_
  
  list(
    LTW = ltw_val,
    RTW = rtw_val
  )
}



#' ============================================================
#' SECTION 10: Q-Q AND CORRELATION FEATURES
#' ============================================================
#' 
#' ---------------------------------------------------------------
#' 15. Orthogonal Cubic Coefficient from a Normal Q-Q Plot
#'
#' Fits a stable cubic orthogonal polynomial of empirical quantiles against 
#' theoretical Normal quantiles. 
qq_cubic_coef <- function(x) {
  # 1.  Standardize to eliminate scale dependencies
  z_sample <- standardize_sample(x)
  z_sample <- as.numeric(na.omit(z_sample))
  n <- length(z_sample)
  
  # Need at least 6 unique points to resolve a clean 3rd-degree polynomial fit
  if (n < 6L || length(unique(z_sample)) < 4L) {
    return(NA_real_)
  }
  
  # Theoretical Normal quantiles
  theo_q <- qnorm(ppoints(n))
  
  # Sort the standardized empirical observations
  emp_q <- sort(z_sample)
  
  # 2. Use orthogonal polynomials (raw = FALSE) to eliminate
  # severe multicollinearity between z and z^3, ensuring numeric stability.
  fit <- tryCatch({
    suppressWarnings(lm(emp_q ~ poly(theo_q, 3, raw = FALSE)))
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(fit)) {
    return(NA_real_)
  }
  
  coefs <- coef(fit)
  
  # Index 4 corresponds to the 3rd-degree (cubic) orthogonal polynomial term
  if (length(coefs) < 4L || !is.finite(coefs[4L])) {
    return(NA_real_)
  }
  
  return(as.numeric(coefs[4L]))
}


#' ----------------------------------------------------------------------
#' 16. de Wet-Venter Weighted Correlation Feature
#'
#' Computes the de Wet-Venter normality statistic, standardized into a 
#' sample-size invariant Z-score. 
dewet_venter_stat <- function(x) {
  x <- sort(as.numeric(na.omit(x)))
  n <- length(x)
  
  # The statistic requires enough order statistics to compute valid weights
  if (n < 8L) {
    return(NA_real_)
  }
  
  # Blom plotting positions
  p <- (seq_len(n) - 0.375) / (n + 0.25)
  m <- qnorm(p)
  
  # Density weights
  w <- dnorm(m)
  w_sum <- sum(w)
  
  if (!is.finite(w_sum) || w_sum <= 0) {
    return(NA_real_)
  }
  
  # Weighted centering
  x_bar <- sum(w * x) / w_sum
  m_bar <- sum(w * m) / w_sum
  x_dev <- x - x_bar
  m_dev <- m - m_bar
  
  # Calculate raw weighted correlation coefficient (r)
  numerator <- sum(w * x_dev * m_dev)
  denom_x <- sum(w * x_dev^2)
  denom_m <- sum(w * m_dev^2)
  
  if (denom_x <= 0 || denom_m <= 0) {
    return(NA_real_)
  }
  
  r <- safe_ratio(num = numerator, den = sqrt(denom_x * denom_m))
  if (!is.finite(r)) {
    return(NA_real_)
  }
  
  # Bound correlation strictly to mathematical limits
  r <- max(-1, min(1, r))
  r2 <- r^2
  
  # Convert raw R^2 into a standardized test statistic.
  # The raw statistic is L = n * (1 - R^2).
  raw_L <- n * (1 - r2)
  
  # De Wet & Venter's asymptotic normal moment corrections for sample size n:
  # E[L] ~ ln(n) + 0.5772 (Euler-Mascheroni constant) - 1
  # Var(L) ~ pi^2 / 6 - 1
  expected_L <- log(n) + 0.5772156649 - 1.0
  var_L <- (pi^2 / 6.0) - 1.0
  
  if (var_L <= 0) {
    return(NA_real_)
  }
  
  # Convert to a standardized scale-invariant Z-score
  standardized_z <- safe_ratio(raw_L - expected_L, sqrt(var_L))
  
  if (!is.finite(standardized_z)) {
    return(NA_real_)
  }
  
  return(standardized_z)
}


#' --------------------------------------------------------------
#' 17. Standardized Ryan-Joiner Candidate Feature
#'
#' Computes the Ryan-Joiner normality statistic, transformed into a 
#' sample-size invariant Z-score. Evaluates global Q-Q plot linearity evenly.
ryan_joiner_stat <- function(x) {
  x <- sort(as.numeric(na.omit(x)))
  n <- length(x)
  
  # Requires adequate order statistics to resolve the correlation structure
  if (n < 8L) {
    return(NA_real_)
  }
  
  # Blom plotting positions and expected Normal order statistics
  p <- (seq_len(n) - 0.375) / (n + 0.25)
  m <- qnorm(p)
  
  # Calculate raw correlation coefficient (r)
  r <- safe_calc(cor(x, m), default = NA_real_)
  
  if (is.na(r) || !is.finite(r)) {
    return(NA_real_)
  }
  
  # Bound numerically and compute R^2
  r <- max(-1, min(1, r))
  r2 <- r^2
  
  # Transform raw R^2 into a size-invariant Z-score.
  # Based on the Shapiro-Francia asymptotic normal properties for log(1 - R^2).
  raw_W <- 1 - r2
  
  # If raw_W is exactly 0 (perfect correlation), bound it to prevent log(0)
  if (raw_W <= 0) raw_W <- 1e-10
  
  log_W <- log(raw_W)
  
  # Asymptotic mean and SD formulas for log(1 - R^2) given sample size n
  expected_log_W <- -log(n) - 0.5772156649 + 1.4808
  sd_log_W       <- sqrt(1.5707963268 / n) # sqrt(pi / 2n)
  
  # Convert to a standardized Z-score
  standardized_z <- safe_ratio(log_W - expected_log_W, sd_log_W)
  
  if (!is.finite(standardized_z)) {
    return(NA_real_)
  }
  
  return(standardized_z)
}


#' ============================================================
#' SECTION 11: ECF AND DISTANCE FEATURES
#' ============================================================
#' 18. Epps-Pulley Geometric Normality Invariant
#'
#' Measures omnibus deviation from a normal distribution shape.
#' Returns a bounded metric: 0 means perfectly normal baseline shape.
epps_pulley_stat <- function(x) {
  x <- as.numeric(na.omit(x))
  n <- length(x)
  if (n < 4L) return(0.0)
  
  s <- sd(x)
  if (!is.finite(s) || s <= 0.0) return(0.0)
  
  # Center and scale data to standard normal space
  z <- (x - mean(x)) / s
  
  # Double-sum loop vectorized using outer-product matrices
  z_diff_sq <- outer(z, z, "-")^2
  z_sum_sq  <- outer(z, z, "+")^2
  
  # Epps-Pulley integration terms
  term_1 <- sum(exp(-z_diff_sq / 2.0)) / (n^2)
  term_2 <- (2.0 / n) * sum(exp(-z^2 / 4.0)) / sqrt(1.5)
  term_3 <- 1.0 / sqrt(2.0)
  
  statistic <- n * (term_1 - term_2 + term_3)
  
  if (!is.finite(statistic) || statistic < 0.0) return(0.0)
  return(statistic)
}

#' --------------------------------------------------------------
#' 19. Empirical Characteristic Function Discrete Deviations
#'
#' Extracts localized spatial shape deviations from the real and imaginary
#' coordinates of the empirical characteristic function.
#' Empirical Characteristic Function Discrete Deviations
ecf_deviations <- function(x, t_vals = c(0.5, 1.0, 2.0)) {
  x <- as.numeric(na.omit(x))
  
  # Generate explicit named slots for automated matrix column binding
  out_names <- c(
    paste0("ecf_re_t", seq_along(t_vals)),
    paste0("ecf_im_t", seq_along(t_vals))
  )
  
  # Initialize the target NA vector layout
  na_return <- setNames(rep(NA_real_, length(out_names)), out_names)
  
  if (length(x) < 3L) {
    return(na_return)
  }
  
  # Swapped '||' for '|' to evaluate vectors 
  if (any(t_vals > 3.0 | t_vals <= 0.0)) {
    return(na_return)
  }
  
  # Standardize sample to mean 0, variance 1
  z <- standardize_sample(x)
  z <- as.numeric(na.omit(z))
  
  # If the row has no variance, standardization returns NA
  if (length(z) < 3L) {
    return(na_return)
  }
  
  # Real deviations: Evaluates localized variance and kurtosis shapes
  real_dev <- vapply(
    t_vals, function(t) {
      val <- mean(cos(t * z)) - exp(-t^2 / 2)
      if (!is.finite(val)) NA_real_ else val
    }, numeric(1L)
  )
  
  # Imaginary deviations: robust asymmetry metric (0.0 for true normal)
  imag_dev <- vapply(
    t_vals, function(t) {
      val <- mean(sin(t * z))
      if (!is.finite(val)) NA_real_ else val
    }, numeric(1L)
  )
  
  final_vector <- c(real_dev, imag_dev)
  
  # Intercept any rare floating-point 
  if (any(!is.finite(final_vector))) {
    return(na_return)
  }
  
  setNames(final_vector, out_names)
}




#' ============================================================
#' SECTION 12: OTHER SHAPE FEATURES
#' ============================================================
#' 20. Excess Moors Kurtosis Feature
#'
#' Evaluates peak flatness and tail weight using octile intervals.
#' Calibrated so that a standard Normal distribution shape maps perfectly to 0.0.
calc_moors_kurtosis <- function(x) {
  x <- as.numeric(na.omit(x))
  if (length(x) < 8L) {
    return(NA_real_)
  }
  
  # Extract the 7 precise octile points
  octiles <- quantile(
    x, probs = seq(1 / 8, 7 / 8, by = 1 / 8), names = FALSE
  )
  
  numerator <- (octiles[7L] - octiles[5L]) + (octiles[3L] - octiles[1L])
  denominator <- octiles[6L] - octiles[2L]
  
  # Guard against zero-variance rows or dense identical ties
  if (!is.finite(numerator) || !is.finite(denominator) || denominator <= 0) {
    return(NA_real_)
  }
  
  raw_moors <- safe_ratio(numerator, denominator)
  
  if (is.na(raw_moors) || !is.finite(raw_moors)) {
    return(NA_real_)
  }
  
  # Subtract the theoretical normal constant
  # This sets the baseline expectation to exactly 0.0 under normality
  excess_moors <- raw_moors - 1.233093
  
  return(excess_moors)
}



#' ===================================================================
#' Master Feature Generation Engine
#'
#' Computes an expansive, highly rigorous pool of location, scale, spacing, 
#' curvature, and robust distribution-shape features for a single sample row.
#' Optimized for high-velocity gating and designed for rigorous feature selection.
#'
#' @param x Numeric vector representing a single distribution sample row.
#' @param feature_set Character vector of allowed feature names. If NULL, all are generated.
#' @param global_range Numeric vector c(min, max) passed to histogram-based entropy functions.
calculate_features <- function(x, feature_set = NULL, global_range = NULL) {
  x <- as.numeric(na.omit(x))
  n <- length(x)
  
  # Absolute baseline threshold to calculate minimum spacing and curvature properties
  if (n < 8L) {
    stop("calculate_features() requires at least 12 observations for stable shape extraction.")
  }
  
  # Fast, lightweight baseline scalars
  mean_x <- mean(x)
  median_x <- median(x)
  s <- sd(x)
  
  if (!is.finite(s) || s <= 0.0) {
    s <- NA_real_
  }
  
  # Dynamic feature-gating helper logic
  need <- function(...) {
    is.null(feature_set) || any(c(...) %in% feature_set)
  }
  
  # Initialize group gating states to optimize computational execution
  need_spacing      <- need("feat_vasicek", "feat_zp", "feat_rao_spacing", "feat_var_log_spacings", "feat_skew_log_spacings")
  need_renyi        <- need("feat_renyi_a0.5", "feat_renyi_a1.0", "feat_renyi_a2.0")
  need_tail_weights <- need("feat_ltw_brys", "feat_rtw_brys")
  need_ecf          <- need("feat_ecf_re_t1", "feat_ecf_re_t2", "feat_ecf_re_t3", "feat_ecf_im_t1", "feat_ecf_im_t2", "feat_ecf_im_t3")
  
  # Execute sub-group functions safely matching your audited signatures
  spacing_stats <- if (need_spacing) {
    calc_spacing_ratio_stats(x)
  } else {
    list(var_log_spacings = NA_real_, skew_log_spacings = NA_real_)
  }
  
  tail_weights <- if (need_tail_weights) {
    calc_tail_weights_brys(x)
  } else {
    list(LTW = NA_real_, RTW = NA_real_)
  }
  
  ecf_vals <- if (need_ecf) {
    ecf_deviations(x, t_vals = c(0.5, 1.0, 2.0))
  } else {
    setNames(rep(NA_real_, 6L), c("ecf_re_t1", "ecf_re_t2", "ecf_re_t3", "ecf_im_t1", "ecf_im_t2", "ecf_im_t3"))
  }
  
  # -------------------------------------------------------------------------
  # Core Feature Assembly (Vetted, Scale-Invariant Candidate Pool)
  # -------------------------------------------------------------------------
  feature_df <- data.frame(
    # Basic Location and Spread Base Features
    Mean                       = mean_x,
    Median                     = median_x,
    Mean_Median_Diff           = mean_x - median_x,
    Studentized_Range          = if (!is.na(s)) diff(range(x)) / s else NA_real_,
    
    # 1. Global Shape Distance & ECF Metrics
    feat_negentropy            = negentropy_approx(x),
    feat_epps_pulley           = epps_pulley_stat(x),
    feat_ecf_re_t1             = ecf_vals["ecf_re_t1"],
    feat_ecf_re_t2             = ecf_vals["ecf_re_t2"],
    feat_ecf_re_t3             = ecf_vals["ecf_re_t3"],
    feat_ecf_im_t1             = ecf_vals["ecf_im_t1"],
    feat_ecf_im_t2             = ecf_vals["ecf_im_t2"],
    feat_ecf_im_t3             = ecf_vals["ecf_im_t3"],
    
    # 2. Central Body Core Layout (Standardized Z-scores)
    feat_dewet_venter          = dewet_venter_stat(x),
    feat_ryan_joiner           = ryan_joiner_stat(x),
    
    # 3. Point Spacing & Order Intervals
    feat_zp                    = if (need_spacing) calculate_zp_statistic(x) else NA_real_,
    feat_rao_spacing           = if (need_spacing) rao_spacing_stat(x) else NA_real_,
    feat_vasicek               = if (need_spacing) calculate_vasicek_kmn(x) else NA_real_,
    feat_greenwood             = greenwood_stat(x),
    feat_var_log_spacings      = spacing_stats$var_log_spacings,
    feat_skew_log_spacings     = spacing_stats$skew_log_spacings,
    
    # 4. Geometric Q-Q Curvature Regression
    feat_qq_cubic_coef         = qq_cubic_coef(x),
    
    # 5. Robust Asymmetry & Skewness
    feat_tail_asymmetry        = calculate_tail_asymmetry(x),
    feat_medcouple             = calc_medcouple(x),
    
    # 6. Robust Tail Heaviness, Peak Flatness & Scale
    feat_moors_kurtosis        = calc_moors_kurtosis(x),
    feat_excess_tail_prop      = calculate_excess_tail_prop(x),
    feat_outlier_proportion    = calculate_outlier_proportion(x),
    feat_qn_robust_spread      = calc_qn_robust_spread(x),
    feat_ltw_brys              = tail_weights$LTW,
    feat_rtw_brys              = tail_weights$RTW,
    
    # 7. Consolidated Multi-Scale Entropy Spectrum Variants
    feat_renyi_a0.5            = if (need_renyi) renyi_entropy(x, alpha = 0.5, global_range = global_range) else NA_real_,
    feat_renyi_a1.0            = if (need_renyi) renyi_entropy(x, alpha = 1.0, global_range = global_range) else NA_real_,
    feat_renyi_a2.0            = if (need_renyi) renyi_entropy(x, alpha = 2.0, global_range = global_range) else NA_real_,
    
    stringsAsFactors           = FALSE
  )
  
  # -------------------------------------------------------------------------
  # Mathematically Defensible Interactions & Transforms
  # -------------------------------------------------------------------------
  # We construct interactions ONLY using vetted, scale-invariant robust metrics
  feature_df$inter_asym_vs_moors <- safe_calc(feature_df$feat_tail_asymmetry * feature_df$feat_moors_kurtosis)
  feature_df$inter_negent_vs_epps <- safe_calc(feature_df$feat_negentropy * feature_df$feat_epps_pulley)
  
  if (!is.null(feature_set)) {
    missing_features <- setdiff(feature_set, names(feature_df))

    if (length(missing_features) > 0L) {
      stop(
        "The following requested features were not generated: ",
        paste(missing_features, collapse = ", ")
      )
    }

    feature_df <- feature_df[, feature_set, drop = FALSE]
  }

  # Return the pure row dataframe with un-mutated NA_real_ components intact
  return(feature_df)
}


#' ============================================================
#' SECTION 14: FEATURE PREPROCESSING
#' ============================================================
#' Streamlined Preprocessing Engine
#'
#' Applies missing-value imputation followed by centering and scaling to 
#' the numeric feature columns, while strictly protecting target labels from leakage.
preprocess_data <- function(train_data,
                            impute_method = c("median", "mean"),
                            scale_method  = c("center", "scale")) {

  impute_method <- match.arg(impute_method)

  if (!"Label" %in% names(train_data)) {
    stop("preprocess_data() requires a Label column.")
  }

  ## Keep the outcome separate to prevent leakage.
  y_train <- train_data$Label

  feature_df <- train_data[, setdiff(names(train_data), "Label"), drop = FALSE]
  numeric_cols <- vapply(feature_df, is.numeric, logical(1L))
  numeric_train <- feature_df[, numeric_cols, drop = FALSE]

  if (ncol(numeric_train) == 0L) {
    stop("No numeric feature columns found for preprocessing.")
  }

  ## Convert Inf and NaN to NA before imputation.
  numeric_train[] <- lapply(numeric_train, function(z) {
    z <- as.numeric(z)
    z[!is.finite(z)] <- NA_real_
    z
  })

  ## Drop columns that are entirely missing or constant.
  keep_cols <- vapply(numeric_train, function(z) {
    any(is.finite(z)) && length(unique(z[is.finite(z)])) > 1L
  }, logical(1L))

  numeric_train <- numeric_train[, keep_cols, drop = FALSE]

  if (ncol(numeric_train) == 0L) {
    stop("All feature columns were empty or constant after cleaning.")
  }

  ## Simple training-only imputation.
  impute_values <- vapply(numeric_train, function(z) {
    val <- switch(
      impute_method,
      median = median(z, na.rm = TRUE),
      mean   = mean(z, na.rm = TRUE)
    )

    if (!is.finite(val)) 0 else val
  }, numeric(1L))

  for (nm in names(numeric_train)) {
    numeric_train[[nm]][is.na(numeric_train[[nm]])] <- impute_values[[nm]]
  }

  ## Centering/scaling is fit on the training features only.
  preproc_obj <- caret::preProcess(
    numeric_train,
    method = scale_method
  )

  train_transformed <- predict(preproc_obj, numeric_train)
  train_transformed$Label <- y_train

  list(
    train         = train_transformed,
    preProc      = preproc_obj,
    impute_values = impute_values,
    kept_features = names(numeric_train),
    feature_names = names(train_transformed)[names(train_transformed) != "Label"]
  )
}


#' Apply saved preprocessing to new feature rows.
apply_preprocess <- function(new_data, prep_obj) {

  if (!all(c("preProc", "impute_values", "kept_features") %in% names(prep_obj))) {
    stop("prep_obj must be the object returned by preprocess_data().")
  }

  x <- new_data[, prep_obj$kept_features, drop = FALSE]

  x[] <- lapply(x, function(z) {
    z <- as.numeric(z)
    z[!is.finite(z)] <- NA_real_
    z
  })

  for (nm in names(x)) {
    x[[nm]][is.na(x[[nm]])] <- prep_obj$impute_values[[nm]]
  }

  predict(prep_obj$preProc, x)
}


#' ============================================================
#' SECTION 15: TRAINING DATA GENERATION
#' ============================================================
#' Balanced Training Data Generation Engine
#'
#' Generates simulated feature matrices from paired moment-matched Normal and 
#' Alternative distribution specifications to force shape-based classification learning.

make_training_data <- function(n,
                               paired_specs,
                               num_sim,
                               feature_set = NULL,
                               standardize_sample = FALSE,
                               center_by = NULL) {
  
  rows <- vector("list", 2L * length(paired_specs))
  row_id <- 0L
  
  for (pair in paired_specs) {
    
    row_id <- row_id + 1L
    rows[[row_id]] <- generate_feature_rows(
      n                  = n,
      n_rep              = num_sim,
      spec               = pair$normal,
      label              = "Normal",
      feature_set         = feature_set,
      standardize_sample  = standardize_sample,
      center_by           = center_by
    )
    
    row_id <- row_id + 1L
    rows[[row_id]] <- generate_feature_rows(
      n                  = n,
      n_rep              = num_sim,
      spec               = pair$alt,
      label              = "Non_Normal",
      feature_set         = feature_set,
      standardize_sample  = standardize_sample,
      center_by           = center_by
    )
  }
  
  train_data <- do.call(rbind, rows)
  
  train_data$Label <- factor(
    train_data$Label,
    levels = c("Non_Normal", "Normal")
  )
  
  rownames(train_data) <- NULL
  
  train_data
}

#' ============================================================
#' SECTION 16: MODEL TRAINING HELPERS
#' ============================================================
#' Caret Cross-Validation Settings
#'
#' Establishes balanced 10-fold cross-validation 
make_cv_control <- function(k_folds = 10) {
  caret::trainControl(
    method          = "cv",
    number          = k_folds,
    classProbs      = TRUE,
    summaryFunction = caret::twoClassSummary,
    savePredictions = "final"
  )
}

#' Train One Robust Classifier
train_one_model <- function(model_name, train_data, cv_ctrl) {
  
  # Determine actual feature dimensions dynamically (subtracting the 'Label' column)
  num_features <- ncol(train_data) - 1L
  
  # -------------------------------------------------------------------------
  # Logistic Regression via Elastic Net (Glmnet)
  # -------------------------------------------------------------------------
  if (model_name == "LR") {
    return(
      suppressWarnings(
        caret::train(
           Label ~ ., 
          data      = train_data, 
          method    = "glmnet", 
          metric    = "ROC", 
          trControl = cv_ctrl,
          tuneGrid  = expand.grid(
            alpha  = c(0.0, 0.5, 1.0),
            lambda = 10^seq(-4, 1, length.out = 20L)
          ),
          family    = "binomial"
        )
      )
    )
  }
  
  # -------------------------------------------------------------------------
  # Random Forest (Dynamic Tuning Bounds)
  # -------------------------------------------------------------------------
  if (model_name == "RF") {
    # Clip mtry choices so they never exceed total feature counts
    mtry_pool <- c(2, 5, 10, 15)
    valid_mtry <- unique(pmin(num_features, mtry_pool[mtry_pool <= num_features]))
    if (length(valid_mtry) == 0L) valid_mtry <- 1L
    
    return(
      caret::train(
         Label ~ ., 
        data      = train_data, 
        method    = "rf", 
        metric    = "ROC", 
        trControl = cv_ctrl,
        tuneGrid  = expand.grid(mtry = valid_mtry),
        ntree     = 500L,
        importance = TRUE
      )
    )
  }
  
  # -------------------------------------------------------------------------
  # Artificial Neural Network (Nnet)
  # -------------------------------------------------------------------------
  if (model_name == "ANN") {
    return(
      caret::train(
         Label ~ ., 
        data      = train_data, 
        method    = "nnet", 
        metric    = "ROC", 
        trControl = cv_ctrl,
        tuneGrid  = expand.grid(
          size  = c(5, 10),
          decay = c(0.001, 0.01, 0.10)
        ),
        MaxNWts   = 5000L,
        maxit     = 300L,
        trace     = FALSE
      )
    )
  }
  
  # -------------------------------------------------------------------------
  # Gradient Boosting Machine (Gbm)
  # -------------------------------------------------------------------------
  if (model_name == "GBM") {
    # Ensure minimum node size scales safely with sample availability
    min_node_size <- if (nrow(train_data) < 50L) c(2, 5) else c(5, 10)
    
    return(
      caret::train(
         Label ~ ., 
        data      = train_data, 
        method    = "gbm", 
        metric    = "ROC", 
        trControl = cv_ctrl,
        tuneGrid  = expand.grid(
          n.trees           = c(100, 300),
          interaction.depth = c(1, 3),
          shrinkage         = c(0.01, 0.10),
          n.minobsinnode    = min_node_size
        ),
        verbose   = FALSE
      )
    )
  }
  
  # -------------------------------------------------------------------------
  # Custom e1071 LIBSVM radial container (Platt Probability Calibrated)
  # -------------------------------------------------------------------------
  if (model_name == "SVM") {
    svm_e1071 <- list(
      label      = "Support Vector Machine (RBF, e1071)",
      type       = "Classification",
      library    = "e1071",
      loop       = NULL,
      parameters = data.frame(
        parameter = c("cost", "gamma"),
        class     = c("numeric", "numeric"),
        label     = c("Cost", "Gamma")
      ),
      grid       = function(x, y, len = NULL, search = "grid") {
        expand.grid(cost = c(0.1, 1, 10), gamma = c(0.01, 0.1, 1))
      },
      fit        = function(x, y, wts, param, lev, last, weights, classProbs, ...) {
        fit <- e1071::svm(
          x           = x,
          y           = y,
          kernel      = "radial",
          cost        = param$cost,
          gamma       = param$gamma,
          probability = TRUE,
          ...
        )

        fit$obsLevels <- lev
        fit
      },
      predict    = function(modelFit, newdata, preProc = NULL, submodels = NULL) {
        predict(modelFit, newdata = newdata)
      },
      prob       = function(modelFit, newdata, preProc = NULL, submodels = NULL) {
        prob_mat <- attr(
          predict(modelFit, newdata = newdata, probability = TRUE),
          "probabilities"
        )

        if (is.null(prob_mat)) {
          stop("SVM probability estimates were not returned.")
        }

        prob_mat[, modelFit$obsLevels, drop = FALSE]
      },
      sort       = function(x) x[order(x$cost, x$gamma), ]
    )
    
    return(
      caret::train(
         Label ~ ., 
        data      = train_data, 
        method    = svm_e1071, 
        metric    = "ROC", 
        trControl = cv_ctrl,
        tuneGrid  = expand.grid(
          cost  = c(0.1, 0.5, 1.0, 5.0, 10.0),
          gamma = c(0.01, 0.05, 0.10, 0.50, 1.00)
        )
      )
    )
  }
  
  # -------------------------------------------------------------------------
  # K-Nearest Neighbors (Knn)
  # -------------------------------------------------------------------------
  if (model_name == "KNN") {
    # Clip neighbor maximum bounds based on total row count limits
    max_k <- max(3L, min(11L, floor(nrow(train_data) / 3.0)))
    k_sequence <- seq(3L, max_k, by = 2L)
    
    return(
      caret::train(
         Label ~ ., 
        data      = train_data, 
        method    = "knn", 
        metric    = "ROC", 
        trControl = cv_ctrl,
        tuneGrid  = expand.grid(k = k_sequence)
      )
    )
  }
  
  stop("Unsupported model name: ", model_name)
}


#' ============================================================
#' SECTION 17: TRAIN MODELS FOR ONE SAMPLE SIZE
#' ============================================================
#' ============================================================
#' SECTION 17: TRAIN MODELS FOR ONE SAMPLE SIZE
#' ============================================================

train_models_one_n <- function(n,
                               paired_specs,
                               num_sim = 100,
                               models_to_train = c("LR", "RF", "ANN", "GBM", "SVM", "KNN"),
                               feature_set = NULL,
                               standardize_sample = FALSE,
                               center_by = NULL,
                               k_folds = 10) {
  
  valid_models <- c("LR", "RF", "ANN", "GBM", "SVM", "KNN")
  bad_models   <- setdiff(models_to_train, valid_models)
  
  if (length(bad_models) > 0L) {
    stop("Unsupported model(s): ", paste(bad_models, collapse = ", "))
  }
  
  center_by <- match.arg(center_by, choices = c("mean", "median", NULL))
  
  cat("\n")
  cat("Training models for n =", n, "\n")
  cat("Distributions:", length(paired_specs), "matched pairs\n")
  cat("Replicates per class per pair:", num_sim, "\n")
  cat("Standardize samples:", standardize_sample, "\n")
  cat("Center by:", center_by, "\n")
  
  raw_train <- make_training_data(
    n                  = n,
    paired_specs       = paired_specs,
    num_sim            = num_sim,
    feature_set        = feature_set,
    standardize_sample = standardize_sample,
    center_by          = center_by
  )
  
  prep <- preprocess_data(
    train_data    = raw_train,
    impute_method = "median",
    scale_method  = c("center", "scale")
  )
  
  train_std <- prep$train
  cv_ctrl   <- make_cv_control(k_folds = k_folds)
  
  models <- vector("list", length(models_to_train))
  names(models) <- models_to_train
  
  for (model_name in models_to_train) {
    
    cat("  Fitting", model_name, "...\n")
    
    models[[model_name]] <- train_one_model(
      model_name = model_name,
      train_data = train_std,
      cv_ctrl    = cv_ctrl
    )
  }
  
  list(
    n                  = n,
    models             = models,
    prep_obj           = prep$preProcStandard,
    feature_names      = prep$feature_names,
    paired_specs       = paired_specs,
    standardize_sample = standardize_sample,
    center_by          = center_by
  )
}

#' ============================================================
#' SECTION 18: TRAIN MODELS ACROSS SAMPLE SIZES
#' ============================================================
#'
#' Trains the requested models for each sample size.
#'
#' Returns:
#'   A named list, where each element is the trained model bundle
#'   for one sample size.
#'
#' Example:
#'   trained[["30"]]$models$RF
#'   trained[["50"]]$models$SVM
#' ============================================================

train_models_all_n <- function(sample_sizes,
                               paired_specs,
                               num_sim = 100,
                               models_to_train = c("LR", "RF", "ANN", "GBM", "SVM", "KNN"),
                               feature_set = NULL,
                               standardize_sample = FALSE,
                               center_by = NULL,
                               k_folds = 10) {
  
  if (length(sample_sizes) == 0L) {
    stop("sample_sizes must contain at least one sample size.")
  }
  
  trained <- vector("list", length(sample_sizes))
  names(trained) <- as.character(sample_sizes)
  
  for (n in sample_sizes) {
    
    cat("\n============================================================\n")
    cat("Training models for sample size n =", n, "\n")
    cat("============================================================\n")
    
    trained[[as.character(n)]] <- train_models_one_n(
      n                  = n,
      paired_specs       = paired_specs,
      num_sim            = num_sim,
      models_to_train    = models_to_train,
      feature_set        = feature_set,
      standardize_sample = standardize_sample,
      center_by          = center_by,
      k_folds            = k_folds
    )
  }
  
  trained
}


#' ============================================================
#' PART 4: PREDICTION AND EVALUATION
#' ============================================================


#' ============================================================
#' SECTION 19: PREDICTION HELPERS
#' ============================================================

#' Resolve centering choice for optional sample standardization.
#'
#' If standardize_sample = FALSE, no centering is applied to the raw sample.
#' If standardize_sample = TRUE and center_by is NULL, mean-centering is used.
resolve_center_by <- function(standardize_sample = FALSE, center_by = NULL) {
  
  if (!isTRUE(standardize_sample)) {
    return("mean")  # ignored by generate_data() when standardize = FALSE
  }
  
  if (is.null(center_by)) {
    return("mean")
  }
  
  match.arg(center_by, choices = c("mean", "median"))
}


#' Prepare one sample for prediction.
#'
#' Computes only the features used during training and applies the
#' preprocessing object fitted on the training data.
prepare_prediction_features <- function(x, model_bundle) {
  
  if (is.null(model_bundle$feature_names)) {
    stop("model_bundle$feature_names is missing.")
  }
  
  feats <- calculate_features(
    x           = x,
    feature_set = model_bundle$feature_names
  )
  
  missing_features <- setdiff(model_bundle$feature_names, names(feats))
  
  if (length(missing_features) > 0L) {
    stop(
      "Missing feature(s) from calculate_features(): ",
      paste(missing_features, collapse = ", ")
    )
  }
  
  feats <- feats[, model_bundle$feature_names, drop = FALSE]
  
  if (!is.null(model_bundle$prep_obj)) {
    feats <- predict(model_bundle$prep_obj, feats)
  }
  
  feats
}


#' Prepare a batch of feature rows for prediction.
prepare_prediction_matrix <- function(feature_df, model_bundle) {
  
  if (is.null(model_bundle$feature_names)) {
    stop("model_bundle$feature_names is missing.")
  }
  
  missing_features <- setdiff(model_bundle$feature_names, names(feature_df))
  
  if (length(missing_features) > 0L) {
    stop(
      "Evaluation feature matrix is missing feature(s): ",
      paste(missing_features, collapse = ", ")
    )
  }
  
  feature_df <- feature_df[, model_bundle$feature_names, drop = FALSE]
  
  feature_df[] <- lapply(feature_df, function(z) {
    z <- as.numeric(z)
    z[!is.finite(z)] <- NA_real_
    z
  })
  
  if (!is.null(model_bundle$prep_obj)) {
    feature_df <- predict(model_bundle$prep_obj, feature_df)
  }
  
  feature_df
}


#' Predict P(Non_Normal) from all trained models.
#'
#' Returns a data frame with one probability column per model.
predict_model_probs <- function(model_bundle, feats_std) {
  
  model_names <- names(model_bundle$models)
  
  prob_df <- lapply(model_names, function(model_name) {
    
    pred <- predict(
      model_bundle$models[[model_name]],
      newdata = feats_std,
      type    = "prob"
    )
    
    if (!"Non_Normal" %in% names(pred)) {
      stop("Model ", model_name, " did not return a 'Non_Normal' probability column.")
    }
    
    as.numeric(pred[, "Non_Normal"])
  })
  
  prob_df <- as.data.frame(prob_df)
  names(prob_df) <- model_names
  
  prob_df
}


#' Fast single-model, single-sample prediction.
#'
#' Returns P(Normal) for one sample and one trained model.
predict_prob_normal_fast <- function(x, model_bundle, model_name = "RF") {
  
  if (!model_name %in% names(model_bundle$models)) {
    stop("Model not found in model_bundle: ", model_name)
  }
  
  feats_std <- prepare_prediction_features(
    x            = x,
    model_bundle = model_bundle
  )
  
  pred_probs <- predict(
    model_bundle$models[[model_name]],
    newdata = feats_std,
    type    = "prob"
  )
  
  if (!"Normal" %in% names(pred_probs)) {
    stop("Model ", model_name, " did not return a 'Normal' probability column.")
  }
  
  as.numeric(pred_probs[, "Normal"])
}


#' Convert P(Non_Normal) to a class label.
class_from_prob <- function(prob, threshold = 0.50) {
  
  ifelse(prob >= threshold, "Non_Normal", "Normal")
}


#' Create a readable distribution label.
make_distribution_label <- function(spec) {
  
  paste0(
    spec$dist,
    "(",
    paste(round(spec$par, 4), collapse = ", "),
    ")"
  )
}


#' ============================================================
#' SECTION 20: EVALUATION DATA GENERATION
#' ============================================================

#' Create flattened evaluation specifications.
#'
#' For each alternative distribution, this returns:
#'   1. the moment-matched Normal distribution;
#'   2. the corresponding alternative distribution.
make_eval_specs <- function(alt_specs) {
  
  pairs <- make_paired_specs(alt_specs)
  
  eval_specs <- vector("list", 2L * length(pairs))
  row_id <- 0L
  
  for (pair in pairs) {
    
    row_id <- row_id + 1L
    eval_specs[[row_id]] <- list(
      spec  = pair$normal,
      label = "Normal"
    )
    
    row_id <- row_id + 1L
    eval_specs[[row_id]] <- list(
      spec  = pair$alt,
      label = "Non_Normal"
    )
  }
  
  eval_specs
}


#' Generate evaluation features for one distribution specification.
#'
#' This returns a feature matrix plus the corresponding true labels
#' and distribution labels.
generate_eval_block <- function(n,
                                n_iter,
                                eval_item,
                                feature_names,
                                standardize_sample = FALSE,
                                center_by = NULL) {
  
  spec       <- eval_item$spec
  true_class <- eval_item$label
  dist_label <- make_distribution_label(spec)
  
  center_by_use <- resolve_center_by(
    standardize_sample = standardize_sample,
    center_by          = center_by
  )
  
  rows <- vector("list", n_iter)
  
  for (i in seq_len(n_iter)) {
    
    x <- generate_data(
      n           = n,
      dist        = spec$dist,
      par         = spec$par,
      standardize = standardize_sample,
      center_by   = center_by_use
    )
    
    rows[[i]] <- calculate_features(
      x           = x,
      feature_set = feature_names
    )
  }
  
  feature_df <- do.call(rbind, rows)
  rownames(feature_df) <- NULL
  
  list(
    features     = feature_df,
    true_class   = rep(true_class, n_iter),
    distribution = rep(dist_label, n_iter)
  )
}


#' ============================================================
#' SECTION 21: MODEL EVALUATION FOR ONE SAMPLE SIZE
#' ============================================================

#' Evaluate trained models for one sample size.
#'
#' For each evaluation distribution:
#'   1. generate n_iter samples;
#'   2. compute the trained feature set;
#'   3. apply the training preprocessing object;
#'   4. predict P(Non_Normal);
#'   5. return prediction data frames for each model and MajorityVote.
evaluate_one_n <- function(n,
                           model_bundle,
                           eval_specs,
                           n_iter = 1000,
                           threshold = 0.50,
                           majority_threshold = 0.50,
                           standardize_sample = FALSE,
                           center_by = NULL,
                           store_holdout_features = FALSE) {
  
  if (is.null(model_bundle$models) || length(model_bundle$models) == 0L) {
    stop("model_bundle$models is missing or empty.")
  }
  
  if (is.null(model_bundle$feature_names)) {
    stop("model_bundle$feature_names is missing.")
  }
  
  model_names <- names(model_bundle$models)
  all_names   <- c(model_names, "MajorityVote")
  
  cat("\n")
  cat("Evaluating n =", n, "\n")
  cat("Evaluation specs:", length(eval_specs), "\n")
  cat("Replicates per spec:", n_iter, "\n")
  cat("Standardize samples:", standardize_sample, "\n")
  
  if (isTRUE(standardize_sample)) {
    cat("Center by:", resolve_center_by(TRUE, center_by), "\n")
  }
  
  block_results <- vector("list", length(eval_specs))
  
  for (j in seq_along(eval_specs)) {
    
    eval_item <- eval_specs[[j]]
    
    cat("  Evaluating", make_distribution_label(eval_item$spec), "\n")
    
    eval_block <- generate_eval_block(
      n                  = n,
      n_iter             = n_iter,
      eval_item          = eval_item,
      feature_names      = model_bundle$feature_names,
      standardize_sample = standardize_sample,
      center_by          = center_by
    )
    
    feats_std <- prepare_prediction_matrix(
      feature_df   = eval_block$features,
      model_bundle = model_bundle
    )
    
    probs <- predict_model_probs(
      model_bundle = model_bundle,
      feats_std    = feats_std
    )
    
    block_results[[j]] <- list(
      probs        = probs,
      true_class   = eval_block$true_class,
      distribution = eval_block$distribution,
      features     = if (store_holdout_features) feats_std else NULL
    )
  }
  
  true_class_all <- unlist(
    lapply(block_results, `[[`, "true_class"),
    use.names = FALSE
  )
  
  distribution_all <- unlist(
    lapply(block_results, `[[`, "distribution"),
    use.names = FALSE
  )
  
  prob_all <- do.call(
    rbind,
    lapply(block_results, `[[`, "probs")
  )
  
  rownames(prob_all) <- NULL
  
  pred_store <- vector("list", length(all_names))
  names(pred_store) <- all_names
  
  for (model_name in model_names) {
    
    prob_nn <- prob_all[[model_name]]
    
    pred_store[[model_name]] <- data.frame(
      True_Class      = true_class_all,
      Predicted_Class = class_from_prob(prob_nn, threshold),
      Prob_Non_Normal = prob_nn,
      Distribution    = distribution_all,
      stringsAsFactors = FALSE
    )
  }
  
  mean_prob <- rowMeans(prob_all[, model_names, drop = FALSE], na.rm = TRUE)
  
  pred_store[["MajorityVote"]] <- data.frame(
    True_Class      = true_class_all,
    Predicted_Class = class_from_prob(mean_prob, majority_threshold),
    Prob_Non_Normal = mean_prob,
    Distribution    = distribution_all,
    stringsAsFactors = FALSE
  )
  
  holdout_features <- NULL
  
  if (isTRUE(store_holdout_features)) {
    holdout_features <- do.call(
      rbind,
      lapply(block_results, `[[`, "features")
    )
    
    rownames(holdout_features) <- NULL
  }
  
  list(
    predictions      = pred_store,
    holdout_features = holdout_features
  )
}

#' The ML classifier returns an estimated probability that a sample belongs 
#' to the non-normal class. For the primary analysis, a fixed threshold of 0.50 
#' was used: samples with estimated probability at least 0.50 were classified as 
#' non-normal, and samples below 0.50 were classified as normal. This choice 
#' avoids tuning the classification rule on the evaluation distributions and 
#' gives a transparent decision rule. Alternative thresholds may be examined as 
#' sensitivity analyses.

#' ============================================================
#' SECTION 22: PERFORMANCE METRICS
#' ============================================================

#' Compute classification metrics from prediction data.
#'
#' The positive class is "Non_Normal".
compute_metrics <- function(pred_df,
                            positive_class = "Non_Normal",
                            digits = NULL) {
  
  required_cols <- c("True_Class", "Predicted_Class", "Prob_Non_Normal")
  missing_cols  <- setdiff(required_cols, names(pred_df))
  
  if (length(missing_cols) > 0L) {
    stop(
      "pred_df is missing required column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  class_levels <- c("Non_Normal", "Normal")
  
  pred_df$True_Class <- factor(
    pred_df$True_Class,
    levels = class_levels
  )
  
  pred_df$Predicted_Class <- factor(
    pred_df$Predicted_Class,
    levels = class_levels
  )
  
  cm <- caret::confusionMatrix(
    data      = pred_df$Predicted_Class,
    reference = pred_df$True_Class,
    positive  = positive_class
  )
  
  accuracy    <- as.numeric(cm$overall["Accuracy"])
  sensitivity <- as.numeric(cm$byClass["Sensitivity"])
  specificity <- as.numeric(cm$byClass["Specificity"])
  precision   <- as.numeric(cm$byClass["Precision"])
  f1          <- as.numeric(cm$byClass["F1"])
  
  auc <- NA_real_
  
  if (length(unique(pred_df$True_Class)) == 2L &&
      all(is.finite(pred_df$Prob_Non_Normal))) {
    
    auc <- safe_calc(
      as.numeric(
        pROC::auc(
          response  = pred_df$True_Class,
          predictor = pred_df$Prob_Non_Normal,
          levels    = c("Normal", "Non_Normal"),
          direction = "<",
          quiet     = TRUE
        )
      ),
      default = NA_real_
    )
  }
  
  metric_values <- c(
    Accuracy    = accuracy,
    Sensitivity = sensitivity,
    Specificity = specificity,
    Precision   = precision,
    F1          = f1,
    AUC         = auc
  )
  
  if (!is.null(digits)) {
    metric_values <- round(metric_values, digits)
  }
  
  list(
    ConfusionMatrix = cm,
    Metrics         = metric_values,
    Predictions     = pred_df
  )
}


#' Summarize metrics for all models at one sample size.
get_metrics_summary <- function(model_metrics,
                                digits = 4) {
  
  model_names <- names(model_metrics)
  
  summary_df <- do.call(
    rbind,
    lapply(model_names, function(model_name) {
      
      metrics <- model_metrics[[model_name]]$Metrics
      
      data.frame(
        Model       = model_name,
        Accuracy    = metrics["Accuracy"],
        Sensitivity = metrics["Sensitivity"],
        Specificity = metrics["Specificity"],
        Precision   = metrics["Precision"],
        F1          = metrics["F1"],
        AUC         = metrics["AUC"],
        stringsAsFactors = FALSE
      )
    })
  )
  
  rownames(summary_df) <- NULL
  
  numeric_cols <- setdiff(names(summary_df), "Model")
  summary_df[numeric_cols] <- lapply(
    summary_df[numeric_cols],
    function(z) round(as.numeric(z), digits)
  )
  
  summary_df
}


#' ============================================================
#' SECTION 23: EVALUATE MODELS ACROSS SAMPLE SIZES
#' ============================================================

#' Evaluate trained models across all sample sizes.
#'
#' Returns one list element per sample size. Each element contains:
#'   - metrics: per-model metric objects;
#'   - summary: compact metric table;
#'   - predictions: raw prediction data frames;
#'   - holdout_features: optional processed holdout feature matrix.
evaluate_models_all_n <- function(trained_models,
                                  eval_specs,
                                  n_iter = 1000,
                                  threshold = 0.50,
                                  majority_threshold = 0.50,
                                  standardize_sample = FALSE,
                                  center_by = NULL,
                                  store_holdout_features = FALSE,
                                  digits = 4) {
  
  if (length(trained_models) == 0L) {
    stop("trained_models is empty.")
  }
  
  results <- vector("list", length(trained_models))
  names(results) <- names(trained_models)
  
  for (n_key in names(trained_models)) {
    
    model_bundle <- trained_models[[n_key]]
    
    if (is.null(model_bundle$n)) {
      stop("trained_models[['", n_key, "']] is missing element 'n'.")
    }
    
    n <- model_bundle$n
    
    eval_raw <- evaluate_one_n(
      n                      = n,
      model_bundle           = model_bundle,
      eval_specs             = eval_specs,
      n_iter                 = n_iter,
      threshold              = threshold,
      majority_threshold     = majority_threshold,
      standardize_sample     = standardize_sample,
      center_by              = center_by,
      store_holdout_features = store_holdout_features
    )
    
    model_metrics <- lapply(
      eval_raw$predictions,
      compute_metrics
    )
    
    summary_df <- get_metrics_summary(
      model_metrics = model_metrics,
      digits        = digits
    )
    
    results[[n_key]] <- list(
      n                = n,
      metrics          = model_metrics,
      summary          = summary_df,
      predictions      = eval_raw$predictions,
      holdout_features = eval_raw$holdout_features
    )
    
    cat("\nPerformance summary for n =", n, "\n")
    print(summary_df)
  }
  
  results
}


#' ============================================================
#' PART 5: ROC PLOTTING AND CLASSICAL-VS-ML COMPARISON
#' ============================================================


#' ============================================================
#' SECTION 24: ROC HELPER FUNCTIONS
#' ============================================================

#' Compute AUC using the trapezoidal rule.
#'
#' This is mainly used for manually constructed ROC curves, such as
#' classical-test ROC curves computed over a grid of alpha cutoffs.
compute_auc <- function(fpr,
                        tpr,
                        ensure_endpoints = TRUE) {
  
  ok <- is.finite(fpr) & is.finite(tpr)
  fpr <- fpr[ok]
  tpr <- tpr[ok]
  
  if (length(fpr) < 2L) {
    return(NA_real_)
  }
  
  ord <- order(fpr, tpr)
  fpr <- fpr[ord]
  tpr <- tpr[ord]
  
  if (isTRUE(ensure_endpoints)) {
    
    if (fpr[1L] > 0 || tpr[1L] > 0) {
      fpr <- c(0, fpr)
      tpr <- c(0, tpr)
    }
    
    last <- length(fpr)
    
    if (fpr[last] < 1 || tpr[last] < 1) {
      fpr <- c(fpr, 1)
      tpr <- c(tpr, 1)
    }
  }
  
  sum(diff(fpr) * (head(tpr, -1L) + tail(tpr, -1L)) / 2)
}


#' Compute ROC curve and AUC from ML prediction data.
#'
#' The positive class is "Non_Normal". Larger values of
#' Prob_Non_Normal indicate stronger evidence against Normality.
get_ml_roc <- function(pred_df) {
  
  required_cols <- c("True_Class", "Prob_Non_Normal")
  missing_cols  <- setdiff(required_cols, names(pred_df))
  
  if (length(missing_cols) > 0L) {
    stop(
      "pred_df is missing required column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  roc_obj <- pROC::roc(
    response = factor(
      pred_df$True_Class,
      levels = c("Normal", "Non_Normal")
    ),
    predictor = pred_df$Prob_Non_Normal,
    levels    = c("Normal", "Non_Normal"),
    direction = "<",
    quiet     = TRUE
  )
  
  list(
    fpr = 1 - roc_obj$specificities,
    tpr = roc_obj$sensitivities,
    auc = as.numeric(pROC::auc(roc_obj))
  )
}


#' Select plotting color.
get_plot_color <- function(name, color_map = NULL, default = "gray40") {
  
  if (!is.null(color_map) && name %in% names(color_map)) {
    return(color_map[[name]])
  }
  
  default
}


#' ============================================================
#' SECTION 25: ML-ONLY ROC PLOTS
#' ============================================================

#' Plot ML ROC curves for one sample size.
#'
#' This function uses stored prediction results. It does not rerun
#' simulations.
plot_ml_roc_one_n <- function(eval_results_n,
                              n,
                              ml_colors = NULL,
                              main_title = NULL) {
  
  if (is.null(eval_results_n$predictions)) {
    stop("eval_results_n must contain a 'predictions' element.")
  }
  
  prediction_list <- eval_results_n$predictions
  model_names     <- names(prediction_list)
  
  if (length(model_names) == 0L) {
    stop("No model prediction data found.")
  }
  
  if (is.null(main_title)) {
    main_title <- paste0("ML ROC Curves (n = ", n, ")")
  }
  
  plot(
    NA,
    xlim = c(0, 1),
    ylim = c(0, 1),
    xlab = "False Positive Rate",
    ylab = "True Positive Rate",
    main = main_title,
    las  = 1
  )
  
  abline(0, 1, lty = 2, col = "gray60")
  
  legend_labels <- character(0)
  legend_cols   <- character(0)
  
  for (model_name in model_names) {
    
    roc_vals <- get_ml_roc(prediction_list[[model_name]])
    col_m    <- get_plot_color(model_name, ml_colors)
    
    lines(
      roc_vals$fpr,
      roc_vals$tpr,
      col = col_m,
      lwd = 2
    )
    
    legend_labels <- c(
      legend_labels,
      sprintf("%s (AUC = %.3f)", model_name, roc_vals$auc)
    )
    
    legend_cols <- c(legend_cols, col_m)
  }
  
  legend(
    "bottomright",
    legend = legend_labels,
    col    = legend_cols,
    lwd    = 2,
    bty    = "o",
    cex    = 0.85
  )
  
  invisible(NULL)
}


#' Save ML ROC plots for all sample sizes.
save_ml_roc_plots <- function(eval_results_by_n,
                              sample_sizes = names(eval_results_by_n),
                              output_dir = "results",
                              ml_colors = NULL) {
  
  make_output_dir(output_dir)
  
  for (n in sample_sizes) {
    
    n_key <- as.character(n)
    
    if (!n_key %in% names(eval_results_by_n)) {
      warning("Skipping n = ", n, ": no evaluation result found.")
      next
    }
    
    pdf_file <- file.path(
      output_dir,
      paste0("ml_roc_n", n_key, ".pdf")
    )
    
    pdf(pdf_file, width = 6, height = 6)
    
    plot_ml_roc_one_n(
      eval_results_n = eval_results_by_n[[n_key]],
      n              = n_key,
      ml_colors      = ml_colors
    )
    
    dev.off()
    
    cat("Saved ML ROC plot:", pdf_file, "\n")
  }
  
  invisible(NULL)
}


#' ============================================================
#' SECTION 26: CLASSICAL NORMALITY-TEST ROC CURVES
#' ============================================================

#' Compute one classical normality-test p-value.
#'
#' Small p-values indicate evidence against Normality.
run_classical_test <- function(x,
                               test,
                               normal_par = NULL) {
  
  test <- toupper(trimws(test))
  
  if (test == "SW") {
    return(safe_calc(shapiro.test(x)$p.value))
  }
  
  if (test == "SF") {
    return(safe_calc(nortest::sf.test(x)$p.value))
  }
  
  if (test == "AD") {
    return(safe_calc(nortest::ad.test(x)$p.value))
  }
  
  if (test == "LF") {
    return(safe_calc(nortest::lillie.test(x)$p.value))
  }
  
  if (test == "JB") {
    return(safe_calc(tseries::jarque.bera.test(x)$p.value))
  }
  
  if (test == "KS") {
    
    if (is.null(normal_par)) {
      mu <- mean(x)
      s  <- sd(x)
    } else {
      mu <- normal_par[1L]
      s  <- normal_par[2L]
    }
    
    if (!is.finite(s) || s <= 0) {
      return(NA_real_)
    }
    
    return(
      safe_calc(
        ks.test(x, "pnorm", mean = mu, sd = s)$p.value
      )
    )
  }
  
  stop("Unsupported classical test: ", test)
}


#' Compute p-values for one classical test under Normal and alternative samples.
simulate_classical_pvalues <- function(n,
                                       normal_spec,
                                       alt_spec,
                                       test,
                                       n_sim = 1000) {
  
  p_norm <- numeric(n_sim)
  p_alt  <- numeric(n_sim)
  
  for (i in seq_len(n_sim)) {
    
    x_norm <- generate_data(
      n    = n,
      dist = normal_spec$dist,
      par  = normal_spec$par
    )
    
    x_alt <- generate_data(
      n    = n,
      dist = alt_spec$dist,
      par  = alt_spec$par
    )
    
    p_norm[i] <- run_classical_test(
      x          = x_norm,
      test       = test,
      normal_par = normal_spec$par
    )
    
    p_alt[i] <- run_classical_test(
      x          = x_alt,
      test       = test,
      normal_par = normal_spec$par
    )
  }
  
  list(
    p_norm = p_norm,
    p_alt  = p_alt
  )
}


#' Convert Normal and alternative p-values into ROC coordinates.
pvalues_to_roc <- function(p_norm,
                           p_alt,
                           alpha_grid) {
  
  fpr <- vapply(
    alpha_grid,
    function(alpha) mean(p_norm < alpha, na.rm = TRUE),
    numeric(1L)
  )
  
  tpr <- vapply(
    alpha_grid,
    function(alpha) mean(p_alt < alpha, na.rm = TRUE),
    numeric(1L)
  )
  
  list(
    fpr = fpr,
    tpr = tpr
  )
}


#' Compute classical ROC curves for one Normal/alternative pair.
#'
#' For each classical normality test:
#'   FPR = P(reject Normality | matched Normal)
#'   TPR = P(reject Normality | alternative)
#'
#' Rejection rule:
#'   p-value < alpha.
compute_classical_roc <- function(n,
                                  normal_spec,
                                  alt_spec,
                                  tests = c("SW", "AD", "JB"),
                                  alpha_grid = seq(0, 1, by = 0.05),
                                  n_sim = 1000) {
  
  tests <- toupper(trimws(tests))
  
  fpr <- matrix(
    NA_real_,
    nrow = length(tests),
    ncol = length(alpha_grid),
    dimnames = list(tests, NULL)
  )
  
  tpr <- fpr
  
  for (test_name in tests) {
    
    pvals <- simulate_classical_pvalues(
      n           = n,
      normal_spec = normal_spec,
      alt_spec    = alt_spec,
      test        = test_name,
      n_sim       = n_sim
    )
    
    roc_vals <- pvalues_to_roc(
      p_norm     = pvals$p_norm,
      p_alt      = pvals$p_alt,
      alpha_grid = alpha_grid
    )
    
    fpr[test_name, ] <- roc_vals$fpr
    tpr[test_name, ] <- roc_vals$tpr
  }
  
  list(
    FPR   = fpr,
    TPR   = tpr,
    alpha = alpha_grid,
    tests = tests
  )
}


#' ============================================================
#' SECTION 27: ML ROC CURVES FOR ONE DISTRIBUTION PAIR
#' ============================================================

#' Generate feature rows for one distribution pair.
generate_pair_prediction_features <- function(n,
                                              normal_spec,
                                              alt_spec,
                                              model_bundle,
                                              n_sim = 1000) {
  
  specs <- list(
    list(spec = normal_spec, label = "Normal"),
    list(spec = alt_spec,    label = "Non_Normal")
  )
  
  rows        <- vector("list", 2L * n_sim)
  true_class  <- character(2L * n_sim)
  row_id      <- 0L
  
  for (item in specs) {
    
    for (i in seq_len(n_sim)) {
      
      row_id <- row_id + 1L
      
      x <- generate_data(
        n    = n,
        dist = item$spec$dist,
        par  = item$spec$par
      )
      
      rows[[row_id]] <- calculate_features(
        x           = x,
        feature_set = model_bundle$feature_names
      )
      
      true_class[row_id] <- item$label
    }
  }
  
  feature_df <- do.call(rbind, rows)
  rownames(feature_df) <- NULL
  
  list(
    features   = feature_df,
    true_class = true_class
  )
}


#' Compute ML probabilities for one Normal/alternative pair.
compute_ml_pair_probabilities <- function(n,
                                          normal_spec,
                                          alt_spec,
                                          model_bundle,
                                          ml_methods = names(model_bundle$models),
                                          n_sim = 1000) {
  
  ml_methods <- intersect(ml_methods, names(model_bundle$models))
  
  if (length(ml_methods) == 0L) {
    stop("None of the requested ML methods are available in model_bundle.")
  }
  
  pair_data <- generate_pair_prediction_features(
    n            = n,
    normal_spec  = normal_spec,
    alt_spec     = alt_spec,
    model_bundle = model_bundle,
    n_sim        = n_sim
  )
  
  feats_std <- prepare_prediction_matrix(
    feature_df   = pair_data$features,
    model_bundle = model_bundle
  )
  
  all_probs <- predict_model_probs(
    model_bundle = model_bundle,
    feats_std    = feats_std
  )
  
  all_probs <- all_probs[, ml_methods, drop = FALSE]
  
  pred_store <- lapply(ml_methods, function(model_name) {
    
    data.frame(
      True_Class      = pair_data$true_class,
      Prob_Non_Normal = all_probs[[model_name]],
      stringsAsFactors = FALSE
    )
  })
  
  names(pred_store) <- ml_methods
  
  pred_store$MajorityVote <- data.frame(
    True_Class      = pair_data$true_class,
    Prob_Non_Normal = rowMeans(all_probs, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  
  pred_store
}


#' Compute ML ROC curves for one Normal/alternative pair.
#'
#' This function reruns pair-specific simulations and obtains classifier
#' probabilities. It is intended for direct classical-vs-ML ROC comparison
#' on the same distribution pair.
compute_ml_roc_pair <- function(n,
                                normal_spec,
                                alt_spec,
                                model_bundle,
                                ml_methods = names(model_bundle$models),
                                n_sim = 1000) {
  
  pred_store <- compute_ml_pair_probabilities(
    n            = n,
    normal_spec  = normal_spec,
    alt_spec     = alt_spec,
    model_bundle = model_bundle,
    ml_methods   = ml_methods,
    n_sim        = n_sim
  )
  
  lapply(pred_store, get_ml_roc)
}


#' ============================================================
#' SECTION 28: COMBINED CLASSICAL-VS-ML ROC PLOT
#' ============================================================

#' Plot classical and ML ROC curves in one panel.
#'
#' Classical tests are plotted with dashed lines.
#' ML classifiers are plotted with solid lines.
plot_classical_ml_roc <- function(classical_roc,
                                  ml_roc,
                                  main_title,
                                  test_colors = NULL,
                                  ml_colors = NULL,
                                  legend_cex = 0.80) {
  
  plot(
    NA,
    xlim = c(0, 1),
    ylim = c(0, 1),
    xlab = "False Positive Rate",
    ylab = "True Positive Rate",
    main = main_title,
    las  = 1
  )
  
  abline(0, 1, lty = 2, col = "gray60")
  
  legend_labels <- character(0)
  legend_cols   <- character(0)
  legend_lty    <- numeric(0)
  
  ## Classical ROC curves
  for (test_name in rownames(classical_roc$FPR)) {
    
    auc_val <- compute_auc(
      fpr = classical_roc$FPR[test_name, ],
      tpr = classical_roc$TPR[test_name, ]
    )
    
    col_t <- get_plot_color(test_name, test_colors)
    
    lines(
      classical_roc$FPR[test_name, ],
      classical_roc$TPR[test_name, ],
      col = col_t,
      lwd = 2,
      lty = 2
    )
    
    legend_labels <- c(
      legend_labels,
      sprintf("%s (AUC = %.3f)", test_name, auc_val)
    )
    
    legend_cols <- c(legend_cols, col_t)
    legend_lty  <- c(legend_lty, 2)
  }
  
  ## ML ROC curves
  for (model_name in names(ml_roc)) {
    
    col_m <- get_plot_color(model_name, ml_colors)
    
    lines(
      ml_roc[[model_name]]$fpr,
      ml_roc[[model_name]]$tpr,
      col = col_m,
      lwd = 2,
      lty = 1
    )
    
    legend_labels <- c(
      legend_labels,
      sprintf("%s (AUC = %.3f)", model_name, ml_roc[[model_name]]$auc)
    )
    
    legend_cols <- c(legend_cols, col_m)
    legend_lty  <- c(legend_lty, 1)
  }
  
  legend(
    "bottomright",
    legend = legend_labels,
    col    = legend_cols,
    lty    = legend_lty,
    lwd    = 2,
    bty    = "o",
    cex    = legend_cex
  )
  
  invisible(NULL)
}


#' ============================================================
#' SECTION 29: RUN ROC COMPARISON ACROSS PAIRS AND SAMPLE SIZES
#' ============================================================

#' Create a safe label for file names.
safe_dist_label <- function(spec) {
  
  label <- paste0(
    spec$dist,
    "_",
    paste(spec$par, collapse = "_")
  )
  
  label <- gsub("[^A-Za-z0-9_]+", "_", label)
  label
}


#' Run classical-vs-ML ROC comparison across distribution pairs.
#'
#' For each alternative distribution, this function saves one PDF.
#' Each PDF contains one panel per sample size.
run_roc_comparison <- function(trained_models,
                               roc_alt_specs,
                               sample_sizes = names(trained_models),
                               classical_tests = c("SW", "AD", "JB"),
                               ml_methods = names(trained_models[[1L]]$models),
                               alpha_grid = seq(0, 1, by = 0.05),
                               n_sim = 1000,
                               output_dir = "results",
                               test_colors = NULL,
                               ml_colors = NULL) {
  
  if (length(trained_models) == 0L) {
    stop("trained_models is empty.")
  }
  
  if (length(roc_alt_specs) == 0L) {
    stop("roc_alt_specs is empty.")
  }
  
  make_output_dir(output_dir)
  
  roc_pairs <- make_paired_specs(roc_alt_specs)
  
  roc_classical <- list()
  roc_ml        <- list()
  
  for (pair_id in seq_along(roc_pairs)) {
    
    pair     <- roc_pairs[[pair_id]]
    alt_name <- safe_dist_label(pair$alt)
    
    pdf_file <- file.path(
      output_dir,
      paste0("roc_comparison_", alt_name, ".pdf")
    )
    
    pdf(
      file   = pdf_file,
      width  = 6 * length(sample_sizes),
      height = 6
    )
    
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par), add = TRUE)
    
    par(
      mfrow = c(1, length(sample_sizes)),
      mar   = c(4, 4, 3, 1),
      oma   = c(1, 1, 2, 1)
    )
    
    for (n in sample_sizes) {
      
      n_key <- as.character(n)
      
      if (!n_key %in% names(trained_models)) {
        warning("Skipping n = ", n, ": no trained model found.")
        next
      }
      
      cat("ROC comparison:", pair$alt$dist, "| n =", n_key, "\n")
      
      classical_res <- compute_classical_roc(
        n           = as.numeric(n_key),
        normal_spec = pair$normal,
        alt_spec    = pair$alt,
        tests       = classical_tests,
        alpha_grid  = alpha_grid,
        n_sim       = n_sim
      )
      
      ml_res <- compute_ml_roc_pair(
        n            = as.numeric(n_key),
        normal_spec  = pair$normal,
        alt_spec     = pair$alt,
        model_bundle = trained_models[[n_key]],
        ml_methods   = ml_methods,
        n_sim        = n_sim
      )
      
      result_key <- paste0(alt_name, "_n", n_key)
      
      roc_classical[[result_key]] <- classical_res
      roc_ml[[result_key]]        <- ml_res
      
      title_text <- paste0(
        tools::toTitleCase(gsub("_", " ", pair$alt$dist)),
        " vs Matched Normal | n = ",
        n_key
      )
      
      plot_classical_ml_roc(
        classical_roc = classical_res,
        ml_roc        = ml_res,
        main_title    = title_text,
        test_colors   = test_colors,
        ml_colors     = ml_colors
      )
    }
    
    mtext(
      paste0(
        "ROC Comparison: ",
        tools::toTitleCase(gsub("_", " ", pair$alt$dist)),
        " vs Moment-Matched Normal"
      ),
      outer = TRUE,
      cex   = 1.1,
      font  = 2
    )
    
    dev.off()
    
    cat("Saved ROC comparison:", pdf_file, "\n")
  }
  
  invisible(
    list(
      classical = roc_classical,
      ml        = roc_ml
    )
  )
}


#' ============================================================
#' PART 6: PERMUTATION AUC VARIABLE IMPORTANCE
#' ============================================================


#' ============================================================
#' SECTION 30: PERMUTATION VIP HELPERS
#' ============================================================

#' Compute AUC from true labels and predicted probabilities.
#'
#' The positive class is "Non_Normal".
compute_prob_auc <- function(true_class,
                             prob_non_normal) {
  
  if (length(unique(true_class)) < 2L) {
    return(NA_real_)
  }
  
  if (!all(is.finite(prob_non_normal))) {
    return(NA_real_)
  }
  
  safe_calc(
    as.numeric(
      pROC::auc(
        pROC::roc(
          response = factor(
            true_class,
            levels = c("Normal", "Non_Normal")
          ),
          predictor = prob_non_normal,
          levels    = c("Normal", "Non_Normal"),
          direction = "<",
          quiet     = TRUE
        )
      )
    ),
    default = NA_real_
  )
}


#' Extract held-out labels from one sample-size evaluation object.
#'
#' All classifiers are evaluated on the same held-out samples, so
#' the labels can be taken from the first prediction data frame.
get_holdout_labels <- function(eval_results_n) {
  
  if (is.null(eval_results_n$predictions)) {
    stop("eval_results_n must contain a 'predictions' element.")
  }
  
  prediction_list <- eval_results_n$predictions
  
  if (length(prediction_list) == 0L) {
    stop("No prediction data found in eval_results_n$predictions.")
  }
  
  first_model <- names(prediction_list)[1L]
  
  prediction_list[[first_model]]$True_Class
}


#' Predict P(Non_Normal) for one trained model.
predict_one_model_prob <- function(model,
                                   new_data) {
  
  pred <- predict(
    model,
    newdata = new_data,
    type    = "prob"
  )
  
  if (!"Non_Normal" %in% names(pred)) {
    stop("Model prediction does not contain a 'Non_Normal' probability column.")
  }
  
  as.numeric(pred[, "Non_Normal"])
}



#' ============================================================
#' SECTION 31: PERMUTATION VIP FOR ONE MODEL
#' ============================================================

#' Compute permutation AUC-drop importance for one model.
#'
#' For each feature:
#'   1. permute the feature in the held-out data;
#'   2. recompute predicted probabilities;
#'   3. recompute AUC;
#'   4. record baseline AUC - permuted AUC.
#'
#' Larger AUC drops indicate greater model reliance on the feature.
compute_vip_one_model <- function(model,
                                  model_name,
                                  holdout_x,
                                  holdout_y,
                                  n_permutations = 20,
                                  digits = 4) {
  
  if (ncol(holdout_x) == 0L) {
    stop("holdout_x has no feature columns.")
  }
  
  feature_names <- names(holdout_x)
  
  base_prob <- predict_one_model_prob(
    model    = model,
    new_data = holdout_x
  )
  
  baseline_auc <- compute_prob_auc(
    true_class      = holdout_y,
    prob_non_normal = base_prob
  )
  
  cat(
    sprintf(
      "  %s baseline AUC = %.4f | %d features x %d permutations\n",
      model_name,
      baseline_auc,
      length(feature_names),
      n_permutations
    )
  )
  
  auc_drop <- setNames(
    numeric(length(feature_names)),
    feature_names
  )
  
  for (feature in feature_names) {
    
    perm_auc <- numeric(n_permutations)
    
    for (b in seq_len(n_permutations)) {
      
      perm_x <- holdout_x
      perm_x[[feature]] <- sample(perm_x[[feature]])
      
      perm_prob <- safe_calc(
        predict_one_model_prob(
          model    = model,
          new_data = perm_x
        ),
        default = rep(NA_real_, nrow(perm_x))
      )
      
      perm_auc[b] <- compute_prob_auc(
        true_class      = holdout_y,
        prob_non_normal = perm_prob
      )
    }
    
    auc_drop[feature] <- baseline_auc - mean(perm_auc, na.rm = TRUE)
  }
  
  imp_df <- data.frame(
    Feature  = names(auc_drop),
    AUC_Drop = round(as.numeric(auc_drop), digits),
    stringsAsFactors = FALSE
  )
  
  imp_df <- imp_df[
    order(imp_df$AUC_Drop, decreasing = TRUE),
    ,
    drop = FALSE
  ]
  
  imp_df$Rank <- seq_len(nrow(imp_df))
  imp_df <- imp_df[, c("Rank", "Feature", "AUC_Drop")]
  rownames(imp_df) <- NULL
  
  list(
    model_name   = model_name,
    baseline_auc = baseline_auc,
    imp_df       = imp_df
  )
}


#' ============================================================
#' SECTION 32: PERMUTATION VIP FOR ALL MODELS AND SAMPLE SIZES
#' ============================================================

#' Compute permutation AUC-drop importance for selected models.
#'
#' This requires evaluate_models_all_n(..., store_holdout_features = TRUE).
#' The function is optional because it can be time-consuming.
compute_permutation_vip <- function(trained_models,
                                    eval_results_by_n,
                                    n_permutations = 20,
                                    models_for_vip = NULL,
                                    digits = 4) {
  
  if (length(trained_models) == 0L) {
    stop("trained_models is empty.")
  }
  
  vip_results <- vector("list", length(trained_models))
  names(vip_results) <- names(trained_models)
  
  for (n_key in names(trained_models)) {
    
    cat("\nPermutation VIP for n =", n_key, "\n")
    
    if (!n_key %in% names(eval_results_by_n)) {
      stop("No evaluation results found for n = ", n_key)
    }
    
    model_bundle <- trained_models[[n_key]]
    eval_n       <- eval_results_by_n[[n_key]]
    
    if (is.null(eval_n$holdout_features)) {
      stop(
        "Missing holdout_features for n = ",
        n_key,
        ". Rerun evaluate_models_all_n(..., store_holdout_features = TRUE)."
      )
    }
    
    holdout_x <- eval_n$holdout_features
    holdout_y <- get_holdout_labels(eval_n)
    
    model_names <- names(model_bundle$models)
    
    if (!is.null(models_for_vip)) {
      model_names <- intersect(model_names, models_for_vip)
    }
    
    if (length(model_names) == 0L) {
      stop("No requested VIP models found for n = ", n_key)
    }
    
    vip_results[[n_key]] <- vector("list", length(model_names))
    names(vip_results[[n_key]]) <- model_names
    
    for (model_name in model_names) {
      
      vip_results[[n_key]][[model_name]] <- compute_vip_one_model(
        model          = model_bundle$models[[model_name]],
        model_name     = model_name,
        holdout_x      = holdout_x,
        holdout_y      = holdout_y,
        n_permutations = n_permutations,
        digits         = digits
      )
      
      vip_results[[n_key]][[model_name]]$n <- model_bundle$n
    }
  }
  
  vip_results
}


#' ============================================================
#' SECTION 33: SUMMARIZE AND SELECT PERMUTATION VIP FEATURES
#' ============================================================

#' Convert permutation VIP results to long format.
vip_results_to_long <- function(vip_results) {
  
  rows <- list()
  
  for (n_key in names(vip_results)) {
    
    vip_n <- vip_results[[n_key]]
    
    if (is.null(vip_n)) {
      next
    }
    
    for (model_name in names(vip_n)) {
      
      entry <- vip_n[[model_name]]
      
      if (is.null(entry)) {
        next
      }
      
      rows[[length(rows) + 1L]] <- data.frame(
        n        = entry$n,
        Model    = model_name,
        Feature  = entry$imp_df$Feature,
        AUC_Drop = entry$imp_df$AUC_Drop,
        stringsAsFactors = FALSE
      )
    }
  }
  
  if (length(rows) == 0L) {
    return(NULL)
  }
  
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  
  out
}


#' Summarize permutation VIP across models and sample sizes.
summarize_permutation_vip <- function(vip_results,
                                      digits = 4) {
  
  vip_long <- vip_results_to_long(vip_results)
  
  if (is.null(vip_long) || nrow(vip_long) == 0L) {
    return(NULL)
  }
  
  split_rows <- split(vip_long, vip_long$Feature)
  
  summary_df <- do.call(
    rbind,
    lapply(names(split_rows), function(feature) {
      
      z <- split_rows[[feature]]$AUC_Drop
      
      data.frame(
        Feature         = feature,
        Mean_AUC_Drop   = round(mean(z, na.rm = TRUE), digits),
        Median_AUC_Drop = round(median(z, na.rm = TRUE), digits),
        Min_AUC_Drop    = round(min(z, na.rm = TRUE), digits),
        Max_AUC_Drop    = round(max(z, na.rm = TRUE), digits),
        Pct_Positive    = round(mean(z > 0, na.rm = TRUE), digits),
        stringsAsFactors = FALSE
      )
    })
  )
  
  summary_df <- summary_df[
    order(summary_df$Mean_AUC_Drop, decreasing = TRUE),
    ,
    drop = FALSE
  ]
  
  rownames(summary_df) <- NULL
  summary_df
}


#' Select top features from a permutation VIP summary table.
select_top_vip_features <- function(vip_summary,
                                    top_n = 30,
                                    rank_by = c("mean", "median", "min"),
                                    min_pct_positive = 0,
                                    require_positive_mean = TRUE) {
  
  rank_by <- match.arg(rank_by)
  
  if (is.null(vip_summary) || nrow(vip_summary) == 0L) {
    stop("vip_summary is empty.")
  }
  
  needed_cols <- c(
    "Feature",
    "Mean_AUC_Drop",
    "Median_AUC_Drop",
    "Min_AUC_Drop",
    "Pct_Positive"
  )
  
  missing_cols <- setdiff(needed_cols, names(vip_summary))
  
  if (length(missing_cols) > 0L) {
    stop(
      "vip_summary is missing required column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  vip_keep <- vip_summary[
    vip_summary$Pct_Positive >= min_pct_positive,
    ,
    drop = FALSE
  ]
  
  if (isTRUE(require_positive_mean)) {
    vip_keep <- vip_keep[
      vip_keep$Mean_AUC_Drop > 0,
      ,
      drop = FALSE
    ]
  }
  
  if (nrow(vip_keep) == 0L) {
    warning("No features met the selection criteria.")
    return(character(0))
  }
  
  rank_col <- switch(
    rank_by,
    mean   = "Mean_AUC_Drop",
    median = "Median_AUC_Drop",
    min    = "Min_AUC_Drop"
  )
  
  vip_keep <- vip_keep[
    order(vip_keep[[rank_col]], decreasing = TRUE),
    ,
    drop = FALSE
  ]
  
  head(vip_keep$Feature, top_n)
}


#' Return selected features together with their VIP metrics.
get_selected_feature_metrics <- function(vip_summary,
                                         selected_features) {
  
  out <- vip_summary[
    match(selected_features, vip_summary$Feature),
    ,
    drop = FALSE
  ]
  
  out <- data.frame(
    Rank = seq_len(nrow(out)),
    out,
    row.names = NULL
  )
  
  out
}


#' ============================================================
#' SECTION 34: PERMUTATION VIP PLOTS
#' ============================================================

#' Plot permutation VIP for one model and one sample size.
plot_permutation_vip <- function(vip_entry,
                                 top_n = 30,
                                 main_title = NULL) {
  
  if (is.null(vip_entry)) {
    stop("vip_entry is NULL.")
  }
  
  imp_df <- head(vip_entry$imp_df, top_n)
  
  if (nrow(imp_df) == 0L) {
    stop("vip_entry$imp_df is empty.")
  }
  
  ## Reverse order so the largest AUC drop appears at the top.
  imp_df <- imp_df[rev(seq_len(nrow(imp_df))), , drop = FALSE]
  
  if (is.null(main_title)) {
    main_title <- paste0(
      "Permutation VIP: ",
      vip_entry$model_name,
      " (n = ",
      vip_entry$n,
      ")"
    )
  }
  
  old_mar <- par("mar")
  on.exit(par(mar = old_mar), add = TRUE)
  
  par(mar = c(5, 12, 4, 2))
  
  barplot(
    height    = imp_df$AUC_Drop,
    names.arg = imp_df$Feature,
    horiz     = TRUE,
    las       = 1,
    border    = NA,
    xlab      = "AUC Drop After Permutation",
    main      = main_title,
    cex.names = 0.80,
    cex.main  = 0.95
  )
  
  abline(v = 0, lty = 2, col = "gray50")
  
  invisible(NULL)
}


#' Save VIP plots for selected models and sample sizes.
save_permutation_vip_plots <- function(vip_results,
                                       sample_sizes = names(vip_results),
                                       models = c("RF"),
                                       top_n = 30,
                                       output_dir = "results") {
  
  make_output_dir(output_dir)
  
  for (n in sample_sizes) {
    
    n_key <- as.character(n)
    
    if (!n_key %in% names(vip_results)) {
      warning("Skipping n = ", n_key, ": no VIP result found.")
      next
    }
    
    for (model_name in models) {
      
      vip_entry <- vip_results[[n_key]][[model_name]]
      
      if (is.null(vip_entry)) {
        next
      }
      
      pdf_file <- file.path(
        output_dir,
        paste0("vip_perm_", model_name, "_n", n_key, ".pdf")
      )
      
      pdf(pdf_file, width = 8, height = 7)
      
      plot_permutation_vip(
        vip_entry = vip_entry,
        top_n     = top_n
      )
      
      dev.off()
      
      cat("Saved permutation VIP plot:", pdf_file, "\n")
    }
  }
  
  invisible(NULL)
}


#' ============================================================
#' PART 7: DENSITY PLOTS FOR MOMENT-MATCHED PAIRS
#' ============================================================


#' ============================================================
#' SECTION 35: DENSITY HELPER FUNCTIONS
#' ============================================================

#' Evaluate the theoretical density for a distribution specification.
#'
#' Used only for plotting.
#'
#' Required structure:
#'   spec$dist
#'   spec$par
density_value <- function(x_grid, spec) {
  
  dist <- tolower(trimws(spec$dist))
  par  <- spec$par
  
  out <- switch(
    dist,
    
    normal = dnorm(
      x_grid,
      mean = par[1L],
      sd   = par[2L]
    ),
    
    gumbel = evd::dgumbel(
      x_grid,
      loc   = par[1L],
      scale = par[2L]
    ),
    
    chi_square = ifelse(
      x_grid >= 0,
      dchisq(x_grid, df = par[1L]),
      0
    ),
    
    gamma = ifelse(
      x_grid >= 0,
      dgamma(x_grid, shape = par[1L], rate = par[2L]),
      0
    ),
    
    exponential = ifelse(
      x_grid >= 0,
      dexp(x_grid, rate = par[1L]),
      0
    ),
    
    weibull = ifelse(
      x_grid >= 0,
      dweibull(x_grid, shape = par[1L], scale = par[2L]),
      0
    ),
    
    laplace = LaplacesDemon::dlaplace(
      x_grid,
      location = par[1L],
      scale    = par[2L]
    ),
    
    beta = ifelse(
      x_grid >= 0 & x_grid <= 1,
      dbeta(x_grid, shape1 = par[1L], shape2 = par[2L]),
      0
    ),
    
    uniform = dunif(
      x_grid,
      min = par[1L],
      max = par[2L]
    ),
    
    logistic = dlogis(
      x_grid,
      location = par[1L],
      scale    = par[2L]
    ),
    
    lognormal = ifelse(
      x_grid > 0,
      dlnorm(x_grid, meanlog = par[1L], sdlog = par[2L]),
      0
    ),
    
    t = dt(
      x_grid,
      df = par[1L]
    ),
    
    f = ifelse(
      x_grid >= 0,
      df(x_grid, df1 = par[1L], df2 = par[2L]),
      0
    ),
    
    pareto = VGAM::dpareto(
      x_grid,
      shape = par[1L],
      scale = par[2L]
    ),
    
    cauchy = dcauchy(
      x_grid,
      location = par[1L],
      scale    = par[2L]
    ),
    
    contaminated = {
      p   <- par[1L]
      mu  <- par[2L]
      sd1 <- par[3L]
      sd2 <- par[4L]
      
      p * dnorm(x_grid, mean = mu, sd = sd1) +
        (1 - p) * dnorm(x_grid, mean = mu, sd = sd2)
    },
    
    stop("Unsupported distribution in density_value(): ", dist)
  )
  
  out[!is.finite(out)] <- NA_real_
  out
}


#' Choose an x-axis range for a moment-matched density comparison.
#'
#' The range combines:
#'   1. a distribution-specific plotting range for the alternative;
#'   2. mean +/- 4 SD for the matched Normal.
density_xlim <- function(alt_spec,
                         normal_spec) {
  
  dist <- tolower(trimws(alt_spec$dist))
  par  <- alt_spec$par
  
  normal_mu <- normal_spec$par[1L]
  normal_sd <- normal_spec$par[2L]
  
  alt_range <- switch(
    dist,
    
    gumbel = c(
      par[1L] - 3 * par[2L],
      par[1L] + 7 * par[2L]
    ),
    
    chi_square = c(
      0,
      qchisq(0.999, df = par[1L])
    ),
    
    gamma = c(
      0,
      qgamma(0.999, shape = par[1L], rate = par[2L])
    ),
    
    exponential = c(
      0,
      qexp(0.999, rate = par[1L])
    ),
    
    weibull = c(
      0,
      qweibull(0.999, shape = par[1L], scale = par[2L])
    ),
    
    laplace = c(
      par[1L] - 8 * par[2L],
      par[1L] + 8 * par[2L]
    ),
    
    beta = c(-0.05, 1.05),
    
    uniform = c(
      par[1L] - 0.05,
      par[2L] + 0.05
    ),
    
    logistic = qlogis(
      c(0.001, 0.999),
      location = par[1L],
      scale    = par[2L]
    ),
    
    lognormal = c(
      0,
      qlnorm(0.999, meanlog = par[1L], sdlog = par[2L])
    ),
    
    t = qt(
      c(0.001, 0.999),
      df = par[1L]
    ),
    
    f = c(
      0,
      qf(0.999, df1 = par[1L], df2 = par[2L])
    ),
    
    pareto = c(
      par[2L],
      VGAM::qpareto(0.999, shape = par[1L], scale = par[2L])
    ),
    
    cauchy = c(
      par[1L] - 8 * par[2L],
      par[1L] + 8 * par[2L]
    ),
    
    contaminated = c(
      par[2L] - 4 * par[4L],
      par[2L] + 4 * par[4L]
    ),
    
    stop("Unsupported distribution in density_xlim(): ", dist)
  )
  
  normal_range <- c(
    normal_mu - 4 * normal_sd,
    normal_mu + 4 * normal_sd
  )
  
  xlim <- range(c(alt_range, normal_range), finite = TRUE)
  
  if (length(xlim) != 2L || !all(is.finite(xlim)) || xlim[1L] >= xlim[2L]) {
    stop("Could not construct a valid x-axis range for ", dist)
  }
  
  xlim
}


#' Create a readable distribution label.
dist_label <- function(spec, digits = 3) {
  
  dist_name <- tools::toTitleCase(gsub("_", " ", spec$dist))
  par_text  <- paste(round(spec$par, digits), collapse = ", ")
  
  paste0(dist_name, "(", par_text, ")")
}



#' ============================================================
#' SECTION 36: DENSITY PLOT FUNCTION
#' ============================================================

#' Plot alternative distributions against their moment-matched Normals.
#'
#' Each panel compares one alternative distribution to its own
#' moment-matched Normal distribution.
#'
#' Input:
#'   paired_specs from make_paired_specs().
#'
#' Output:
#'   A saved PDF file.
plot_density_pairs <- function(paired_specs,
                               output_dir = "results",
                               file_name = "density_pairs.pdf",
                               n_grid = 2000,
                               n_col = 4,
                               alt_col = "brown",
                               normal_col = "steelblue") {
  
  if (length(paired_specs) == 0L) {
    stop("paired_specs is empty.")
  }
  
  if (n_grid < 100L) {
    stop("n_grid should be at least 100.")
  }
  
  make_output_dir(output_dir)
  
  n_pairs <- length(paired_specs)
  n_col   <- min(n_col, n_pairs)
  n_row   <- ceiling(n_pairs / n_col)
  
  pdf_file <- file.path(output_dir, file_name)
  
  pdf(
    file   = pdf_file,
    width  = 4.2 * n_col,
    height = 4.0 * n_row
  )
  
  old_par <- par(no.readonly = TRUE)
  
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)
  
  par(
    mfrow = c(n_row, n_col),
    mar   = c(4.2, 4.2, 3.0, 1.2),
    oma   = c(1.0, 1.0, 3.0, 0.5),
    mgp   = c(2.5, 0.7, 0)
  )
  
  for (pair in paired_specs) {
    
    normal_spec <- pair$normal
    alt_spec    <- pair$alt
    
    xlim <- density_xlim(
      alt_spec    = alt_spec,
      normal_spec = normal_spec
    )
    
    x_grid <- seq(
      from       = xlim[1L],
      to         = xlim[2L],
      length.out = n_grid
    )
    
    d_alt <- density_value(
      x_grid = x_grid,
      spec   = alt_spec
    )
    
    d_norm <- density_value(
      x_grid = x_grid,
      spec   = normal_spec
    )
    
    ymax <- max(c(d_alt, d_norm), na.rm = TRUE)
    
    if (!is.finite(ymax) || ymax <= 0) {
      ymax <- 1
    }
    
    plot(
      NA,
      xlim = xlim,
      ylim = c(0, 1.10 * ymax),
      xlab = "x",
      ylab = "Density",
      main = paste0(
        tools::toTitleCase(gsub("_", " ", alt_spec$dist)),
        " vs Matched Normal"
      ),
      las      = 1,
      cex.main = 1.2
    )
    
    polygon(
      x      = c(x_grid, rev(x_grid)),
      y      = c(d_alt, rep(0, length(x_grid))),
      col    = adjustcolor(alt_col, alpha.f = 0.18),
      border = NA
    )
    
    polygon(
      x      = c(x_grid, rev(x_grid)),
      y      = c(d_norm, rep(0, length(x_grid))),
      col    = adjustcolor(normal_col, alpha.f = 0.18),
      border = NA
    )
    
    lines(
      x_grid,
      d_alt,
      col = alt_col,
      lwd = 3
    )
    
    lines(
      x_grid,
      d_norm,
      col = normal_col,
      lwd = 3,
      lty = 2
    )
    
    abline(
      v   = normal_spec$par[1L],
      col = "gray50",
      lty = 3,
      lwd = 1.5
    )
    
    legend(
      "topright",
      legend = c(
        dist_label(alt_spec),
        dist_label(normal_spec)
      ),
      col = c(
        alt_col,
        normal_col
      ),
      lwd = 3,
      lty = c(1, 2),
      bty = "o",
      cex = 0.95
    )
  }
  
  mtext(
    "Alternative Distributions vs Moment-Matched Normal Distributions",
    outer = TRUE,
    cex   = 1.2,
    font  = 3,
    line  = 1.0
  )
  
  cat("Saved density plots:", pdf_file, "\n")
  
  invisible(pdf_file)
}


#' ============================================================
#' PART 8: MAIN RUN FUNCTION
#' ============================================================


#' ============================================================
#' SECTION 37: MAIN PIPELINE WRAPPER
#' ============================================================

#' Run the ML-based normality testing framework.
#'
#' This is the main entry point for the full workflow.
#'
#' The pipeline can:
#'   1. train or load ML models;
#'   2. evaluate models on held-out distributions using a fixed threshold;
#'   3. save ML ROC plots;
#'   4. optionally compute permutation AUC VIP;
#'   5. optionally create density plots;
#'   6. optionally run classical-vs-ML ROC comparisons.
#'
#' Expensive steps are optional so the framework can be debugged quickly.
run_normality_framework <- function(
    
  ## Core simulation settings
  sample_sizes = c(10, 50),
  num_sim      = 100,
  n_iter_eval  = 1000,
  
  ## Models
  models_to_train = c("LR", "RF", "ANN", "GBM", "SVM", "KNN"),
  
  ## Fixed probability thresholds
  threshold          = 0.50,
  majority_threshold = 0.50,
  
  ## Distribution settings
  train_alt_specs = default_training_specs(),
  eval_alt_specs  = default_eval_specs(),
  roc_alt_specs   = default_eval_specs(),
  
  ## Feature settings
  feature_set          = NULL,
  standardize_sample   = FALSE,
  center_by            = NULL,
  
  ## Model training/loading
  train_models    = TRUE,
  model_save_path = NULL,
  k_folds         = 10,
  
  ## Optional steps
  make_density_plots  = TRUE,
  run_permutation_vip = TRUE,
  run_roc_analysis    = TRUE,
  
  ## Permutation VIP settings
  n_permutations = 20,
  models_for_vip = c("RF"),
  vip_top_n      = 30,
  vip_rank_by    = c("mean", "median", "min"),
  vip_min_pct_positive      = 0,
  vip_require_positive_mean = TRUE,
  
  ## ROC settings
  classical_tests = c("SW", "AD", "JB"),
  ml_methods_roc  = c("RF", "GBM", "ANN"),
  alpha_grid      = seq(0, 1, by = 0.05),
  n_sim_roc       = 1000,
  
  ## Plot colors
  test_colors = c(
    SW = "#CC79A7",
    AD = "#FF7F00",
    JB = "#56B4E9"
  ),
  
  ml_colors = c(
    LR           = "black",
    RF           = "red",
    GBM          = "blue",
    ANN          = "forestgreen",
    SVM          = "green",
    KNN          = "purple",
    MajorityVote = "brown"
  ),
  
  ## Output
  output_dir = "results"
) {
  
  ## ------------------------------------------------------------
  ## Initial setup
  ## ------------------------------------------------------------
  
  make_output_dir(output_dir)
  vip_rank_by <- match.arg(vip_rank_by)
  
  if (is.null(model_save_path)) {
    model_save_path <- file.path(output_dir, "trained_models.RData")
  }
  
  train_pairs <- make_paired_specs(train_alt_specs)
  eval_specs  <- make_eval_specs(eval_alt_specs)
  
  run_settings <- list(
    sample_sizes      = sample_sizes,
    num_sim           = num_sim,
    n_iter_eval       = n_iter_eval,
    models_to_train   = models_to_train,
    threshold         = threshold,
    majority_threshold = majority_threshold,
    standardize_sample = standardize_sample,
    center_by          = center_by,
    k_folds            = k_folds,
    train_alt_specs    = train_alt_specs,
    eval_alt_specs     = eval_alt_specs,
    roc_alt_specs      = roc_alt_specs,
    feature_set        = feature_set
  )
  
  ## ------------------------------------------------------------
  ## Optional density plots
  ## ------------------------------------------------------------
  
  density_plot_file <- NULL
  
  if (isTRUE(make_density_plots)) {
    
    cat("\n", strrep("=", 60), "\n", sep = "")
    cat("DENSITY PLOTS\n")
    cat(strrep("=", 60), "\n")
    
    density_plot_file <- plot_density_pairs(
      paired_specs = make_paired_specs(eval_alt_specs),
      output_dir   = output_dir,
      file_name    = "density_pairs.pdf"
    )
  }
  
  ## ------------------------------------------------------------
  ## Step 1: Train or load models
  ## ------------------------------------------------------------
  
  cat("\n", strrep("=", 60), "\n", sep = "")
  cat("STEP 1: MODEL TRAINING\n")
  cat(strrep("=", 60), "\n")
  
  if (isTRUE(train_models)) {
    
    trained_models <- train_models_all_n(
      sample_sizes       = sample_sizes,
      paired_specs       = train_pairs,
      num_sim            = num_sim,
      models_to_train    = models_to_train,
      feature_set        = feature_set,
      standardize_sample = standardize_sample,
      center_by          = center_by,
      k_folds            = k_folds
    )
    
    save(
      trained_models,
      run_settings,
      file = model_save_path
    )
    
    cat("\nSaved trained models to:", model_save_path, "\n")
    
  } else {
    
    if (!file.exists(model_save_path)) {
      stop(
        "Model file not found: ",
        model_save_path,
        "\nSet train_models = TRUE to train models first."
      )
    }
    
    load(model_save_path)
    
    if (!exists("trained_models")) {
      stop("Loaded file does not contain an object named 'trained_models'.")
    }
    
    cat("Loaded trained models from:", model_save_path, "\n")
  }
  
  ## ------------------------------------------------------------
  ## Step 2: Held-out evaluation
  ## ------------------------------------------------------------
  
  cat("\n", strrep("=", 60), "\n", sep = "")
  cat("STEP 2: HELD-OUT EVALUATION\n")
  cat(strrep("=", 60), "\n")
  
  eval_results_by_n <- evaluate_models_all_n(
    trained_models         = trained_models,
    eval_specs             = eval_specs,
    n_iter                 = n_iter_eval,
    threshold              = threshold,
    majority_threshold     = majority_threshold,
    standardize_sample     = standardize_sample,
    center_by              = center_by,
    store_holdout_features = FALSE
  )
  
  eval_save_path <- file.path(output_dir, "eval_results.RData")
  
  save(
    eval_results_by_n,
    run_settings,
    file = eval_save_path
  )
  
  cat("\nSaved evaluation results to:", eval_save_path, "\n")
  
  ## ------------------------------------------------------------
  ## Step 3: ML-only ROC plots
  ## ------------------------------------------------------------
  
  cat("\n", strrep("=", 60), "\n", sep = "")
  cat("STEP 3: ML ROC PLOTS\n")
  cat(strrep("=", 60), "\n")
  
  save_ml_roc_plots(
    eval_results_by_n = eval_results_by_n,
    sample_sizes      = sample_sizes,
    output_dir        = output_dir,
    ml_colors         = ml_colors
  )
  
  ## ------------------------------------------------------------
  ## Step 4: Optional permutation VIP
  ## ------------------------------------------------------------
  
  vip_results              <- NULL
  vip_summary              <- NULL
  top_features             <- NULL
  selected_feature_metrics <- NULL
  
  if (isTRUE(run_permutation_vip)) {
    
    cat("\n", strrep("=", 60), "\n", sep = "")
    cat("STEP 4: PERMUTATION AUC VARIABLE IMPORTANCE\n")
    cat(strrep("=", 60), "\n")
    
    ## VIP needs stored holdout features. To keep the main evaluation
    ## object lightweight, we rerun evaluation here only when VIP is requested.
    eval_results_for_vip <- evaluate_models_all_n(
      trained_models         = trained_models,
      eval_specs             = eval_specs,
      n_iter                 = n_iter_eval,
      threshold              = threshold,
      majority_threshold     = majority_threshold,
      standardize_sample     = standardize_sample,
      center_by              = center_by,
      store_holdout_features = TRUE
    )
    
    vip_results <- compute_permutation_vip(
      trained_models     = trained_models,
      eval_results_by_n  = eval_results_for_vip,
      n_permutations     = n_permutations,
      models_for_vip     = models_for_vip
    )
    
    vip_summary <- summarize_permutation_vip(vip_results)
    
    top_features <- select_top_vip_features(
      vip_summary           = vip_summary,
      top_n                 = vip_top_n,
      rank_by               = vip_rank_by,
      min_pct_positive      = vip_min_pct_positive,
      require_positive_mean = vip_require_positive_mean
    )
    
    selected_feature_metrics <- get_selected_feature_metrics(
      vip_summary       = vip_summary,
      selected_features = top_features
    )
    
    cat("\nSelected features from permutation VIP:\n")
    print(selected_feature_metrics)
    
    selected_metrics_path <- file.path(
      output_dir,
      "selected_feature_metrics.csv"
    )
    
    write.csv(
      selected_feature_metrics,
      file      = selected_metrics_path,
      row.names = FALSE
    )
    
    cat("\nSaved selected feature metrics to:", selected_metrics_path, "\n")
    
    vip_save_path <- file.path(output_dir, "permutation_vip_results.RData")
    
    save(
      vip_results,
      vip_summary,
      top_features,
      selected_feature_metrics,
      run_settings,
      file = vip_save_path
    )
    
    cat("\nSaved permutation VIP results to:", vip_save_path, "\n")
    
    save_permutation_vip_plots(
      vip_results  = vip_results,
      sample_sizes = sample_sizes,
      models       = models_for_vip,
      top_n        = vip_top_n,
      output_dir   = output_dir
    )
    
  } else {
    
    cat("\nSkipping permutation VIP because run_permutation_vip = FALSE.\n")
  }
  
  ## ------------------------------------------------------------
  ## Step 5: Optional classical-vs-ML ROC comparison
  ## ------------------------------------------------------------
  
  roc_results <- NULL
  
  if (isTRUE(run_roc_analysis)) {
    
    cat("\n", strrep("=", 60), "\n", sep = "")
    cat("STEP 5: CLASSICAL VS ML ROC COMPARISON\n")
    cat(strrep("=", 60), "\n")
    
    roc_results <- run_roc_comparison(
      trained_models  = trained_models,
      roc_alt_specs   = roc_alt_specs,
      sample_sizes    = sample_sizes,
      classical_tests = classical_tests,
      ml_methods      = ml_methods_roc,
      alpha_grid      = alpha_grid,
      n_sim           = n_sim_roc,
      output_dir      = output_dir,
      test_colors     = test_colors,
      ml_colors       = ml_colors
    )
    
    roc_save_path <- file.path(output_dir, "roc_comparison_results.RData")
    
    save(
      roc_results,
      run_settings,
      classical_tests,
      ml_methods_roc,
      file = roc_save_path
    )
    
    cat("\nSaved ROC comparison results to:", roc_save_path, "\n")
    
  } else {
    
    cat("\nSkipping classical-vs-ML ROC comparison because run_roc_analysis = FALSE.\n")
  }
  
  ## ------------------------------------------------------------
  ## Completion
  ## ------------------------------------------------------------
  
  cat("\n", strrep("=", 60), "\n", sep = "")
  cat("PIPELINE COMPLETE\n")
  cat(strrep("=", 60), "\n")
  cat("Output directory:", output_dir, "\n")
  
  invisible(
    list(
      trained_models           = trained_models,
      eval_results_by_n        = eval_results_by_n,
      vip_results              = vip_results,
      vip_summary              = vip_summary,
      top_features             = top_features,
      selected_feature_metrics = selected_feature_metrics,
      roc_results              = roc_results,
      density_plot_file        = density_plot_file,
      run_settings             = run_settings
    )
  )
}



vip_test <- run_normality_framework(
  sample_sizes        = c(10, 50),
  num_sim             = 100,
  n_iter_eval         = 100,
  models_to_train     = c("LR", "RF", "ANN", "GBM", "SVM", "KNN"),
  classical_tests = c("SW", "AD", "JB"),
  ml_methods_roc  = c("RF", "GBM", "ANN"),
  run_permutation_vip = TRUE,
  n_permutations      = 10,
  models_for_vip      = c("RF"),
  run_roc_analysis    = TRUE,
  make_density_plots  = TRUE,
  k_folds             = 10,
  output_dir          = "quick_vip_test_results"
)
