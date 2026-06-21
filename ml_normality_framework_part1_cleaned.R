#' ============================================================
#' ML-BASED NORMALITY TESTING FRAMEWORK
#' Part 1: Setup, Helpers, Distribution Matching, Data Generation
#' ============================================================

## Optional: set this path manually before sourcing this file, if needed.
## Do not hard-code setwd() in a reusable script.
## Example:
## project_dir <- "~/Desktop/OSU/Research/Pretest-Simulation/User_framework_Rpkg/ML_application"
## if (dir.exists(project_dir)) setwd(project_dir)

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
                          center_by   = c("mean", "median")) {
  
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
                                  center_by = "mean") {
  
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
  if (n < 12L) {
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
                               center_by = "mean") {
  
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
                               center_by = "mean",
                               k_folds = 10) {
  
  valid_models <- c("LR", "RF", "ANN", "GBM", "SVM", "KNN")
  bad_models   <- setdiff(models_to_train, valid_models)
  
  if (length(bad_models) > 0L) {
    stop("Unsupported model(s): ", paste(bad_models, collapse = ", "))
  }
  
  center_by <- match.arg(center_by, choices = c("mean", "median"))
  
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
                               center_by = "mean",
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


trained_test <- train_models_all_n(
  sample_sizes    = c(20, 30),
  paired_specs    = paired_specs,
  num_sim         = 10,
  models_to_train = c("LR", "RF", "KNN"),
  k_folds         = 3
)

names(trained_test)

lapply(trained_test, function(obj) names(obj$models))

lapply(trained_test, function(obj) {
  lapply(obj$models, function(m) m$bestTune)
})