#' ============================================================
#' ML-BASED NORMALITY TESTING FRAMEWORK
#' Part 1: Setup, Helpers, Distribution Matching, Data Generation
#' ============================================================

setwd("~/Desktop/OSU/Research/Pretest-Simulation/User_framework_Rpkg/one_sample")

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
        iqr    = qexp(0.75, rate = rate) -
          qexp(0.25, rate = rate)
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
          sqrt(
            2 * df2^2 * (df1 + df2 - 2) / (df1 * (df2 - 2)^2 * (df2 - 4))
          )
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
        iqr    = xm * (0.25)^(-1 / alpha) -
          xm * (0.75)^(-1 / alpha)
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
        iqr    = qnorm(0.75, mu, sd_mix) -
          qnorm(0.25, mu, sd_mix)
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
# The goal is to give the ML models many complementary
#' signals about shape, tails, spacing, entropy, and Q-Q fit.
#' ============================================================

#' ============================================================
#' SECTION 6: BASIC FEATURE HELPERS
#' ============================================================

#' Standardize a numeric sample using its sample mean and sample SD.
standardize_sample <- function(x) {
  
  x <- as.numeric(na.omit(x))
  s <- sd(x)
  
  if (!is.finite(s) || s <= 0) {
    return(rep(0, length(x)))
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
  
  s <- sd(x)
  
  data.frame(
    SW_stat  = safe_calc(as.numeric(shapiro.test(x)$statistic)),
    SF_stat  = safe_calc(as.numeric(nortest::sf.test(x)$statistic)),
    AD_stat  = safe_calc(as.numeric(nortest::ad.test(x)$statistic)),
    LF_stat  = safe_calc(as.numeric(nortest::lillie.test(x)$statistic)),
    KS_stat  = safe_calc(as.numeric(ks.test(x, "pnorm", mean(x), s)$statistic)),
    JB_stat  = safe_calc(as.numeric(tseries::jarque.bera.test(x)$statistic))
  )
}


#' ============================================================
#' SECTION 8: ENTROPY AND SPACING FEATURES
#' ============================================================

#' Lin-Mudholkar Zp statistic.
#'
#' The statistic is based on the correlation between the ordered
#' sample and leave-one-out variance terms. It is designed to
#' capture departures from Normality through the relationship
#' between order statistics and local dispersion.
calculate_zp_statistic <- function(x) {
  
  x <- sort(as.numeric(na.omit(x)))
  n <- length(x)
  
  if (n < 4L) {
    return(NA_real_)
  }
  
  n1 <- n - 1L
  
  # Leave-one-out variance for each observation.
  loo_var <- (sum(x^2) - x^2) / n1 -
    ((sum(x) - x) / n1)^2
  
  # The statistic uses the cube root of the leave-one-out variances.
  safe_calc(atanh(cor(x, loo_var^(1 / 3))))
}


#' Vasicek entropy-based statistic.
#'
#' This statistic estimates entropy using order-statistic spacings.
#' For smooth unimodal distributions, the spacing structure is
#' informative about shape. Normal samples tend to produce a stable
#' spacing pattern after scaling by the sample SD.
calculate_vasicek_kmn <- function(x) {
  
  x <- sort(as.numeric(na.omit(x)))
  n <- length(x)
  
  if (n < 4L) {
    return(NA_real_)
  }
  
  # A common window choice for Vasicek-type entropy estimators.
  m <- floor(sqrt(n))
  
  s <- sd(x)
  
  if (!is.finite(s) || s <= 0) {
    return(NA_real_)
  }
  
  # Windowed spacings x_(i+m) - x_(i-m), with boundary correction.
  upper <- pmin(n, seq_len(n) + m)
  lower <- pmax(1L, seq_len(n) - m)
  
  spacings <- x[upper] - x[lower]
  
  if (any(spacings <= 0, na.rm = TRUE)) {
    return(NA_real_)
  }
  
  (n / (2 * m * s)) * exp(mean(log(spacings)))
}



#' Rényi entropy using a histogram approximation.
#'
#' This summarizes how concentrated or dispersed the empirical
#' distribution is across bins. It is included as a broad shape
#' feature rather than as a formal Normality test.
renyi_entropy <- function(x, alpha = 2, bins = 10) {
  
  x <- as.numeric(na.omit(x))
  
  if (length(x) < 3L) {
    return(NA_real_)
  }
  
  # Histogram probabilities.
  p <- table(cut(x, breaks = bins, include.lowest = TRUE))
  p <- as.numeric(p) / sum(p)
  
  # Empty bins do not contribute to entropy.
  p <- p[p > 0]
  
  (1 / (1 - alpha)) * log(sum(p^alpha))
}



#' Shannon entropy using a histogram approximation.
#'
#' This is the usual empirical entropy, -sum p log(p), computed
#' from binned data. It captures concentration and spread of the
#' sample distribution.
shannon_entropy_stat <- function(x, bins = 10) {
  
  x <- as.numeric(na.omit(x))
  
  if (length(x) < 3L) {
    return(NA_real_)
  }
  
  p <- table(cut(x, breaks = bins, include.lowest = TRUE))
  p <- as.numeric(p) / sum(p)
  p <- p[p > 0]
  
  -sum(p * log(p))
}


#' Sample entropy.
#'
#' Sample entropy measures the regularity of a sequence. Although
#' it is more common in time-series settings, it can still act as
#' a nonlinear summary of repeated local patterns in the sample.
#' Here it is used only as a descriptive ML feature.
sample_entropy_stat <- function(x) {
  
  x <- as.numeric(na.omit(x))
  
  if (length(x) < 5L) {
    return(NA_real_)
  }
  
  s <- sd(x)
  
  if (!is.finite(s) || s <= 0) {
    return(NA_real_)
  }
  
  # Tolerance is set relative to the sample SD.
  safe_calc(pracma::sample_entropy(x, edim = 2, r = 0.2 * s))
}

#' Ratio of two Vasicek-type entropy estimates.
#'
#' The same sample is evaluated using two spacing-window widths.
#' A ratio far from one suggests that the estimated entropy is
#' unstable across scales, which can occur under non-Normal shapes.
entropy_ratio <- function(x, m1 = 2, m2 = 4) {
  
  x <- sort(as.numeric(na.omit(x)))
  n <- length(x)
  
  if (n <= m2 + 1L) {
    return(NA_real_)
  }
  
  vasicek_entropy <- function(z, m) {
    
    n <- length(z)
    
    # Positive m-step spacings are required before taking logs.
    d <- z[(m + 1):n] - z[1:(n - m)]
    d <- d[d > 0]
    
    if (length(d) == 0L) {
      return(NA_real_)
    }
    
    mean(log((n / m) * d))
  }
  
  v1 <- vasicek_entropy(x, m1)
  v2 <- vasicek_entropy(x, m2)
  
  safe_ratio(v1, v2)
}



#' Approximate negentropy.
#'
#' Negentropy measures departure from Gaussianity. The exact
#' quantity is difficult to compute, so this uses the common
#' approximation based on G(u) = log(cosh(u)).
#'
#' Larger values suggest stronger non-Normality.
negentropy_approx <- function(x) {
  
  z <- standardize_sample(x)
  
  if (length(z) < 3L) {
    return(NA_real_)
  }
  
  g <- function(u) {
    log(cosh(u))
  }
  
  # A fixed reference sample approximates E[G(Z)] for Z ~ N(0, 1).
  z0 <- rnorm(10000)
  
  safe_calc((mean(g(z)) - mean(g(z0)))^2)
}


#' Variance and skewness of log spacings.
#'
#' Consecutive order-statistic spacings describe how observations
#' are distributed across the support. The variance and skewness
#' of log spacings can detect clustering, gaps, and tail behavior.
calc_spacing_ratio_stats <- function(x) {
  
  x <- sort(as.numeric(na.omit(x)))
  
  # Consecutive spacings between ordered observations.
  d <- diff(x)
  d <- d[d > 0]
  
  if (length(d) < 2L) {
    return(list(
      var_log_spacings  = NA_real_,
      skew_log_spacings = NA_real_
    ))
  }
  
  log_d <- log(d)
  
  list(
    var_log_spacings  = var(log_d),
    skew_log_spacings = safe_calc(as.numeric(moments::skewness(log_d)))
  )
}

#' Greenwood spacing statistic.
#'
#' The Greenwood statistic is based on squared normalized spacings.
#' Large values occur when spacings are uneven, which may indicate
#' clustering, gaps, or heavy-tailed behavior.
greenwood_stat <- function(x) {
  
  x <- sort(as.numeric(na.omit(x)))
  n <- length(x)
  
  if (n < 4L) {
    return(NA_real_)
  }
  
  d <- diff(x)
  total_d <- sum(d)
  
  if (!is.finite(total_d) || total_d <= 0) {
    return(NA_real_)
  }
  
  # Greenwood statistic based on normalized spacings.
  g <- sum(d^2) / total_d^2
  
  # Approximate mean and variance used for standardization.
  expected_g <- 2 / n
  var_g      <- 4 * (n - 1) / (n^2 * (n + 1))
  
  safe_ratio(g - expected_g, sqrt(var_g))
}


#' Rao spacing statistic.
#'
#' This statistic measures the total absolute deviation of
#' normalized spacings from equal spacing. It is useful as a
#' general spacing irregularity feature.
rao_spacing_stat <- function(x) {
  
  x <- sort(as.numeric(na.omit(x)))
  n <- length(x)
  
  if (n < 4L) {
    return(NA_real_)
  }
  
  d   <- diff(x)
  rng <- x[n] - x[1L]
  
  if (!is.finite(rng) || rng <= 0) {
    return(NA_real_)
  }
  
  m <- n - 1L
  
  # Normalize spacings so that they sum to one.
  s <- d / rng
  
  # Total absolute deviation from equal spacing.
  cn <- (n / 2) * sum(abs(s - 1 / m))
  
  # Standardized version used as an ML feature.
  cn_scaled <- cn / (n / 2)
  mean0     <- 2 * (m - 1) / m
  sd0       <- sqrt(2 / (n * pi))
  
  safe_ratio(cn_scaled - mean0, sd0)
}

#' ============================================================
#' SECTION 9: TAIL, OUTLIER, AND ROBUSTNESS FEATURES
#' ============================================================

#' Tail asymmetry ratio.
#'
#' This compares the average upper-tail excess to the average
#' lower-tail deficit, relative to the sample mean.
#'
#' For symmetric distributions, this ratio should be close to 1.
#' Values far from 1 indicate asymmetric tail behavior.
calculate_tail_asymmetry <- function(x) {
  
  x <- as.numeric(na.omit(x))
  
  if (length(x) < 6L) {
    return(NA_real_)
  }
  
  mu <- mean(x)
  
  q10 <- quantile(x, 0.10, names = FALSE)
  q90 <- quantile(x, 0.90, names = FALSE)
  
  upper_tail <- x[x > q90]
  lower_tail <- x[x < q10]
  
  if (length(upper_tail) == 0L || length(lower_tail) == 0L) {
    return(NA_real_)
  }
  
  upper_excess <- mean(upper_tail) - mu
  lower_deficit <- mu - mean(lower_tail)
  
  safe_ratio(upper_excess, lower_deficit)
}

#' Hill tail-index ratio.
#'
#' This estimates tail heaviness separately in the upper and lower
#' tails using a Hill-type estimator, then returns their ratio.
#'
#' A value near 1 suggests roughly balanced tails. Values far from
#' 1 suggest tail asymmetry.
tail_index_ratio <- function(x, k = NULL) {
  
  x <- as.numeric(na.omit(x))
  n <- length(x)
  
  if (n < 6L) {
    return(NA_real_)
  }
  
  if (is.null(k)) {
    k <- max(2L, floor(sqrt(n)))
  }
  
  hill_one_side <- function(z, k_side) {
    
    z <- sort(z, decreasing = TRUE)
    n_side <- length(z)
    
    k_side <- min(k_side, n_side - 1L)
    
    if (k_side < 1L) {
      return(NA_real_)
    }
    
    # Shift values if needed so logs are well-defined.
    if (min(z, na.rm = TRUE) <= 0) {
      z <- z - min(z, na.rm = TRUE) + 1
    }
    
    threshold <- z[k_side + 1L]
    
    if (!is.finite(threshold) || threshold <= 0) {
      return(NA_real_)
    }
    
    hill_val <- mean(log(z[seq_len(k_side)]) - log(threshold))
    
    if (!is.finite(hill_val) || hill_val <= 0) {
      return(NA_real_)
    }
    
    1 / hill_val
  }
  
  upper_index <- safe_calc(hill_one_side(x,  k))
  lower_index <- safe_calc(hill_one_side(-x, k))
  
  safe_ratio(upper_index, lower_index)
}

#' Pickands tail-index estimator.
#'
#' Pickands' estimator is based on high order statistics. It is
#' useful for detecting tail behavior, especially heavy-tail or
#' bounded-tail departures from Normality.
#'
#' Positive values suggest heavier tails; negative values suggest
#' lighter or bounded tails.
pickands_tail_index <- function(x, k = NULL) {
  
  x <- sort(as.numeric(na.omit(x)), decreasing = TRUE)
  n <- length(x)
  
  if (n < 6L) {
    return(NA_real_)
  }
  
  if (is.null(k)) {
    k <- max(1L, floor(n / 5L))
  }
  
  # Pickands requires indices k, 2k, and 4k to exist.
  k <- min(k, floor(n / 4L))
  
  if (k < 1L || 4L * k > n) {
    return(NA_real_)
  }
  
  numerator   <- x[k]      - x[2L * k]
  denominator <- x[2L * k] - x[4L * k]
  
  if (!is.finite(numerator) || !is.finite(denominator)) {
    return(NA_real_)
  }
  
  if (numerator <= 0 || denominator <= 0) {
    return(NA_real_)
  }
  
  log(numerator / denominator) / log(2)
}

#' Mean excess slope.
#'
#' The mean excess function is
#'   e(u) = E[X - u | X > u].
#'
#' Heavy-tailed distributions tend to show increasing mean excess
#' as the threshold increases. This function estimates the slope
#' of the empirical mean excess curve over upper-tail thresholds.
mean_excess_slope <- function(x) {
  
  x <- as.numeric(na.omit(x))
  n <- length(x)
  
  if (n < 8L) {
    return(NA_real_)
  }
  
  # Use upper-half thresholds while keeping enough exceedances.
  probs <- seq(
    from       = 0.50,
    to         = 1 - 2 / n,
    length.out = min(8L, floor(n / 3L))
  )
  
  thresholds <- quantile(x, probs = probs, names = FALSE)
  
  mean_excess <- vapply(thresholds, function(u) {
    
    excess <- x[x > u] - u
    
    if (length(excess) < 2L) {
      return(NA_real_)
    }
    
    mean(excess)
  }, numeric(1L))
  
  ok <- is.finite(thresholds) & is.finite(mean_excess)
  
  if (sum(ok) < 3L) {
    return(NA_real_)
  }
  
  safe_calc(as.numeric(coef(lm(mean_excess[ok] ~ thresholds[ok]))[2L]))
}

#' Normalized tail-weight ratio.
#'
#' This compares the observed 5%-to-95% spread of the standardized
#' sample to the corresponding standard Normal spread.
#'
#' Values greater than 1 suggest heavier or more dispersed tails
#' than a Normal distribution.
calculate_twr_normalized <- function(x) {
  
  z <- standardize_sample(x)
  
  if (length(z) < 4L) {
    return(NA_real_)
  }
  
  observed_spread <- quantile(z, 0.95, names = FALSE) -
    quantile(z, 0.05, names = FALSE)
  
  normal_spread <- qnorm(0.95) - qnorm(0.05)
  
  safe_ratio(observed_spread, normal_spread)
}

#' Excess tail proportion.
#'
#' This computes the proportion of standardized observations beyond
#' +/- 2 SD and scales it by the standard Normal reference value.
#'
#' The Normal reference is approximately P(|Z| > 2) = 0.0455.
calculate_excess_tail_prop <- function(x) {
  
  z <- standardize_sample(x)
  
  if (length(z) < 4L) {
    return(NA_real_)
  }
  
  mean(abs(z) > 2) / 0.0455
}

#' Tukey outlier proportion.
#'
#' This is the proportion of observations outside the usual Tukey
#' boxplot fences:
#'   Q1 - 1.5 IQR  and  Q3 + 1.5 IQR.
#'
#' Heavy-tailed distributions tend to produce larger values.
calculate_outlier_proportion <- function(x) {
  
  x <- as.numeric(na.omit(x))
  
  if (length(x) < 4L) {
    return(NA_real_)
  }
  
  q <- quantile(x, probs = c(0.25, 0.75), names = FALSE)
  h <- 1.5 * IQR(x)
  
  mean(x < q[1] - h | x > q[2] + h)
}

#' Qn robust spread ratio.
#'
#' Qn is a robust scale estimator. This feature compares Qn to
#' the usual sample SD.
#'
#' Under Normality, the ratio should be relatively stable. Heavy
#' tails and outliers often affect the SD more strongly than Qn.
calc_qn_robust_spread <- function(x) {
  
  x <- as.numeric(na.omit(x))
  
  if (length(x) < 5L) {
    return(NA_real_)
  }
  
  s <- sd(x)
  
  if (!is.finite(s) || s <= 0) {
    return(NA_real_)
  }
  
  qn_val <- safe_calc(
    as.numeric(robustbase::Qn(x, finite.corr = TRUE))
  )
  
  safe_ratio(qn_val, s)
}

#' Medcouple robust skewness.
#'
#' The medcouple is a robust measure of skewness. It is less
#' sensitive to extreme observations than ordinary moment skewness.
#'
#' Values near 0 suggest symmetry. Positive values suggest right
#' skewness; negative values suggest left skewness.
calc_medcouple <- function(x) {
  
  x <- as.numeric(na.omit(x))
  
  if (length(x) < 5L) {
    return(NA_real_)
  }
  
  safe_calc(as.numeric(robustbase::mc(x)))
}

#' Left and right tail weights from medcouple ideas.
#'
#' This computes robust tail-weight summaries separately for the
#' left and right tails. The tails are measured relative to the
#' sample median.
#'
#' LTW:
#'   robust left-tail weight
#'
#' RTW:
#'   robust right-tail weight
#'
#' These are useful because two samples can have similar skewness
#' but different tail behavior.
calc_tail_weights_brys <- function(x) {
  
  x <- as.numeric(na.omit(x))
  
  if (length(x) < 6L) {
    return(list(
      LTW = NA_real_,
      RTW = NA_real_
    ))
  }
  
  m <- median(x)
  
  # Left-tail distances from the median.
  left_weight <- safe_calc({
    
    left_dist <- m - x[x <= m]
    left_dist <- left_dist[left_dist > 0]
    
    if (length(left_dist) < 3L) {
      return(NA_real_)
    }
    
    -as.numeric(robustbase::mc(left_dist))
  })
  
  # Right-tail distances from the median.
  right_weight <- safe_calc({
    
    right_dist <- x[x >= m] - m
    right_dist <- right_dist[right_dist > 0]
    
    if (length(right_dist) < 3L) {
      return(NA_real_)
    }
    
    as.numeric(robustbase::mc(right_dist))
  })
  
  list(
    LTW = left_weight,
    RTW = right_weight
  )
}


#' ============================================================
#' SECTION 10: Q-Q AND CORRELATION FEATURES
#' ============================================================
#' 

#' Cubic coefficient from a Normal Q-Q polynomial fit.
#'
#' A Normal Q-Q plot should be approximately linear. This function
#' fits a cubic polynomial of empirical quantiles on theoretical
#' Normal quantiles and returns the cubic coefficient.
#'
#' A large cubic term captures systematic S-shaped or reverse
#' S-shaped curvature in the Q-Q plot.
qq_cubic_coef <- function(x) {
  
  x <- as.numeric(na.omit(x))
  n <- length(x)
  
  if (n < 6L) {
    return(NA_real_)
  }
  
  # Theoretical Normal quantiles.
  theo_q <- qnorm(ppoints(n))
  
  # Empirical order statistics.
  emp_q <- sort(x)
  
  fit <- suppressWarnings(
    lm(emp_q ~ poly(theo_q, 3, raw = TRUE))
  )
  
  coefs <- coef(fit)
  
  # Intercept, linear, quadratic, cubic.
  if (length(coefs) < 4L || !is.finite(coefs[4L])) {
    return(NA_real_)
  }
  
  as.numeric(coefs[4L])
}

#' Standard deviation of Q-Q residuals.
#'
#' A Normal Q-Q plot should be close to a straight line. This
#' feature fits a linear model between empirical quantiles and
#' theoretical Normal quantiles, then returns the residual SD.
#'
#' Larger values suggest stronger Q-Q plot departures from
#' linearity.
calculate_qq_resid_sd <- function(x) {
  
  x <- as.numeric(na.omit(x))
  n <- length(x)
  
  if (n < 4L) {
    return(NA_real_)
  }
  
  theo_q <- qnorm(ppoints(n))
  emp_q  <- sort(x)
  
  fit <- lm(emp_q ~ theo_q)
  
  sd(residuals(fit))
}

#' Pearson's second skewness coefficient.
#'
#' This is defined as:
#'   3 * (mean - median) / SD.
#'
#' It is a simple location-based skewness measure. For symmetric
#' distributions, the mean and median should be close, so the
#' value should be near 0.
pearson_cor_skew <- function(x) {
  
  x <- as.numeric(na.omit(x))
  
  if (length(x) < 3L) {
    return(NA_real_)
  }
  
  safe_ratio(
    num = 3 * (mean(x) - median(x)),
    den = sd(x)
  )
}

#' de Wet-Venter weighted correlation statistic.
#'
#' This statistic is a weighted squared correlation between the
#' ordered sample and expected Normal scores. The weights are
#' based on the Normal density at the plotting positions.
#'
#' Values close to 1 indicate that the ordered sample aligns well
#' with expected Normal order-statistic behavior.
dewet_venter_stat <- function(x) {
  
  x <- sort(as.numeric(na.omit(x)))
  n <- length(x)
  
  if (n < 5L) {
    return(NA_real_)
  }
  
  # Blom-type plotting positions.
  p <- (seq_len(n) - 3 / 8) / (n + 1 / 4)
  
  # Expected Normal scores.
  m <- qnorm(p)
  
  # Density weights: central Normal scores receive more weight.
  w <- dnorm(m)
  
  w_sum <- sum(w)
  
  if (!is.finite(w_sum) || w_sum <= 0) {
    return(NA_real_)
  }
  
  # Weighted means.
  x_bar <- sum(w * x) / w_sum
  m_bar <- sum(w * m) / w_sum
  
  # Weighted centered values.
  x_dev <- x - x_bar
  m_dev <- m - m_bar
  
  numerator <- sum(w * x_dev * m_dev)
  denom_x   <- sum(w * x_dev^2)
  denom_m   <- sum(w * m_dev^2)
  
  r <- safe_ratio(
    num = numerator,
    den = sqrt(denom_x * denom_m)
  )
  
  if (!is.finite(r)) {
    return(NA_real_)
  }
  
  # Bound numerically before squaring.
  r <- max(-1, min(1, r))
  
  r^2
}


#' Ryan-Joiner statistic.
#'
#' The Ryan-Joiner statistic is the squared correlation between
#' ordered observations and expected Normal scores.
#'
#' It is closely related to Shapiro-Francia. Values close to 1
#' support Normal-like Q-Q behavior.
ryan_joiner_stat <- function(x) {
  
  x <- sort(as.numeric(na.omit(x)))
  n <- length(x)
  
  if (n < 5L) {
    return(NA_real_)
  }
  
  p <- (seq_len(n) - 3 / 8) / (n + 1 / 4)
  m <- qnorm(p)
  
  r <- safe_calc(cor(x, m))
  
  if (!is.finite(r)) {
    return(NA_real_)
  }
  
  r^2
}

#' Modified Shapiro-Francia statistic.
#'
#' This is a squared correlation between ordered observations
#' and expected Normal scores using Blom plotting positions.
#'
#' It is included as a smooth correlation-based Normality feature.
modified_shapiro_francia <- function(x) {
  
  x <- sort(as.numeric(na.omit(x)))
  n <- length(x)
  
  if (n < 5L) {
    return(NA_real_)
  }
  
  p <- (seq_len(n) - 3 / 8) / (n + 1 / 4)
  expected <- qnorm(p)
  
  r <- safe_calc(cor(x, expected))
  
  if (!is.finite(r)) {
    return(NA_real_)
  }
  
  r^2
}


#' Robust biweight midcorrelation with Normal scores.
#'
#' This is a robust correlation between the ordered sample and
#' expected Normal scores. The biweight weights reduce the effect
#' of extreme points.
#'
#' This can be useful because ordinary Q-Q correlations may be
#' overly influenced by a few extreme observations.
biweight_midcorrelation_norm <- function(x) {
  
  x <- sort(as.numeric(na.omit(x)))
  n <- length(x)
  
  if (n < 8L) {
    return(NA_real_)
  }
  
  expected <- qnorm(ppoints(n))
  
  med_x <- median(x)
  med_e <- median(expected)
  
  mad_x <- mad(x)
  mad_e <- mad(expected)
  
  if (!is.finite(mad_x) || !is.finite(mad_e)) {
    return(NA_real_)
  }
  
  if (mad_x <= 0 || mad_e <= 0) {
    return(NA_real_)
  }
  
  # Standard biweight tuning constant.
  c_val <- 9
  
  u_x <- (x - med_x) / (c_val * mad_x)
  u_e <- (expected - med_e) / (c_val * mad_e)
  
  # Tukey biweight weights.
  w_x <- ifelse(abs(u_x) < 1, (1 - u_x^2)^2, 0)
  w_e <- ifelse(abs(u_e) < 1, (1 - u_e^2)^2, 0)
  
  if (sum(w_x) <= 0 || sum(w_e) <= 0) {
    return(NA_real_)
  }
  
  x_centered <- (x - med_x) * w_x
  e_centered <- (expected - med_e) * w_e
  
  numerator <- sum(x_centered * e_centered)
  denominator <- sqrt(sum(x_centered^2) * sum(e_centered^2))
  
  safe_ratio(numerator, denominator)
}

#' ============================================================
#' SECTION 11: ECF AND DISTANCE FEATURES
#' ============================================================
#' 
#' Epps-Pulley statistic.
#'
#' The Epps-Pulley statistic compares the empirical characteristic
#' function of the standardized sample to the characteristic
#' function of a standard Normal distribution.
#'
#' It is sensitive to broad departures from Normality, including
#' skewness, kurtosis, and other shape differences.
epps_pulley_stat <- function(x) {
  
  x <- as.numeric(na.omit(x))
  n <- length(x)
  
  if (n < 3L) {
    return(NA_real_)
  }
  
  z <- standardize_sample(x)
  
  # Pairwise squared distances among standardized observations.
  d <- outer(z, z, function(a, b) (a - b)^2)
  
  # Terms from the Epps-Pulley integrated squared distance formula.
  term1 <- (1 / n) * sum(exp(-d / 4))
  term2 <- (2 / n) * sqrt(2) * sum(exp(-z^2 / 4))
  term3 <- n / sqrt(3)
  
  safe_calc(term1 - term2 + term3)
}

#' Henze-Zirkler statistic.
#'
#' This statistic is also based on the empirical characteristic
#' function, using a weighted L2 distance between the sample ECF
#' and the Normal characteristic function.
#'
#' Larger values generally indicate stronger departure from
#' Normality.
henze_zirkler_stat <- function(x) {
  
  x <- as.numeric(na.omit(x))
  n <- length(x)
  
  if (n < 5L) {
    return(NA_real_)
  }
  
  z <- standardize_sample(x)
  
  # Smoothing parameter commonly used in the Henze-Zirkler statistic.
  beta <- (1 / sqrt(2)) * ((2 * n + 1) / 4)^(1 / 7)
  b    <- beta^2
  
  # Pairwise squared distances among standardized observations.
  d <- outer(z, z, function(a, bb) (a - bb)^2)
  
  term1 <- (1 / n) * sum(exp(-b * d / 2))
  term2 <- 2 * (1 + b)^(-0.5) * sum(exp(-b * z^2 / (2 * (1 + b))))
  term3 <- n * (1 + 2 * b)^(-0.5)
  
  safe_calc(term1 - term2 + term3)
}


#' Empirical characteristic function deviations.
#'
#' For a standardized sample, this computes the difference between
#' the empirical characteristic function and the standard Normal
#' characteristic function at fixed t values.
#'
#' The real part compares mean cos(tX) with exp(-t^2 / 2).
#' The imaginary part uses mean sin(tX), which is zero under a
#' symmetric standard Normal distribution.
ecf_deviations <- function(x, t_vals = c(0.5, 1.0, 2.0)) {
  
  x <- as.numeric(na.omit(x))
  
  out_names <- c(
    paste0("ecf_re_t", seq_along(t_vals)),
    paste0("ecf_im_t", seq_along(t_vals))
  )
  
  if (length(x) < 3L) {
    return(setNames(rep(NA_real_, length(out_names)), out_names))
  }
  
  z <- standardize_sample(x)
  
  # Real deviations from the standard Normal characteristic function.
  real_dev <- vapply(
    t_vals,
    function(t) {
      mean(cos(t * z)) - exp(-t^2 / 2)
    },
    numeric(1L)
  )
  
  # Imaginary deviations. These are zero for symmetric Normal data.
  imag_dev <- vapply(
    t_vals,
    function(t) {
      mean(sin(t * z))
    },
    numeric(1L)
  )
  
  setNames(c(real_dev, imag_dev), out_names)
}

#' Empirical moment-ratio distance.
#'
#' This combines skewness and kurtosis deviations from their
#' Normal-reference values into one standardized distance.
#'
#' Under Normality:
#'   skewness = 0
#'   kurtosis = 3
#'
#' The statistic uses approximate variances 6/n and 24/n for
#' sample skewness and kurtosis.
empirical_moment_ratios <- function(x) {
  
  x <- as.numeric(na.omit(x))
  n <- length(x)
  
  if (n < 10L) {
    return(NA_real_)
  }
  
  centered_x <- x - mean(x)
  m2 <- mean(centered_x^2)
  
  if (!is.finite(m2) || m2 <= 0) {
    return(NA_real_)
  }
  
  # Moment skewness and kurtosis.
  skew_val <- mean(centered_x^3) / m2^(3 / 2)
  kurt_val <- mean(centered_x^4) / m2^2
  
  # Distance from the Normal-reference vector (0, 3).
  diff_vec <- c(skew_val, kurt_val - 3)
  
  # Approximate covariance matrix under Normality.
  cov_mat <- diag(c(6 / n, 24 / n))
  
  safe_calc(
    as.numeric(sqrt(t(diff_vec) %*% solve(cov_mat) %*% diff_vec))
  )
}


#' Energy distance to a fitted Normal sample.
#'
#' This approximates the energy distance between the empirical
#' sample distribution and a Normal distribution with the same
#' sample mean and sample SD.
#'
#' Energy distance is sensitive to general distributional
#' differences, not only moment differences.
energy_distance_to_normal <- function(x) {
  
  x <- as.numeric(na.omit(x))
  n <- length(x)
  
  if (n < 5L) {
    return(NA_real_)
  }
  
  mu <- mean(x)
  s  <- sd(x)
  
  if (!is.finite(s) || s <= 0) {
    return(NA_real_)
  }
  
  # Simulated fitted Normal comparison sample.
  z <- rnorm(n, mean = mu, sd = s)
  
  # Average pairwise distances.
  d_xx <- mean(abs(outer(x, x, "-")))
  d_zz <- mean(abs(outer(z, z, "-")))
  d_xz <- mean(abs(outer(x, z, "-")))
  
  safe_calc(2 * d_xz - d_xx - d_zz)
}

#' ============================================================
#' SECTION 12: OTHER SHAPE FEATURES
#' ============================================================
#' 
#' Moors kurtosis.
#'
#' Moors kurtosis is based on octiles rather than ordinary fourth
#' moments. This makes it less sensitive to extreme observations
#' than classical kurtosis.
#'
#' Larger values usually indicate heavier tails or more peaked
#' behavior relative to the middle spread of the sample.
calc_moors_kurtosis <- function(x) {
  
  x <- as.numeric(na.omit(x))
  
  if (length(x) < 8L) {
    return(NA_real_)
  }
  
  # Octiles: O1, ..., O7.
  octiles <- quantile(
    x,
    probs = seq(1 / 8, 7 / 8, by = 1 / 8),
    names = FALSE
  )
  
  numerator <- (octiles[7] - octiles[5]) +
    (octiles[3] - octiles[1])
  
  denominator <- octiles[6] - octiles[2]
  
  safe_ratio(numerator, denominator)
}

#' D'Agostino omnibus statistic.
#'
#' This combines skewness and kurtosis information into one
#' omnibus departure-from-Normality statistic.
#'
#' It is included as a shape feature, not as a final hypothesis
#' test decision rule.
dagostino_omnibus <- function(x) {
  
  x <- as.numeric(na.omit(x))
  n <- length(x)
  
  if (n < 20L) {
    return(NA_real_)
  }
  
  centered_x <- x - mean(x)
  m2 <- mean(centered_x^2)
  
  if (!is.finite(m2) || m2 <= 0) {
    return(NA_real_)
  }
  
  g1 <- mean(centered_x^3) / m2^(3 / 2)
  g2 <- mean(centered_x^4) / m2^2
  
  # Skewness transformation.
  y <- g1 * sqrt((n + 1) * (n + 3) / (6 * (n - 2)))
  
  beta2 <- 3 * (n^2 + 27 * n - 70) * (n + 1) * (n + 3) /
    ((n - 2) * (n + 5) * (n + 7) * (n + 9))
  
  w2 <- sqrt(2 * (beta2 - 1))
  
  if (!is.finite(w2) || w2 <= 1) {
    return(NA_real_)
  }
  
  delta <- 1 / sqrt(log(sqrt(w2)))
  alpha <- sqrt(2 / (w2 - 1))
  
  z1 <- delta * log(y / alpha + sqrt((y / alpha)^2 + 1))
  
  # Kurtosis transformation.
  if (!is.finite(g2) || g2 <= 0) {
    return(NA_real_)
  }
  
  a <- 6 + (8 / g2) * (2 / g2 + sqrt(1 + 4 / g2^2))
  
  if (!is.finite(a) || a <= 4) {
    return(NA_real_)
  }
  
  numerator <- 1 - 2 / (9 * a)
  
  denominator <- 1 + sqrt(2 / (a - 4)) *
    (g2 - 3 - 6 / (n + 1)) / sqrt(24 / (n + 1))
  
  if (!is.finite(denominator) || abs(denominator) < .Machine$double.eps) {
    return(NA_real_)
  }
  
  z2 <- (numerator - ((1 - 2 / a) / denominator)^(1 / 3)) /
    sqrt(2 / (9 * a))
  
  safe_calc(as.numeric(z1^2 + z2^2))
}

#' Lorenz asymmetry coefficient.
#'
#' This measures asymmetry using the empirical Lorenz curve.
#' In this ML setting, it is used as a distribution-shape feature.
#'
#' Because the Lorenz curve depends on cumulative sums, this
#' feature is most stable when the sample values are not centered
#' around a total close to zero.
lorenz_asymmetry_coef <- function(x) {
  
  x <- sort(as.numeric(na.omit(x)))
  n <- length(x)
  
  if (n < 10L) {
    return(NA_real_)
  }
  
  total_x <- sum(x)
  
  if (!is.finite(total_x) || abs(total_x) < .Machine$double.eps) {
    return(NA_real_)
  }
  
  p <- seq_len(n) / n
  l <- cumsum(x) / total_x
  
  # Find where p + L(p) is closest to 1.
  idx <- which.min(abs(p + l - 1))
  
  2 * abs(p[idx] - 0.5)
}


#' Density derivative ratio.
#'
#' This uses a kernel density estimate and compares the maximum
#' absolute slope to the maximum density height.
#'
#' Sharply peaked or irregular distributions often have larger
#' derivative-to-height ratios than smoother Normal-like samples.
density_derivative_ratio <- function(x) {
  
  x <- as.numeric(na.omit(x))
  
  if (length(x) < 10L) {
    return(NA_real_)
  }
  
  dens <- safe_calc(
    density(x, n = 512),
    default = NULL
  )
  
  if (is.null(dens)) {
    return(NA_real_)
  }
  
  dx <- diff(dens$x)[1]
  dy <- diff(dens$y)
  
  if (!is.finite(dx) || dx <= 0) {
    return(NA_real_)
  }
  
  max_abs_slope <- max(abs(dy / dx), na.rm = TRUE)
  max_density   <- max(dens$y, na.rm = TRUE)
  
  safe_ratio(max_abs_slope, max_density)
}

#' Kolmogorov-Smirnov distance from fitted standard Normal.
#'
#' The sample is first standardized using its own mean and SD.
#' The statistic then measures the maximum distance between the
#' empirical CDF of the standardized sample and the standard
#' Normal CDF.
ks_distance <- function(x) {
  
  x <- as.numeric(na.omit(x))
  
  if (length(x) < 4L) {
    return(NA_real_)
  }
  
  z <- standardize_sample(x)
  
  safe_calc(as.numeric(ks.test(z, "pnorm")$statistic))
}


#' ============================================================
#' SECTION 13: MASTER FEATURE WRAPPER
#' ============================================================
#' 
#' Compute all features for one sample.
#'
#' This is the main feature-generation function used by the ML
#' pipeline. It takes one numeric sample and returns a one-row
#' data frame containing location, spread, shape, tail, spacing,
#' entropy, Q-Q, ECF, and interaction features.
#'
#' Important:
#'   - Each column is one feature.
#'   - The function keeps all features that were previously
#'     computed but not always returned in the original code.
#'   - NA values are imputed at the end so caret can train models
#'     without failing.
calculate_features <- function(x) {
  
  x <- as.numeric(na.omit(x))
  n <- length(x)
  
  if (n < 4L) {
    stop("calculate_features() requires at least 4 observations.")
  }
  
  s <- sd(x)
  
  if (!is.finite(s) || s <= 0) {
    s <- NA_real_
  }
  
  # Standardized sample used by several moment and tail features.
  z <- standardize_sample(x)
  
  
  #' ------------------------------------------------------------
  #' Basic location and spread
  #' ------------------------------------------------------------
  
  mean_x   <- mean(x)
  median_x <- median(x)
  
  mean_median_diff <- mean_x - median_x
  mean_to_median   <- safe_ratio(mean_x, median_x)
  
  iqr_sd_ratio <- safe_ratio(IQR(x), s) / 1.3490
  
  mad_sd_ratio <- safe_ratio(mad(x), s) / 0.6745
  
  mean_abs_dev_ratio <- safe_ratio(
    mean(abs(x - mean_x)),
    s
  ) / sqrt(2 / pi)
  
  geary_ratio <- safe_ratio(
    mean(abs(x - mean_x)),
    s
  )
  
  studentized_range <- safe_ratio(
    diff(range(x)),
    s
  )
  
  
  #' ------------------------------------------------------------
  #' Ordinary moment and L-moment features
  #' ------------------------------------------------------------
  
  skewness_val <- safe_calc(
    as.numeric(e1071::skewness(x))
  )
  
  kurtosis_val <- safe_calc(
    as.numeric(e1071::kurtosis(x))
  )
  
  # Higher standardized moments.
  m5_squared <- safe_calc(mean(z^5)^2)
  m6_normalized <- safe_calc(mean(z^6) / 15)
  
  lmom <- safe_calc(
    Lmoments::Lmoments(x),
    default = rep(NA_real_, 4)
  )
  
  l_skewness <- lmom[3]
  l_kurtosis <- lmom[4]
  
  
  #' ------------------------------------------------------------
  #' Classical normality-test statistics
  #' ------------------------------------------------------------
  
  test_stats <- normality_test_stats(x)
  
  
  #' ------------------------------------------------------------
  #' Entropy and spacing features
  #' ------------------------------------------------------------
  
  spacing_stats <- safe_calc(
    calc_spacing_ratio_stats(x),
    default = list(
      var_log_spacings  = NA_real_,
      skew_log_spacings = NA_real_
    )
  )
  
  
  #' ------------------------------------------------------------
  #' Tail and robust features
  #' ------------------------------------------------------------
  
  tail_weights <- safe_calc(
    calc_tail_weights_brys(x),
    default = list(
      LTW = NA_real_,
      RTW = NA_real_
    )
  )
  
  
  #' ------------------------------------------------------------
  #' Empirical characteristic function features
  #' ------------------------------------------------------------
  
  ecf_vals <- safe_calc(
    ecf_deviations(x),
    default = setNames(
      rep(NA_real_, 6),
      c(
        "ecf_re_t1", "ecf_re_t2", "ecf_re_t3",
        "ecf_im_t1", "ecf_im_t2", "ecf_im_t3"
      )
    )
  )
  
  
  #' ------------------------------------------------------------
  #' Assemble all base features
  #' ------------------------------------------------------------
  
  feature_df <- data.frame(
    
    ## Basic location and spread
    Mean                    = mean_x,
    Median                  = median_x,
    Mean_Median_Diff        = mean_median_diff,
    mean_to_median          = mean_to_median,
    Studentized_Range       = studentized_range,
    IQR_SD_Ratio            = iqr_sd_ratio,
    MAD_SD_Ratio            = mad_sd_ratio,
    MeanAbsDev_Ratio        = mean_abs_dev_ratio,
    Geary_Ratio             = geary_ratio,
    qn_robust_spread        = safe_calc(calc_qn_robust_spread(x)),
    
    ## Moment and L-moment features
    Skewness                = skewness_val,
    Kurtosis                = kurtosis_val,
    M5_Squared              = m5_squared,
    M6_Normalized           = m6_normalized,
    L_Skewness              = l_skewness,
    L_Kurtosis              = l_kurtosis,
    empirical_MR            = safe_calc(empirical_moment_ratios(x)),
    #dagostino_om            = safe_calc(dagostino_omnibus(x)),
    moors_kurtosis          = safe_calc(calc_moors_kurtosis(x)),
    
    ## Classical normality-test statistics
    SW_stat                 = test_stats$SW_stat,
    SF_stat                 = test_stats$SF_stat,
    AD_stat                 = test_stats$AD_stat,
    LF_stat                 = test_stats$LF_stat,
    KS_stat                 = test_stats$KS_stat,
    JB_stat                 = test_stats$JB_stat,
    Zp_stat                 = safe_calc(calculate_zp_statistic(x)),
    
    ## Entropy and spacing features
    Vasicek_stat            = safe_calc(calculate_vasicek_kmn(x)),
    Mean_Log_Spacing        = safe_calc(mean_log_spacing(x)),
    renyi_entropy           = safe_calc(renyi_entropy(x)),
    shannon_entropy         = safe_calc(shannon_entropy_stat(x)),
    #sample_entropy          = safe_calc(sample_entropy_stat(x)),
    entropy_ratio           = safe_calc(entropy_ratio(x)),
    #negentropy              = safe_calc(negentropy_approx(x)),
    greenwood_stat          = safe_calc(greenwood_stat(x)),
    rao_spacing             = safe_calc(rao_spacing_stat(x)),
    var_log_spacings        = spacing_stats$var_log_spacings,
    skew_log_spacings       = spacing_stats$skew_log_spacings,
    
    ## Tail, outlier, and robust features
    Tail_Asymmetry          = safe_calc(calculate_tail_asymmetry(x)),
    tail_index              = safe_calc(tail_index_ratio(x)),
    pickands_tail           = safe_calc(pickands_tail_index(x)),
    mean_excess             = safe_calc(mean_excess_slope(x)),
    TWR_Normalized          = safe_calc(calculate_twr_normalized(x)),
    Excess_Tail_Prop        = safe_calc(calculate_excess_tail_prop(x)),
    Outlier_Prop            = safe_calc(calculate_outlier_proportion(x)),
    Medcouple               = safe_calc(calc_medcouple(x)),
    LTW                     = tail_weights$LTW,
    RTW                     = tail_weights$RTW,
    
    ## Q-Q and correlation features
    QQ_Corr                 = safe_calc(cor(sort(x), qnorm(ppoints(n)))),
    QQ_Resid_SD             = safe_calc(calculate_qq_resid_sd(x)),
    Pearson_Skew            = safe_calc(pearson_cor_skew(x)),
    cubic_coef              = safe_calc(qq_cubic_coef(x)),
    dewet_venter            = safe_calc(dewet_venter_stat(x)),
    ryan_joiner             = safe_calc(ryan_joiner_stat(x)),
    modified_sf             = safe_calc(modified_shapiro_francia(x)),
    biweight_mid_cor        = safe_calc(biweight_midcorrelation_norm(x)),
    
    ## ECF and distance features
    epps_pulley             = safe_calc(epps_pulley_stat(x)),
    henze_zirkler           = safe_calc(henze_zirkler_stat(x)),
    ecf_re_t1               = ecf_vals["ecf_re_t1"],
    ecf_re_t2               = ecf_vals["ecf_re_t2"],
    ecf_re_t3               = ecf_vals["ecf_re_t3"],
    ecf_im_t1               = ecf_vals["ecf_im_t1"],
    ecf_im_t2               = ecf_vals["ecf_im_t2"],
    ecf_im_t3               = ecf_vals["ecf_im_t3"],
    D_ks                    = safe_calc(ks_distance(x)),
    energy_dist_to_normal   = safe_calc(energy_distance_to_normal(x)),
    
    ## Inequality and density-shape features
    lorenz_asymmetry        = safe_calc(lorenz_asymmetry_coef(x)),
    density_der_ratio       = safe_calc(density_derivative_ratio(x)),
    
    stringsAsFactors = FALSE
  )
  
  
  #' ------------------------------------------------------------
  #' Transformations of selected features
  #' ------------------------------------------------------------
  #' These are included because some relationships may be more
  #' useful after stabilizing scale or using absolute magnitude.
  #' ------------------------------------------------------------
  
  feature_df$log_vasicek_stat <- safe_calc(
    log(abs(feature_df$Vasicek_stat))
  )
  
  feature_df$log_studentized_range <- safe_calc(
    log(abs(feature_df$Studentized_Range))
  )
  
  feature_df$sq_tail_index <- safe_calc(
    sqrt(abs(feature_df$tail_index))
  )
  
  feature_df$abs_tail_asymmetry <- safe_calc(
    abs(feature_df$Tail_Asymmetry)
  )
  
  
  #' ------------------------------------------------------------
  #' Interaction features
  #' ------------------------------------------------------------
  #' These allow the model to use combinations of shape signals,
  #' rather than relying only on one feature at a time.
  #' ------------------------------------------------------------
  
  feature_df$vasicek_vs_lorenz <- safe_calc(
    feature_df$Vasicek_stat * feature_df$lorenz_asymmetry
  )
  
  feature_df$tail_asymmetry_vs_L_kurt <- safe_calc(
    feature_df$Tail_Asymmetry * feature_df$L_Kurtosis
  )
  
  feature_df$zp_stat_vs_studentized_range <- safe_calc(
    feature_df$Zp_stat * feature_df$Studentized_Range
  )
  
  feature_df$ecf_re_t3_vs_henze <- safe_calc(
    feature_df$ecf_re_t3 * feature_df$henze_zirkler
  )
  
  
  #' ------------------------------------------------------------
  #' Final cleanup
  #' ------------------------------------------------------------
  #' Some features can be undefined for small n or unusual samples.
  #' Imputation keeps the ML pipeline from failing.
  #' ------------------------------------------------------------
  
  impute_na(feature_df, method = "median")
}

#' ============================================================
#' PART 3: PREPROCESSING AND MODEL TRAINING
#' ============================================================
#' 
#' ============================================================
#' SECTION 14: FEATURE PREPROCESSING
#' ============================================================

#' Preprocess training features.
#'
#' This function prepares the feature matrix before model fitting.
#'
#' It:
#'   1. separates features from the class label;
#'   2. keeps numeric feature columns only;
#'   3. converts non-finite values to NA;
#'   4. removes features with too few observed values;
#'   5. imputes missing values;
#'   6. removes zero-variance and near-zero-variance features;
#'   7. centers and scales the remaining features.
#'
#' This is important for small n. Some features, such as
#' dagostino_om and sample_entropy, may be undefined or constant
#' when n is small.
preprocess_train_data <- function(train_data,
                                  impute_method = c("median", "mean"),
                                  scale_features = TRUE,
                                  min_nonmissing_prop = 0.50) {
  
  impute_method <- match.arg(impute_method)
  
  #' ------------------------------------------------------------
  #' Separate features and labels
  #' ------------------------------------------------------------
  
  y_train <- train_data$Label
  
  x_train <- train_data[
    ,
    setdiff(names(train_data), "Label"),
    drop = FALSE
  ]
  
  ## Keep numeric features only.
  x_train <- x_train[
    ,
    vapply(x_train, is.numeric, logical(1L)),
    drop = FALSE
  ]
  
  
  #' ------------------------------------------------------------
  #' Convert non-finite values to NA
  #' ------------------------------------------------------------
  #' caret does not handle Inf, -Inf, or NaN cleanly. Treat them
  #' as missing values before any screening or imputation.
  
  x_train[] <- lapply(x_train, function(z) {
    z <- as.numeric(z)
    z[!is.finite(z)] <- NA_real_
    z
  })
  
  
  #' ------------------------------------------------------------
  #' Drop columns with too few usable values
  #' ------------------------------------------------------------
  
  nonmissing_prop <- vapply(
    x_train,
    function(z) mean(!is.na(z)),
    numeric(1L)
  )
  
  too_missing <- nonmissing_prop < min_nonmissing_prop
  
  if (any(too_missing)) {
    message(
      "Dropping mostly-missing feature(s): ",
      paste(names(x_train)[too_missing], collapse = ", ")
    )
    
    x_train <- x_train[, !too_missing, drop = FALSE]
  }
  
  
  #' ------------------------------------------------------------
  #' Drop columns with zero variance before imputation
  #' ------------------------------------------------------------
  
  zero_var_before <- vapply(
    x_train,
    function(z) {
      z <- z[!is.na(z)]
      
      if (length(z) == 0L) {
        return(TRUE)
      }
      
      length(unique(z)) <= 1L
    },
    logical(1L)
  )
  
  if (any(zero_var_before)) {
    message(
      "Dropping zero-variance feature(s): ",
      paste(names(x_train)[zero_var_before], collapse = ", ")
    )
    
    x_train <- x_train[, !zero_var_before, drop = FALSE]
  }
  
  
  #' ------------------------------------------------------------
  #' Simple imputation
  #' ------------------------------------------------------------
  #' Do this manually so we can inspect and drop bad columns before
  #' applying center/scale.
  
  for (nm in names(x_train)) {
    
    z <- x_train[[nm]]
    
    if (!anyNA(z)) {
      next
    }
    
    fill_value <- switch(
      impute_method,
      median = median(z, na.rm = TRUE),
      mean   = mean(z, na.rm = TRUE)
    )
    
    if (!is.finite(fill_value)) {
      fill_value <- 0
    }
    
    z[is.na(z)] <- fill_value
    x_train[[nm]] <- z
  }
  
  
  #' ------------------------------------------------------------
  #' Drop zero-variance or near-zero-variance features after imputation
  #' ------------------------------------------------------------
  
  nzv_index <- caret::nearZeroVar(
    x_train,
    saveMetrics = FALSE
  )
  
  dropped_nzv <- character(0)
  
  if (length(nzv_index) > 0L) {
    dropped_nzv <- names(x_train)[nzv_index]
    
    message(
      "Dropping near-zero-variance feature(s): ",
      paste(dropped_nzv, collapse = ", ")
    )
    
    x_train <- x_train[, -nzv_index, drop = FALSE]
  }
  
  
  #' ------------------------------------------------------------
  #' Center and scale remaining features
  #' ------------------------------------------------------------
  
  if (scale_features) {
    
    prep_obj <- caret::preProcess(
      x_train,
      method = c("center", "scale")
    )
    
    x_train_std <- predict(prep_obj, x_train)
    
  } else {
    
    prep_obj <- NULL
    x_train_std <- x_train
  }
  
  x_train_std$Label <- y_train
  
  
  #' ------------------------------------------------------------
  #' Return preprocessing information
  #' ------------------------------------------------------------
  
  list(
    train_data            = x_train_std,
    prep_obj              = prep_obj,
    feature_names         = names(x_train),
    dropped_mostly_missing = names(nonmissing_prop)[too_missing],
    dropped_zero_var       = names(zero_var_before)[zero_var_before],
    dropped_near_zero_var  = dropped_nzv
  )
}

#' ============================================================
#' SECTION 15: TRAINING DATA GENERATION
#' ============================================================
#' 
#' Generate training data from paired Normal/non-Normal specs.
#'
#' Each alternative distribution is paired with its own
#' moment-matched Normal distribution.
#'
#' For each pair:
#'   - num_sim Normal samples are generated;
#'   - num_sim non-Normal samples are generated;
#'   - features are computed for each sample.
#'
#' This creates a balanced binary training set.
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
      feature_set        = feature_set,
      standardize_sample = standardize_sample,
      center_by          = center_by
    )
    
    row_id <- row_id + 1L
    
    rows[[row_id]] <- generate_feature_rows(
      n                  = n,
      n_rep              = num_sim,
      spec               = pair$alt,
      label              = "Non_Normal",
      feature_set        = feature_set,
      standardize_sample = standardize_sample,
      center_by          = center_by
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
#' 
#' caret cross-validation settings.
#'
#' The positive class is "Non_Normal", so the ROC metric measures
#' how well the model separates non-Normal samples from Normal
#' samples.
make_cv_control <- function(k_folds = 10) {
  
  caret::trainControl(
    method          = "cv",
    number          = k_folds,
    classProbs      = TRUE,
    summaryFunction = caret::twoClassSummary,
    savePredictions = "final"
  )
}

#' Train one classifier.
#'
#' This wrapper keeps the main training function clean. Each model
#' has a modest tuning grid so the code is usable without becoming
#' too slow.
train_one_model <- function(model_name,
                            train_data,
                            cv_ctrl) {
  
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
            alpha  = c(0, 0.5, 1),
            lambda = 10^seq(-4, 1, length.out = 20)
          ),
          family = "binomial"
        )
      )
    )
  }
  
  if (model_name == "RF") {
    
    return(
      caret::train(
        Label ~ .,
        data       = train_data,
        method     = "rf",
        metric     = "ROC",
        trControl  = cv_ctrl,
        tuneGrid   = expand.grid(mtry = c(2, 5, 10, 15)),
        ntree      = 500,
        importance = TRUE
      )
    )
  }
  
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
        MaxNWts = 5000,
        maxit   = 300,
        trace   = FALSE
      )
    )
  }
  
  if (model_name == "GBM") {
    
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
          n.minobsinnode    = c(5, 10)
        ),
        verbose = FALSE
      )
    )
  }
  
  if (model_name == "SVM") {
    
    return(
      caret::train(
        Label ~ .,
        data      = train_data,
        method    = "svmRadialCost",
        metric    = "ROC",
        trControl = cv_ctrl,
        
        ## tune only C; sigma is estimated internally
        tuneGrid = expand.grid(
          C = c(0.1, 0.5, 1, 5, 10)
        )
      )
    )
  }
  
  if (model_name == "KNN") {
    
    return(
      caret::train(
        Label ~ .,
        data      = train_data,
        method    = "knn",
        metric    = "ROC",
        trControl = cv_ctrl,
        tuneGrid  = expand.grid(k = c(3, 5, 7, 9, 11))
      )
    )
  }
  
  stop("Unsupported model name: ", model_name)
}


#' ============================================================
#' SECTION 17: TRAIN MODELS FOR ONE SAMPLE SIZE
#' ============================================================
#' 
#' Train all selected classifiers for one sample size.
#'
#' Output is a model bundle containing:
#'   - sample size;
#'   - trained caret models;
#'   - preprocessing object;
#'   - names of the feature columns used;
#'   - paired distribution specs used for training.
train_models_one_n <- function(n,
                               paired_specs,
                               num_sim = 100,
                               models_to_train = c("LR", "RF", "ANN", "GBM", "SVM", "KNN"),
                               feature_set = NULL,
                               standardize_sample = FALSE,
                               center_by = "mean",
                               k_folds = 10) {
  
  valid_models <- c("LR", "RF", "ANN", "GBM", "SVM", "KNN")
  
  bad_models <- setdiff(models_to_train, valid_models)
  
  if (length(bad_models) > 0L) {
    stop(
      "Unsupported model(s): ",
      paste(bad_models, collapse = ", ")
    )
  }
  
  cat("\n")
  cat("Training models for n =", n, "\n")
  cat("Distributions:", length(paired_specs), "matched pairs\n")
  cat("Replicates per class per pair:", num_sim, "\n")
  
  raw_train <- make_training_data(
    n                  = n,
    paired_specs       = paired_specs,
    num_sim            = num_sim,
    feature_set        = feature_set,
    standardize_sample = standardize_sample,
    center_by          = center_by
  )
  
  prep <- preprocess_train_data(raw_train)
  
  train_std <- prep$train_data
  
  cv_ctrl <- make_cv_control(k_folds = k_folds)
  
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
    n                     = n,
    models                = models,
    prep_obj              = prep$prep_obj,
    feature_names         = prep$feature_names,
    dropped_mostly_missing = prep$dropped_mostly_missing,
    dropped_zero_var       = prep$dropped_zero_var,
    dropped_near_zero_var  = prep$dropped_near_zero_var,
    paired_specs           = paired_specs
  )
}

#' ============================================================
#' SECTION 18: TRAIN MODELS ACROSS SAMPLE SIZES
#' ============================================================
#' 
#' Train models for all requested sample sizes.
#'
#' This is the main training wrapper used by run_normality_framework().
#' It returns a named list, with one model bundle per sample size.
train_models_all_n <- function(sample_sizes,
                               paired_specs,
                               num_sim = 100,
                               models_to_train = c("LR", "RF", "ANN", "GBM", "SVM", "KNN"),
                               feature_set = NULL,
                               standardize_sample = FALSE,
                               center_by = "mean",
                               k_folds = 10) {
  
  trained <- vector("list", length(sample_sizes))
  names(trained) <- as.character(sample_sizes)
  
  for (n in sample_sizes) {
    
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
#' 
#' ============================================================
#' SECTION 19: PREDICTION HELPERS
#' ============================================================
#' Prepare one sample for prediction.
#'
#' This function:
#'   1. computes all features for one sample;
#'   2. keeps only the feature columns retained during training;
#'   3. applies the training preprocessing object.
#'
#' This is important because some features may be dropped for one
#' sample size but kept for another.
prepare_prediction_features <- function(x, model_bundle) {
  
  feats <- calculate_features(x)
  
  missing_features <- setdiff(model_bundle$feature_names, names(feats))
  
  if (length(missing_features) > 0L) {
    stop(
      "The following training features are missing from calculate_features(): ",
      paste(missing_features, collapse = ", ")
    )
  }
  
  feats <- feats[, model_bundle$feature_names, drop = FALSE]
  
  if (!is.null(model_bundle$prep_obj)) {
    feats <- predict(model_bundle$prep_obj, feats)
  }
  
  feats
}

#' Predict class probabilities from all trained models.
#'
#' Returns one row containing P(Non_Normal) for each classifier.
predict_model_probs <- function(model_bundle, feats_std) {
  
  model_names <- names(model_bundle$models)
  
  probs <- vapply(
    model_names,
    function(model_name) {
      pred <- predict(
        model_bundle$models[[model_name]],
        newdata = feats_std,
        type    = "prob"
      )
      
      as.numeric(pred[, "Non_Normal"])
    },
    numeric(1L)
  )
  
  probs
}

#' Convert a probability to a class label.
#'
#' The positive class is "Non_Normal".
class_from_prob <- function(prob, threshold = 0.50) {
  
  ifelse(prob >= threshold, "Non_Normal", "Normal")
}

#' ============================================================
#' SECTION 20: EVALUATION DATA GENERATION
#' ============================================================
#' 
#' Create moment-matched evaluation specifications.
#'
#' Given a list of alternative distributions, this function returns
#' a list containing both the alternative distributions and their
#' matched Normal counterparts.
#'
#' The output is flattened so evaluation includes both classes.
make_eval_specs <- function(alt_specs) {
  
  pairs <- make_paired_specs(alt_specs)
  
  eval_specs <- list()
  
  for (pair in pairs) {
    
    eval_specs[[length(eval_specs) + 1L]] <- list(
      spec  = pair$normal,
      label = "Normal"
    )
    
    eval_specs[[length(eval_specs) + 1L]] <- list(
      spec  = pair$alt,
      label = "Non_Normal"
    )
  }
  
  eval_specs
}

#' ============================================================
#' SECTION 21: MODEL EVALUATION FOR ONE SAMPLE SIZE
#' ============================================================
#' 
#' Evaluate trained models for one sample size.
#'
#' For each evaluation distribution:
#'   - draw n_iter samples;
#'   - compute features;
#'   - predict P(Non_Normal) for each model;
#'   - store predictions for metrics and ROC curves.
#'
#' MajorityVote is based on the average P(Non_Normal) across
#' individual classifiers.
evaluate_one_n <- function(n,
                           model_bundle,
                           eval_specs,
                           n_iter = 1000,
                           threshold = 0.50,
                           majority_threshold = 0.50,
                           standardize_sample = FALSE,
                           center_by = "mean") {
  
  model_names <- names(model_bundle$models)
  all_names   <- c(model_names, "MajorityVote")
  
  pred_store <- lapply(all_names, function(.) {
    data.frame(
      True_Class      = character(0),
      Predicted_Class = character(0),
      Prob_Non_Normal = numeric(0),
      Distribution    = character(0),
      stringsAsFactors = FALSE
    )
  })
  
  names(pred_store) <- all_names
  
  holdout_features <- vector(
    "list",
    length(eval_specs) * n_iter
  )
  
  row_id <- 0L
  
  cat("\n")
  cat("Evaluating n =", n, "\n")
  cat("Evaluation specs:", length(eval_specs), "\n")
  cat("Replicates per spec:", n_iter, "\n")
  
  for (eval_item in eval_specs) {
    
    spec       <- eval_item$spec
    true_class <- eval_item$label
    
    dist_label <- paste0(
      spec$dist,
      "(",
      paste(round(spec$par, 4), collapse = ", "),
      ")"
    )
    
    for (i in seq_len(n_iter)) {
      
      row_id <- row_id + 1L
      
      x <- generate_data(
        n           = n,
        dist        = spec$dist,
        par         = spec$par,
        standardize = standardize_sample,
        center_by   = center_by
      )
      
      feats_std <- prepare_prediction_features(
        x            = x,
        model_bundle = model_bundle
      )
      
      holdout_features[[row_id]] <- feats_std
      
      probs <- predict_model_probs(
        model_bundle = model_bundle,
        feats_std    = feats_std
      )
      
      ## Individual model predictions
      for (model_name in model_names) {
        
        prob_nn <- probs[model_name]
        
        pred_store[[model_name]] <- rbind(
          pred_store[[model_name]],
          data.frame(
            True_Class      = true_class,
            Predicted_Class = class_from_prob(prob_nn, threshold),
            Prob_Non_Normal = prob_nn,
            Distribution    = dist_label,
            stringsAsFactors = FALSE
          )
        )
      }
      
      ## Majority vote prediction based on average probability
      mean_prob <- mean(probs, na.rm = TRUE)
      
      pred_store[["MajorityVote"]] <- rbind(
        pred_store[["MajorityVote"]],
        data.frame(
          True_Class      = true_class,
          Predicted_Class = class_from_prob(mean_prob, majority_threshold),
          Prob_Non_Normal = mean_prob,
          Distribution    = dist_label,
          stringsAsFactors = FALSE
        )
      )
    }
  }
  
  holdout_features <- do.call(rbind, holdout_features)
  
  list(
    predictions      = pred_store,
    holdout_features = holdout_features
  )
}


#' ============================================================
#' SECTION 22: PERFORMANCE METRICS
#' ============================================================
#'
#'#' Compute classification metrics from prediction data.
#'
#' The positive class is "Non_Normal".
compute_metrics <- function(pred_df) {
  
  pred_df$True_Class <- factor(
    pred_df$True_Class,
    levels = c("Non_Normal", "Normal")
  )
  
  pred_df$Predicted_Class <- factor(
    pred_df$Predicted_Class,
    levels = c("Non_Normal", "Normal")
  )
  
  cm <- caret::confusionMatrix(
    data      = pred_df$Predicted_Class,
    reference = pred_df$True_Class,
    positive  = "Non_Normal"
  )
  
  list(
    ConfusionMatrix = cm,
    Accuracy        = as.numeric(cm$overall["Accuracy"]),
    Sensitivity     = as.numeric(cm$byClass["Sensitivity"]),
    Specificity     = as.numeric(cm$byClass["Specificity"]),
    Precision       = as.numeric(cm$byClass["Precision"]),
    F1              = as.numeric(cm$byClass["F1"]),
    Predictions     = pred_df
  )
}

#' Summarize metrics for all models at one sample size.
get_metrics_summary <- function(eval_results_n, digits = 4) {
  
  model_names <- names(eval_results_n)
  
  out <- do.call(rbind, lapply(model_names, function(model_name) {
    
    res <- eval_results_n[[model_name]]
    
    data.frame(
      Model       = model_name,
      Accuracy    = round(res$Accuracy, digits),
      Sensitivity = round(res$Sensitivity, digits),
      Specificity = round(res$Specificity, digits),
      Precision   = round(res$Precision, digits),
      F1          = round(res$F1, digits),
      stringsAsFactors = FALSE
    )
  }))
  
  rownames(out) <- NULL
  out
}

#' ============================================================
#' SECTION 23: EVALUATE MODELS ACROSS SAMPLE SIZES
#' ============================================================
#' 
#' Evaluate trained models across all sample sizes.
#'
#' Returns a named list with one element per sample size. Each
#' element contains per-model metrics plus the held-out feature
#' matrix needed later for permutation AUC importance.
evaluate_models_all_n <- function(trained_models,
                                  eval_specs,
                                  n_iter = 1000,
                                  threshold = 0.50,
                                  majority_threshold = 0.50,
                                  standardize_sample = FALSE,
                                  center_by = "mean") {
  
  results <- vector("list", length(trained_models))
  names(results) <- names(trained_models)
  
  for (n_key in names(trained_models)) {
    
    model_bundle <- trained_models[[n_key]]
    n <- model_bundle$n
    
    eval_raw <- evaluate_one_n(
      n                  = n,
      model_bundle       = model_bundle,
      eval_specs         = eval_specs,
      n_iter             = n_iter,
      threshold          = threshold,
      majority_threshold = majority_threshold,
      standardize_sample = standardize_sample,
      center_by          = center_by
    )
    
    model_results <- lapply(
      eval_raw$predictions,
      compute_metrics
    )
    
    model_results$holdout_features <- eval_raw$holdout_features
    
    results[[n_key]] <- model_results
    
    cat("\nPerformance summary for n =", n, "\n")
    print(get_metrics_summary(model_results[names(eval_raw$predictions)]))
  }
  
  results
}

#' ============================================================
#' PART 5: ROC PLOTTING AND CLASSICAL-VS-ML COMPARISON
#' ============================================================
#' 
#' ============================================================
#' SECTION 24: ROC HELPER FUNCTIONS
#' ============================================================
#' 
#' Compute AUC using the trapezoidal rule.
#'
#' This is used for classical-test ROC curves where we manually
#' compute FPR and TPR across a grid of alpha cutoffs.
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
  
  if (ensure_endpoints) {
    
    if (fpr[1] > 0 || tpr[1] > 0) {
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

#' Compute an ML ROC curve from prediction data.
#'
#' The positive class is "Non_Normal", so larger probabilities
#' should correspond to stronger evidence against Normality.
get_ml_roc <- function(pred_df) {
  
  roc_obj <- pROC::roc(
    response  = factor(
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

#' ============================================================
#' SECTION 25: ML-ONLY ROC PLOTS
#' ============================================================
#'
#' Plot ML ROC curves for one sample size.
#'
#' This uses the prediction data already stored in eval_results_n.
#' It does not rerun simulations.
plot_ml_roc_one_n <- function(eval_results_n,
                              n,
                              ml_colors = NULL,
                              main_title = NULL) {
  
  model_names <- names(eval_results_n)[
    vapply(
      eval_results_n,
      function(z) !is.null(z$Predictions),
      logical(1L)
    )
  ]
  
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
    
    roc_vals <- get_ml_roc(eval_results_n[[model_name]]$Predictions)
    
    col_m <- if (!is.null(ml_colors) && model_name %in% names(ml_colors)) {
      ml_colors[model_name]
    } else {
      "gray40"
    }
    
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
                              sample_sizes,
                              output_dir = "results",
                              ml_colors = NULL) {
  
  make_output_dir(output_dir)
  
  for (n in sample_sizes) {
    
    n_key <- as.character(n)
    
    pdf_file <- file.path(
      output_dir,
      paste0("ml_roc_n", n, ".pdf")
    )
    
    pdf(pdf_file, width = 6, height = 6)
    
    plot_ml_roc_one_n(
      eval_results_n = eval_results_by_n[[n_key]],
      n              = n,
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
#'
#' Compute one classical-test p-value.
#'
#' Each classical test returns a p-value. Small p-values indicate
#' evidence against Normality.
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
      mu <- normal_par[1]
      s  <- normal_par[2]
    }
    
    return(safe_calc(ks.test(x, "pnorm", mean = mu, sd = s)$p.value))
  }
  
  stop("Unsupported classical test: ", test)
}

#' Compute classical ROC curves for one Normal/alternative pair.
#'
#' For each classical normality test:
#'   FPR = P(reject Normality | matched Normal)
#'   TPR = P(reject Normality | alternative)
#'
#' The rejection rule is p-value < alpha.
compute_classical_roc <- function(n,
                                  normal_spec,
                                  alt_spec,
                                  tests = c("SW", "AD", "JB"),
                                  alpha_grid = seq(0, 1, by = 0.05),
                                  n_sim = 1000) {
  
  fpr <- matrix(
    NA_real_,
    nrow = length(tests),
    ncol = length(alpha_grid),
    dimnames = list(tests, NULL)
  )
  
  tpr <- fpr
  
  for (test_name in tests) {
    
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
        x           = x_norm,
        test        = test_name,
        normal_par  = normal_spec$par
      )
      
      p_alt[i] <- run_classical_test(
        x           = x_alt,
        test        = test_name,
        normal_par  = normal_spec$par
      )
    }
    
    for (j in seq_along(alpha_grid)) {
      fpr[test_name, j] <- mean(p_norm < alpha_grid[j], na.rm = TRUE)
      tpr[test_name, j] <- mean(p_alt  < alpha_grid[j], na.rm = TRUE)
    }
  }
  
  list(
    FPR   = fpr,
    TPR   = tpr,
    alpha = alpha_grid
  )
}

#' ============================================================
#' SECTION 27: ML ROC CURVES FOR ONE DISTRIBUTION PAIR
#' ============================================================
#' 
#' Compute ML ROC curves for one Normal/alternative pair.
#'
#' This reruns pair-specific simulations and obtains classifier
#' probabilities. It is used for direct classical-vs-ML ROC
#' comparisons on the same distribution pair.
compute_ml_roc_pair <- function(n,
                                normal_spec,
                                alt_spec,
                                model_bundle,
                                ml_methods = c("RF", "GBM", "ANN"),
                                n_sim = 1000) {
  
  ml_methods <- intersect(ml_methods, names(model_bundle$models))
  
  if (length(ml_methods) == 0L) {
    stop("None of the requested ML methods are available in model_bundle.")
  }
  
  pred_store <- lapply(ml_methods, function(.) {
    data.frame(
      True_Class      = character(0),
      Prob_Non_Normal = numeric(0),
      stringsAsFactors = FALSE
    )
  })
  
  names(pred_store) <- ml_methods
  
  specs <- list(
    list(spec = normal_spec, label = "Normal"),
    list(spec = alt_spec,    label = "Non_Normal")
  )
  
  for (item in specs) {
    
    spec  <- item$spec
    label <- item$label
    
    for (i in seq_len(n_sim)) {
      
      x <- generate_data(
        n    = n,
        dist = spec$dist,
        par  = spec$par
      )
      
      feats_std <- prepare_prediction_features(
        x            = x,
        model_bundle = model_bundle
      )
      
      probs <- predict_model_probs(
        model_bundle = model_bundle,
        feats_std    = feats_std
      )
      
      for (m in ml_methods) {
        
        pred_store[[m]] <- rbind(
          pred_store[[m]],
          data.frame(
            True_Class      = label,
            Prob_Non_Normal = probs[m],
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }
  
  roc_list <- lapply(pred_store, get_ml_roc)
  
  ## MajorityVote for the selected ML methods.
  mv_probs <- rowMeans(
    do.call(cbind, lapply(pred_store, `[[`, "Prob_Non_Normal")),
    na.rm = TRUE
  )
  
  mv_df <- data.frame(
    True_Class      = pred_store[[ml_methods[1]]]$True_Class,
    Prob_Non_Normal = mv_probs,
    stringsAsFactors = FALSE
  )
  
  roc_list$MajorityVote <- get_ml_roc(mv_df)
  
  roc_list
}

#' ============================================================
#' SECTION 28: COMBINED CLASSICAL-VS-ML ROC PLOT
#' ============================================================
#' 
#' Plot classical and ML ROC curves in one panel.
#'
#' Classical tests are dashed. ML classifiers are solid.
plot_classical_ml_roc <- function(classical_roc,
                                  ml_roc,
                                  main_title,
                                  test_colors = NULL,
                                  ml_colors = NULL) {
  
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
  
  ## Classical curves
  for (test_name in rownames(classical_roc$FPR)) {
    
    auc_val <- compute_auc(
      classical_roc$FPR[test_name, ],
      classical_roc$TPR[test_name, ]
    )
    
    col_t <- if (!is.null(test_colors) && test_name %in% names(test_colors)) {
      test_colors[test_name]
    } else {
      "gray40"
    }
    
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
  
  ## ML curves
  for (model_name in names(ml_roc)) {
    
    col_m <- if (!is.null(ml_colors) && model_name %in% names(ml_colors)) {
      ml_colors[model_name]
    } else {
      "gray40"
    }
    
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
    cex    = 0.80
  )
  
  invisible(NULL)
}

#' ============================================================
#' SECTION 29: RUN ROC COMPARISON ACROSS PAIRS AND SAMPLE SIZES
#' ============================================================
#' 
#' Run classical-vs-ML ROC comparison.
#'
#' This function builds pair-specific ROC curves and saves one PDF
#' per alternative distribution. Each PDF contains one panel per
#' sample size.
run_roc_comparison <- function(trained_models,
                               roc_alt_specs,
                               sample_sizes,
                               classical_tests = c("SW", "AD", "JB"),
                               ml_methods = c("RF", "GBM", "ANN"),
                               alpha_grid = seq(0, 1, by = 0.05),
                               n_sim = 1000,
                               output_dir = "results",
                               test_colors = NULL,
                               ml_colors = NULL) {
  
  make_output_dir(output_dir)
  
  roc_pairs <- make_paired_specs(roc_alt_specs)
  
  roc_classical <- list()
  roc_ml        <- list()
  
  for (pair_id in seq_along(roc_pairs)) {
    
    pair <- roc_pairs[[pair_id]]
    
    alt_name <- paste0(
      pair$alt$dist,
      "_",
      paste(pair$alt$par, collapse = "_")
    )
    
    pdf_file <- file.path(
      output_dir,
      paste0("roc_comparison_", alt_name, ".pdf")
    )
    
    pdf(
      pdf_file,
      width  = 6 * length(sample_sizes),
      height = 6
    )
    
    par(
      mfrow = c(1, length(sample_sizes)),
      mar   = c(4, 4, 3, 1),
      oma   = c(1, 1, 2, 1)
    )
    
    for (n in sample_sizes) {
      
      n_key <- as.character(n)
      
      cat("ROC comparison:", pair$alt$dist,"| n =", n, "\n")
      
      classical_res <- compute_classical_roc(
        n           = n,
        normal_spec = pair$normal,
        alt_spec    = pair$alt,
        tests       = classical_tests,
        alpha_grid  = alpha_grid,
        n_sim       = n_sim
      )
      
      ml_res <- compute_ml_roc_pair(
        n            = n,
        normal_spec  = pair$normal,
        alt_spec     = pair$alt,
        model_bundle = trained_models[[n_key]],
        ml_methods   = ml_methods,
        n_sim        = n_sim
      )
      
      roc_classical[[paste0(alt_name, "_n", n)]] <- classical_res
      roc_ml[[paste0(alt_name, "_n", n)]]        <- ml_res
      
      title_text <- paste0(
        tools::toTitleCase(gsub("_", " ", pair$alt$dist)),
        " vs Matched Normal | n = ",
        n
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
      paste0("ROC Comparison: ", tools::toTitleCase(gsub("_", " ", pair$alt$dist)),
             " vs Moment-Matched Normal"
      ),
      outer = TRUE,
      cex   = 1.1,
      font  = 2
    )
    
    dev.off()
    
    cat("Saved ROC comparison:", pdf_file, "\n")
  }
  
  list(
    classical = roc_classical,
    ml        = roc_ml
  )
}

#' ============================================================
#' PART 6: PERMUTATION AUC VARIABLE IMPORTANCE
#' ============================================================
#' 
#' ============================================================
#' SECTION 30: PERMUTATION VIP HELPERS
#' ============================================================
#' 
#' Compute AUC from true labels and predicted probabilities.
#'
#' The positive class is "Non_Normal". Higher predicted probabilities
#' should indicate stronger evidence against Normality.
compute_prob_auc <- function(true_class, prob_non_normal) {
  
  roc_obj <- pROC::roc(
    response  = factor(true_class, levels = c("Normal", "Non_Normal")),
    predictor = prob_non_normal,
    levels    = c("Normal", "Non_Normal"),
    direction = "<",
    quiet     = TRUE
  )
  
  as.numeric(pROC::auc(roc_obj))
}

#' Extract the held-out labels for one sample size.
#'
#' All classifiers are evaluated on the same held-out samples, so
#' the labels can be taken from any model's prediction data.
get_holdout_labels <- function(eval_results_n) {
  
  model_names <- names(eval_results_n)[
    vapply(
      eval_results_n,
      function(z) !is.null(z$Predictions),
      logical(1L)
    )
  ]
  
  if (length(model_names) == 0L) {
    stop("No model prediction data found in eval_results_n.")
  }
  
  eval_results_n[[model_names[1L]]]$Predictions$True_Class
}

#' Predict P(Non_Normal) for one trained model.
#'
#' This small helper keeps the permutation loop easy to read.
predict_one_model_prob <- function(model, new_data) {
  
  pred <- predict(
    model,
    newdata = new_data,
    type    = "prob"
  )
  
  as.numeric(pred[, "Non_Normal"])
}

#' ============================================================
#' SECTION 31: PERMUTATION VIP FOR ONE MODEL
#' ============================================================
#' 
#' Compute permutation AUC-drop importance for one model.
#'
#' For each feature:
#'   1. randomly permute that feature in the held-out data;
#'   2. recompute predicted probabilities;
#'   3. recompute AUC;
#'   4. record baseline AUC - permuted AUC.
#'
#' A larger AUC drop means the model relied more on that feature.
compute_vip_one_model <- function(model,
                                  model_name,
                                  holdout_x,
                                  holdout_y,
                                  n_permutations = 20,
                                  digits = 4) {
  
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
  
  auc_drop <- numeric(length(feature_names))
  names(auc_drop) <- feature_names
  
  for (feature in feature_names) {
    
    perm_auc <- numeric(n_permutations)
    
    for (b in seq_len(n_permutations)) {
      
      perm_x <- holdout_x
      
      # Break the relationship between this feature and the label.
      perm_x[[feature]] <- sample(perm_x[[feature]])
      
      perm_prob <- safe_calc(
        predict_one_model_prob(
          model    = model,
          new_data = perm_x
        ),
        default = rep(NA_real_, nrow(perm_x))
      )
      
      perm_auc[b] <- safe_calc(
        compute_prob_auc(
          true_class      = holdout_y,
          prob_non_normal = perm_prob
        )
      )
    }
    
    auc_drop[feature] <- baseline_auc - mean(perm_auc, na.rm = TRUE)
  }
  
  imp_df <- data.frame(
    Feature  = names(auc_drop),
    AUC_Drop = round(as.numeric(auc_drop), digits),
    stringsAsFactors = FALSE
  )
  
  imp_df <- imp_df[order(imp_df$AUC_Drop, decreasing = TRUE), ]
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
#' 
#' Compute permutation AUC-drop importance for all models.
#'
#' This function uses the held-out feature matrix from
#' evaluate_models_all_n(). It should be optional in the full
#' pipeline because it can be time-consuming.
compute_permutation_vip <- function(trained_models,
                                    eval_results_by_n,
                                    n_permutations = 20,
                                    models_for_vip = NULL,
                                    digits = 4) {
  
  vip_results <- vector("list", length(trained_models))
  names(vip_results) <- names(trained_models)
  
  for (n_key in names(trained_models)) {
    
    cat("\nPermutation VIP for n =", n_key, "\n")
    
    model_bundle <- trained_models[[n_key]]
    eval_n       <- eval_results_by_n[[n_key]]
    
    if (is.null(eval_n$holdout_features)) {
      stop("Missing holdout_features for n = ", n_key)
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
#' SECTION 33: SUMMARIZE PERMUTATION VIP RESULTS
#' ============================================================
#' 
#' Summarize permutation VIP across models and sample sizes.
#'
#' The summary ranks features by their mean AUC drop across all
#' available model/sample-size combinations.
#'
#' Useful columns:
#'   Mean_AUC_Drop   average importance
#'   Median_AUC_Drop robust average importance
#'   Min_AUC_Drop    worst-case contribution
#'   Pct_Positive    proportion of times the feature helped
summarize_permutation_vip <- function(vip_results,
                                      digits = 4) {
  
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
  
  vip_long <- do.call(rbind, rows)
  
  split_rows <- split(vip_long, vip_long$Feature)
  
  summary_df <- do.call(rbind, lapply(names(split_rows), function(feature) {
    
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
  }))
  
  summary_df <- summary_df[
    order(summary_df$Mean_AUC_Drop, decreasing = TRUE),
  ]
  
  rownames(summary_df) <- NULL
  
  summary_df
}

#' Select features from permutation VIP results.
#'
#' This function selects the top features from a VIP summary table.
#'
#' It works whether the VIP summary came from:
#'   - one model only, such as RF;
#'   - multiple models, such as RF, GBM, ANN;
#'   - one sample size;
#'   - multiple sample sizes.
#'
#' Selection rules:
#'   "mean"   : rank by Mean_AUC_Drop
#'   "median" : rank by Median_AUC_Drop
#'   "min"    : rank by Min_AUC_Drop, conservative stability rule
#'
#' min_pct_positive:
#'   Minimum proportion of model/sample-size combinations where
#'   the feature must have positive AUC drop.
#'
#' Example:
#'   min_pct_positive = 0.50 means the feature must help in at
#'   least half of the model/sample-size combinations.
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
  
  vip_keep <- vip_summary
  
  ## Keep features that are consistently positive enough.
  vip_keep <- vip_keep[
    vip_keep$Pct_Positive >= min_pct_positive,
    ,
    drop = FALSE
  ]
  
  ## Optional: remove features with non-positive mean AUC drop.
  if (require_positive_mean) {
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

#' Return selected features with their VIP metrics.
#'
#' This is useful for reporting and checking why each feature was
#' selected.
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
#' 
#' Plot permutation VIP for one model and one sample size.
#'
#' The plot shows the top features ranked by AUC drop.
plot_permutation_vip <- function(vip_entry,
                                 top_n = 30,
                                 main_title = NULL,
                                 bar_color = "steelblue") {
  
  if (is.null(vip_entry)) {
    stop("vip_entry is NULL.")
  }
  
  imp_df <- vip_entry$imp_df
  
  imp_df <- head(imp_df, top_n)
  
  # Reverse order so the most important feature appears at the top.
  imp_df <- imp_df[rev(seq_len(nrow(imp_df))), ]
  
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
    imp_df$AUC_Drop,
    names.arg = imp_df$Feature,
    horiz     = TRUE,
    las       = 1,
    col       = bar_color,
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
                                       sample_sizes,
                                       models = c("RF"),
                                       top_n = 30,
                                       output_dir = "results") {
  
  make_output_dir(output_dir)
  
  for (n in sample_sizes) {
    
    n_key <- as.character(n)
    
    for (model_name in models) {
      
      vip_entry <- vip_results[[n_key]][[model_name]]
      
      if (is.null(vip_entry)) {
        next
      }
      
      pdf_file <- file.path(
        output_dir,
        paste0("vip_perm_", model_name, "_n", n, ".pdf")
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
#' 
#' ============================================================
#' SECTION 35: DENSITY HELPERS
#' ============================================================
#'
#' Evaluate a theoretical density.
#'
#' This function is used only for plotting. It evaluates the
#' density of a distribution specification on a supplied x-grid.
#'
#' The distribution specification must have:
#'   spec$dist
#'   spec$par
density_value <- function(x_grid, spec) {
  
  dist <- tolower(trimws(spec$dist))
  par  <- spec$par
  
  switch(
    dist,
    
    normal = {
      dnorm(x_grid, mean = par[1], sd = par[2])
    },
    
    gumbel = {
      evd::dgumbel(x_grid, loc = par[1], scale = par[2])
    },
    
    chi_square = {
      ifelse(x_grid >= 0, dchisq(x_grid, df = par[1]), 0)
    },
    
    gamma = {
      ifelse(x_grid >= 0, dgamma(x_grid, shape = par[1], rate = par[2]), 0)
    },
    
    exponential = {
      ifelse(x_grid >= 0, dexp(x_grid, rate = par[1]), 0)
    },
    
    weibull = {
      ifelse(x_grid >= 0, dweibull(x_grid, shape = par[1], scale = par[2]), 0)
    },
    
    laplace = {
      LaplacesDemon::dlaplace(
        x_grid,
        location = par[1],
        scale    = par[2]
      )
    },
    
    beta = {
      ifelse(
        x_grid >= 0 & x_grid <= 1,
        dbeta(x_grid, shape1 = par[1], shape2 = par[2]),
        0
      )
    },
    
    uniform = {
      dunif(x_grid, min = par[1], max = par[2])
    },
    
    logistic = {
      dlogis(x_grid, location = par[1], scale = par[2])
    },
    
    lognormal = {
      ifelse(x_grid > 0, dlnorm(x_grid, meanlog = par[1], sdlog = par[2]), 0)
    },
    
    t = {
      dt(x_grid, df = par[1])
    },
    
    f = {
      ifelse(x_grid >= 0, df(x_grid, df1 = par[1], df2 = par[2]), 0)
    },
    
    pareto = {
      VGAM::dpareto(x_grid, shape = par[1], scale = par[2])
    },
    
    cauchy = {
      dcauchy(x_grid, location = par[1], scale = par[2])
    },
    
    contaminated = {
      p   <- par[1]
      mu  <- par[2]
      sd1 <- par[3]
      sd2 <- par[4]
      
      p * dnorm(x_grid, mean = mu, sd = sd1) +
        (1 - p) * dnorm(x_grid, mean = mu, sd = sd2)
    },
    
    stop("Unsupported distribution in density_value(): ", dist)
  )
}

#' Choose a useful x-axis range for density plots.
#'
#' The range combines:
#'   1. a distribution-specific range for the alternative;
#'   2. mean +/- 4 SD for the matched Normal.
#'
#' This prevents the Normal density from being cut off and keeps
#' heavy-tailed alternatives from dominating the plotting window.
density_xlim <- function(alt_spec, normal_spec) {
  
  dist <- tolower(trimws(alt_spec$dist))
  par  <- alt_spec$par
  
  normal_mu <- normal_spec$par[1]
  normal_sd <- normal_spec$par[2]
  
  alt_range <- switch(
    dist,
    
    gumbel = {
      c(par[1] - 3 * par[2], par[1] + 7 * par[2])
    },
    
    chi_square = {
      c(0, qchisq(0.999, df = par[1]))
    },
    
    gamma = {
      c(0, qgamma(0.999, shape = par[1], rate = par[2]))
    },
    
    exponential = {
      c(0, qexp(0.999, rate = par[1]))
    },
    
    weibull = {
      c(0, qweibull(0.999, shape = par[1], scale = par[2]))
    },
    
    laplace = {
      c(par[1] - 8 * par[2], par[1] + 8 * par[2])
    },
    
    beta = {
      c(-0.05, 1.05)
    },
    
    uniform = {
      c(par[1] - 0.05, par[2] + 0.05)
    },
    
    logistic = {
      c(
        qlogis(0.001, location = par[1], scale = par[2]),
        qlogis(0.999, location = par[1], scale = par[2])
      )
    },
    
    lognormal = {
      c(0, qlnorm(0.999, meanlog = par[1], sdlog = par[2]))
    },
    
    t = {
      c(qt(0.001, df = par[1]), qt(0.999, df = par[1]))
    },
    
    f = {
      c(0, qf(0.999, df1 = par[1], df2 = par[2]))
    },
    
    pareto = {
      c(par[2], VGAM::qpareto(0.999, shape = par[1], scale = par[2]))
    },
    
    cauchy = {
      c(par[1] - 8 * par[2], par[1] + 8 * par[2])
    },
    
    contaminated = {
      c(par[2] - 4 * par[4], par[2] + 4 * par[4])
    },
    
    stop("Unsupported distribution in density_xlim(): ", dist)
  )
  
  normal_range <- c(
    normal_mu - 4 * normal_sd,
    normal_mu + 4 * normal_sd
  )
  
  range(c(alt_range, normal_range), finite = TRUE)
}

#' Create a readable distribution label.
#'
#' Example:
#'   list(dist = "beta", par = c(8, 2))
#' becomes:
#'   Beta(8, 2)
dist_label <- function(spec) {
  
  dist_name <- tools::toTitleCase(gsub("_", " ", spec$dist))
  par_text  <- paste(round(spec$par, 3), collapse = ", ")
  
  paste0(dist_name, "(", par_text, ")")
}

#' ============================================================
#' SECTION 36: DENSITY PLOT FUNCTION
#' ============================================================
#' 
#' Plot alternative distributions against their matched Normals.
#'
#' This function takes paired specifications created by
#' make_paired_specs(). Each panel compares one alternative
#' distribution to its own moment-matched Normal distribution.
#'
#' The function saves a PDF and returns its path invisibly.
plot_density_pairs <- function(paired_specs,
                               output_dir = "results",
                               file_name = "density_pairs.pdf",
                               n_grid = 2000,
                               n_col = 4) {
  
  make_output_dir(output_dir)
  
  n_pairs <- length(paired_specs)
  n_row   <- ceiling(n_pairs / n_col)
  
  pdf_file <- file.path(output_dir, file_name)
  
  pdf(
    pdf_file,
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
  
  alt_col  <- "#E05C5C"
  norm_col <- "#4472C4"
  
  for (pair in paired_specs) {
    
    normal_spec <- pair$normal
    alt_spec    <- pair$alt
    
    xlim <- density_xlim(
      alt_spec    = alt_spec,
      normal_spec = normal_spec
    )
    
    x_grid <- seq(
      from       = xlim[1],
      to         = xlim[2],
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
    
    ymax <- max(d_alt, d_norm, na.rm = TRUE)
    
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
      cex.main = 0.95
    )
    
    polygon(
      c(x_grid, rev(x_grid)),
      c(d_alt, rep(0, length(x_grid))),
      col    = adjustcolor(alt_col, alpha.f = 0.18),
      border = NA
    )
    
    polygon(
      c(x_grid, rev(x_grid)),
      c(d_norm, rep(0, length(x_grid))),
      col    = adjustcolor(norm_col, alpha.f = 0.18),
      border = NA
    )
    
    lines(x_grid, d_alt,  col = alt_col,  lwd = 2)
    lines(x_grid, d_norm, col = norm_col, lwd = 2, lty = 2)
    
    abline(
      v   = normal_spec$par[1],
      col = "gray50",
      lty = 3,
      lwd = 1.3
    )
    
    legend(
      "topright",
      legend = c(dist_label(alt_spec), dist_label(normal_spec)),
      col    = c(alt_col, norm_col),
      lwd    = 2,
      lty    = c(1, 2),
      bty    = "o",
      cex    = 0.75
    )
  }
  
  mtext(
    "Alternative Distributions vs Moment-Matched Normal Distributions",
    outer = TRUE,
    cex   = 1.2,
    font  = 2,
    line  = 1.0
  )
  
  cat("Saved density plots:", pdf_file, "\n")
  
  invisible(pdf_file)
}

#' ============================================================
#' PART 8: MAIN RUN FUNCTION
#' ============================================================
#' 
#' ============================================================
#' SECTION 37: MAIN PIPELINE WRAPPER
#' ============================================================
#' 
#' Run the ML-based normality testing framework.
#'
#' This is the main entry point for the complete workflow.
#'
#' The function can:
#'   1. train or load ML models;
#'   2. evaluate models on held-out distributions;
#'   3. save ML ROC plots;
#'   4. optionally compute permutation AUC VIP;
#'   5. optionally create density plots;
#'   6. optionally run classical-vs-ML ROC comparisons.
#'
#' Expensive steps are optional so the user can debug quickly.
run_normality_framework <- function(
    
  #' ------------------------------------------------------------
  #' Core simulation settings
  #' ------------------------------------------------------------
  sample_sizes = c(10, 50),
  num_sim      = 100,
  n_iter_eval  = 1000,
  
  models_to_train = c("LR", "RF", "ANN", "GBM", "SVM", "KNN"),
  
  threshold          = 0.50,
  majority_threshold = 0.40,
  
  #' ------------------------------------------------------------
  #' Distribution settings
  #' ------------------------------------------------------------
  train_alt_specs = default_training_specs(),
  eval_alt_specs  = default_eval_specs(),
  roc_alt_specs   = default_eval_specs(),
  
  #' ------------------------------------------------------------
  #' Feature settings
  #' ------------------------------------------------------------
  feature_set = NULL,
  standardize_sample = FALSE,
  center_by = "mean",
  
  #' ------------------------------------------------------------
  #' Model training / loading
  #' ------------------------------------------------------------
  train_models = TRUE,
  model_save_path = NULL,
  k_folds = 10,
  
  #' ------------------------------------------------------------
  #' Optional steps
  #' ------------------------------------------------------------
  make_density_plots  = TRUE,
  run_permutation_vip = TRUE,
  run_roc_analysis    = TRUE,
  
  #' ------------------------------------------------------------
  #' Permutation VIP settings
  #' ------------------------------------------------------------
  n_permutations = 20,
  
  ## Use one model, e.g. c("RF"), or multiple models, e.g. c("RF", "GBM", "ANN")
  models_for_vip = c("RF"),
  
  vip_top_n = 30,
  
  ## Feature-selection rule based on the VIP summary
  vip_rank_by = c("mean", "median", "min"),
  
  ## Minimum proportion of model/sample-size combinations
  ## where the feature must have positive AUC drop.
  vip_min_pct_positive = 0,
  
  ## If TRUE, drop features whose mean AUC drop is not positive.
  vip_require_positive_mean = TRUE,
  
  #' ------------------------------------------------------------
  #' ROC settings
  #' ------------------------------------------------------------
  classical_tests = c("SW", "AD", "JB"),
  ml_methods_roc  = c("RF", "GBM", "ANN"),
  alpha_grid      = seq(0, 1, by = 0.05),
  n_sim_roc       = 1000,
  
  #' ------------------------------------------------------------
  #' Plot colors
  #' ------------------------------------------------------------
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
  
  #' ------------------------------------------------------------
  #' Output
  #' ------------------------------------------------------------
  output_dir = "results"
) {
  
  make_output_dir(output_dir)
  
  vip_rank_by <- match.arg(vip_rank_by)
  
  if (is.null(model_save_path)) {
    model_save_path <- file.path(output_dir, "trained_models.RData")
  }
  
  #' ------------------------------------------------------------
  #' Build paired distribution specifications
  #' ------------------------------------------------------------
  
  train_pairs <- make_paired_specs(train_alt_specs)
  eval_specs  <- make_eval_specs(eval_alt_specs)
  
  
  #' ------------------------------------------------------------
  #' Optional density plots
  #' ------------------------------------------------------------
  
  density_plot_file <- NULL
  
  if (make_density_plots) {
    
    cat("\n")
    cat(strrep("=", 60), "\n")
    cat("DENSITY PLOTS\n")
    cat(strrep("=", 60), "\n")
    
    density_pairs <- make_paired_specs(eval_alt_specs)
    
    density_plot_file <- plot_density_pairs(
      paired_specs = density_pairs,
      output_dir   = output_dir,
      file_name    = "density_pairs.pdf"
    )
  }
  
  
  #' ------------------------------------------------------------
  #' Step 1: Train or load models
  #' ------------------------------------------------------------
  
  cat("\n")
  cat(strrep("=", 60), "\n")
  cat("STEP 1: MODEL TRAINING\n")
  cat(strrep("=", 60), "\n")
  
  if (train_models) {
    
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
      sample_sizes,
      num_sim,
      train_alt_specs,
      feature_set,
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
    
    cat("Loaded trained models from:", model_save_path, "\n")
  }
  
  
  #' ------------------------------------------------------------
  #' Step 2: Held-out evaluation
  #' ------------------------------------------------------------
  
  cat("\n")
  cat(strrep("=", 60), "\n")
  cat("STEP 2: HELD-OUT EVALUATION\n")
  cat(strrep("=", 60), "\n")
  
  eval_results_by_n <- evaluate_models_all_n(
    trained_models      = trained_models,
    eval_specs          = eval_specs,
    n_iter              = n_iter_eval,
    threshold           = threshold,
    majority_threshold  = majority_threshold,
    standardize_sample  = standardize_sample,
    center_by           = center_by
  )
  
  eval_save_path <- file.path(output_dir, "eval_results.RData")
  
  save(
    eval_results_by_n,
    eval_alt_specs,
    file = eval_save_path
  )
  
  cat("\nSaved evaluation results to:", eval_save_path, "\n")
  
  
  #' ------------------------------------------------------------
  #' Step 3: ML-only ROC plots
  #' ------------------------------------------------------------
  
  cat("\n")
  cat(strrep("=", 60), "\n")
  cat("STEP 3: ML ROC PLOTS\n")
  cat(strrep("=", 60), "\n")
  
  save_ml_roc_plots(
    eval_results_by_n = eval_results_by_n,
    sample_sizes      = sample_sizes,
    output_dir        = output_dir,
    ml_colors         = ml_colors
  )
  
  
  #' ------------------------------------------------------------
  #' Step 4: Optional permutation VIP
  #' ------------------------------------------------------------
  
  vip_results              <- NULL
  vip_summary              <- NULL
  top_features             <- NULL
  selected_feature_metrics <- NULL
  
  if (run_permutation_vip) {
    
    cat("\n")
    cat(strrep("=", 60), "\n")
    cat("STEP 4: PERMUTATION AUC VARIABLE IMPORTANCE\n")
    cat(strrep("=", 60), "\n")
    
    vip_results <- compute_permutation_vip(
      trained_models    = trained_models,
      eval_results_by_n = eval_results_by_n,
      n_permutations    = n_permutations,
      models_for_vip    = models_for_vip
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
  
  
  #' ------------------------------------------------------------
  #' Step 5: Optional classical-vs-ML ROC comparison
  #' ------------------------------------------------------------
  
  roc_results <- NULL
  
  if (run_roc_analysis) {
    
    cat("\n")
    cat(strrep("=", 60), "\n")
    cat("STEP 5: CLASSICAL VS ML ROC COMPARISON\n")
    cat(strrep("=", 60), "\n")
    
    roc_results <- run_roc_comparison(
      trained_models   = trained_models,
      roc_alt_specs    = roc_alt_specs,
      sample_sizes     = sample_sizes,
      classical_tests  = classical_tests,
      ml_methods       = ml_methods_roc,
      alpha_grid       = alpha_grid,
      n_sim            = n_sim_roc,
      output_dir       = output_dir,
      test_colors      = test_colors,
      ml_colors        = ml_colors
    )
    
    roc_save_path <- file.path(output_dir, "roc_comparison_results.RData")
    
    save(
      roc_results,
      roc_alt_specs,
      classical_tests,
      ml_methods_roc,
      file = roc_save_path
    )
    
    cat("\nSaved ROC comparison results to:", roc_save_path, "\n")
    
  } else {
    
    cat("\nSkipping classical-vs-ML ROC comparison because run_roc_analysis = FALSE.\n")
  }
  
  
  #' ------------------------------------------------------------
  #' Completion
  #' ------------------------------------------------------------
  
  cat("\n")
  cat(strrep("=", 60), "\n")
  cat("PIPELINE COMPLETE\n")
  cat(strrep("=", 60), "\n")
  cat("Output directory:", output_dir, "\n")
  
  invisible(
    list(
      trained_models            = trained_models,
      eval_results_by_n         = eval_results_by_n,
      vip_results               = vip_results,
      vip_summary               = vip_summary,
      top_features              = top_features,
      selected_feature_metrics  = selected_feature_metrics,
      roc_results               = roc_results,
      density_plot_file         = density_plot_file
    )
  )
}


selected_features <- c(
  ## Strongest overall RF VIP features
  "vasicek_vs_lorenz",
  "lorenz_asymmetry",
  "sq_tail_index",
  "L_Kurtosis",
  
  ## Classical / Q-Q normality signals
  "SW_stat",
  "henze_zirkler",
  "IQR_SD_Ratio",
  
  ## Tail / spacing / shape complements
  "Tail_Asymmetry",
  "Studentized_Range",
  "ecf_re_t3_vs_henze"
)

# 
# results_rf_selected <- run_normality_framework(
#   sample_sizes        = c(10, 20, 30, 40, 50),
#   #sample_sizes        = c(10, 50),
#   num_sim             = 500,
#   n_iter_eval         = 10,
#   models_to_train     = c("RF", "ANN", "GBM", "SVM"),#, "LR", "KNN"),
#   k_folds             = 10,
# 
#   ## Train using selected features only
#   feature_set         = selected_features,
# 
#   threshold           = 0.50,
#   majority_threshold  = 0.40,
# 
#   ## These are not needed for final model training
#   make_density_plots  = FALSE,
#   run_permutation_vip = FALSE,
#   run_roc_analysis    = FALSE,
# 
#   output_dir          = "results_selected_feat"
# )
# 
