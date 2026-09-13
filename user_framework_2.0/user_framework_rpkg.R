#' =============================================================================
#' USER FRAMEWORK FOR NORMALITY PRETESTING AND ADAPTIVE INFERENCE
#' =============================================================================
#'
#' This script contains the reusable R functions for simulating normality
#' pretesting procedures and evaluating their effects on downstream inference.
#' It supports one-sample, two-sample, ANOVA, and regression comparisons,
#' classical and machine-learning normality assessments, threshold selection,
#' power and Type I error evaluation, ROC analysis, and result generation.
#'
#' Main components
#' 1. Validate inputs, packages, files, and simulation controls.
#' 2. Generate and standardize data from supported distributions.
#' 3. Apply normality tests and downstream inferential procedures.
#' 4. Run Monte Carlo simulations and threshold trade-off analyses.
#' 5. Summarize, plot, and save reproducible analysis results.
#'
#' The file defines functions only and does not run an analysis when sourced.
#' =============================================================================

#' =============================================================================
#' 1. CORE UTILITIES AND VALIDATION
#' =============================================================================

#' -----------------------------------------------------------------------------
#' Return the right-hand value only when the left-hand value is NULL.
#' -----------------------------------------------------------------------------
`%||%` <- function(lhs, rhs) {
  if (!is.null(lhs)) lhs else rhs
}

#' -----------------------------------------------------------------------------
#' Stop with an informative message when a required package is unavailable.
#' -----------------------------------------------------------------------------
require_namespace <- function(package, purpose = NULL) {
  if (!requireNamespace(package, quietly = TRUE)) {
    detail <- if (is.null(purpose)) "" else paste0(" for ", purpose)
    stop(
      "Package '", package, "' is required", detail, ". ",
      "Install it with install.packages('", package, "').",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

#' -----------------------------------------------------------------------------
#' Validate that an input is one positive integer.
#' -----------------------------------------------------------------------------
validate_positive_integer <- function(x, name = deparse(substitute(x))) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0 || x != round(x)) {
    stop(name, " must be one positive integer.", call. = FALSE)
  }

  as.integer(x)
}

#' -----------------------------------------------------------------------------
#' Return TRUE only when finite normality scores all exceed the routing threshold.
#' -----------------------------------------------------------------------------
use_parametric_test <- function(decision_values, threshold) {
  decision_values <- as.numeric(decision_values)
  finite_values <- decision_values[is.finite(decision_values)]

  #' Route failed normality assessments to test 2 rather than defaulting to test 1.
  if (length(finite_values) == 0L) {
    return(FALSE)
  }

  all(finite_values > threshold)
}


#' -----------------------------------------------------------------------------
#' Resolve readable labels for internal downstream-test method names.
#' -----------------------------------------------------------------------------
resolve_method_labels <- function(methods, method_labels = NULL) {
  default_labels <- setNames(
    c("Test 1", "Test 2", "Adaptive"),
    c("test_1", "test_2", "adaptive")
  )

  labels <- default_labels[methods]
  labels[is.na(labels)] <- methods[is.na(labels)]

  if (is.null(method_labels)) {
    #' Return the completed result.
    return(labels)
  }

  if (is.null(names(method_labels))) {
    if (length(method_labels) != length(methods)) {
      stop(
        "An unnamed method_labels vector must match the number of methods.",
        call. = FALSE
      )
    }

    names(method_labels) <- methods
  }

  #' Compute recognized methods.
  recognized_methods <- intersect(methods, names(method_labels))
  labels[recognized_methods] <- method_labels[recognized_methods]
  labels
}

#' -----------------------------------------------------------------------------
#' Open a PDF device, run one plotting function, and always close the device.
#' -----------------------------------------------------------------------------
save_pdf_plot <- function(file, width, height, plot_function) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(file, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  result <- plot_function()
  invisible(result)
}

#' -----------------------------------------------------------------------------
#' Compute a two-sided Monte Carlo p-value with the standard plus-one correction.
#' -----------------------------------------------------------------------------
monte_carlo_p_value <- function(simulated_statistics, observed_statistic) {
  valid_statistics <- simulated_statistics[is.finite(simulated_statistics)]

  if (!is.finite(observed_statistic) || length(valid_statistics) == 0L) {
    return(NA_real_)
  }

  (sum(abs(valid_statistics) >= abs(observed_statistic)) + 1) /
    (length(valid_statistics) + 1)
}

#' =============================================================================
#' 2. DISTRIBUTION GENERATION AND STANDARDIZATION
#' =============================================================================

#' -----------------------------------------------------------------------------
#' Return default parameters for a supported named distribution.
#' -----------------------------------------------------------------------------
default_distribution_parameters <- function(dist) {
  switch(
    dist,
    normal       = c(0, 1),
    t            = 3,
    uniform      = c(0, 1),
    laplace      = c(0, 1),
    cauchy       = c(0, 1),
    chi_square   = 3,
    gamma        = c(3, 0.1),
    exponential  = 1,
    f            = c(6, 15),
    gumbel       = c(0, 1),
    log_gamma    = c(1.25, 1),
    weibull      = c(1, 2),
    lognormal    = c(0, 1),
    beta         = c(2, 5),
    logistic     = c(0, 1),
    pareto       = c(3, 1),
    contaminated = c(0.65, 0, 1, 5),
    stop("Unsupported distribution: '", dist, "'.", call. = FALSE)
  )
}

#' Compute the exact IQR of a symmetric contaminated-Normal mixture and cache it.
contaminated_normal_iqr <- local({
  cache <- new.env(parent = emptyenv())

  function(p, location, sd1, sd2) {
    key <- paste(format(c(p, location, sd1, sd2), digits = 16), collapse = "|")

    if (exists(key, envir = cache, inherits = FALSE)) {
      return(get(key, envir = cache, inherits = FALSE))
    }

    mixture_cdf <- function(q) {
      p * stats::pnorm(q, mean = location, sd = sd1) 
      + (1 - p) * stats::pnorm(q, mean = location, sd = sd2)
    }

    mixture_quantile <- function(probability) {
      span <- 12 * max(sd1, sd2)
      lower <- location - span
      upper <- location + span

      stats::uniroot(
        function(q) mixture_cdf(q) - probability,
        lower = lower,
        upper = upper,
        tol = 1e-10
      )$root
    }

    value <- mixture_quantile(0.75) - mixture_quantile(0.25)
    assign(key, value, envir = cache)
    value
  }
})

#' -----------------------------------------------------------------------------
#' Return the theoretical center and scale values used for standardization.
#' -----------------------------------------------------------------------------
distribution_reference_values <- function(dist, par) {
  dist <- tolower(trimws(dist))

  switch(
    dist,
    #' Normal moments.
    normal = {
      list(
        mean   = par[1],
        median = par[1],
        sd     = par[2],
        iqr    = stats::qnorm(0.75, par[1], par[2]) - stats::qnorm(0.25, par[1], par[2])
      )
    },
    #' Student t moments.
    t = {
      df <- par[1]
      list(
        mean   = if (df > 1) 0 else NA_real_,
        median = 0,
        sd     = if (df > 2) sqrt(df / (df - 2)) else NA_real_,
        iqr    = stats::qt(0.75, df) - stats::qt(0.25, df)
      )
    },
    #' Uniform moments.
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
    #' Laplace moments.
    laplace = {
      location <- par[1]
      scale <- par[2]
      list(
        mean   = location,
        median = location,
        sd     = sqrt(2) * scale,
        iqr    = 2 * scale * log(2)
      )
    },
    #' Cauchy robust moments.
    cauchy = {
      location <- par[1]
      scale <- par[2]
      list(
        mean   = NA_real_,
        median = location,
        sd     = NA_real_,
        iqr    = 2 * scale
      )
    },
    #' Chi-square moments.
    chi_square = {
      df <- par[1]
      list(
        mean   = df,
        median = stats::qchisq(0.50, df),
        sd     = sqrt(2 * df),
        iqr    = stats::qchisq(0.75, df) - stats::qchisq(0.25, df)
      )
    },
    #' Gamma moments.
    gamma = {
      shape <- par[1]
      rate <- par[2]
      list(
        mean   = shape / rate,
        median = stats::qgamma(0.50, shape = shape, rate = rate),
        sd     = sqrt(shape) / rate,
        iqr    = stats::qgamma(0.75, shape = shape, rate = rate) -
          stats::qgamma(0.25, shape = shape, rate = rate)
      )
    },
    #' Exponential moments.
    exponential = {
      rate <- par[1]
      list(
        mean   = 1 / rate,
        median = log(2) / rate,
        sd     = 1 / rate,
        iqr    = stats::qexp(0.75, rate) - stats::qexp(0.25, rate)
      )
    },
    #' F-distribution moments.
    f = {
      df1 <- par[1]
      df2 <- par[2]
      list(
        mean   = if (df2 > 2) df2 / (df2 - 2) else NA_real_,
        median = stats::qf(0.50, df1, df2),
        sd = if (df2 > 4) {
          sqrt(2 * df2^2 * (df1 + df2 - 2) / (df1 * (df2 - 2)^2 * (df2 - 4)))
        } else {
          NA_real_
        },
        iqr = stats::qf(0.75, df1, df2) - stats::qf(0.25, df1, df2)
      )
    },
    #' Gumbel moments.
    gumbel = {
      location <- par[1]
      scale <- par[2]
      q25 <- location - scale * log(-log(0.25))
      q75 <- location - scale * log(-log(0.75))
      list(
        mean   = location + scale * 0.5772156649015329,
        median = location - scale * log(log(2)),
        sd     = pi * scale / sqrt(6),
        iqr    = q75 - q25
      )
    },
    #' Log-Gamma moments.
    log_gamma = {
      shape <- par[1]
      rate <- par[2]
      q25 <- -log(stats::qgamma(0.75, shape = shape, rate = rate))
      q50 <- -log(stats::qgamma(0.50, shape = shape, rate = rate))
      q75 <- -log(stats::qgamma(0.25, shape = shape, rate = rate))
      list(
        mean   = log(rate) - digamma(shape),
        median = q50,
        sd     = sqrt(trigamma(shape)),
        iqr    = q75 - q25
      )
    },
    #' Weibull moments.
    weibull = {
      shape <- par[1]
      scale <- par[2]
      list(
        mean   = scale * gamma(1 + 1 / shape),
        median = scale * log(2)^(1 / shape),
        sd = scale * sqrt(gamma(1 + 2 / shape) - gamma(1 + 1 / shape)^2
        ),
        iqr = stats::qweibull(0.75, shape, scale) - stats::qweibull(0.25, shape, scale)
      )
    },
    #' Lognormal moments.
    lognormal = {
      meanlog <- par[1]
      sdlog <- par[2]
      list(
        mean   = exp(meanlog + sdlog^2 / 2),
        median = exp(meanlog),
        sd = sqrt((exp(sdlog^2) - 1) * exp(2 * meanlog + sdlog^2)),
        iqr = stats::qlnorm(0.75, meanlog, sdlog) - stats::qlnorm(0.25, meanlog, sdlog)
      )
    },
    #' Beta moments.
    beta = {
      a <- par[1]
      b <- par[2]
      list(
        mean   = a / (a + b),
        median = stats::qbeta(0.50, a, b),
        sd     = sqrt(a * b / ((a + b)^2 * (a + b + 1))),
        iqr    = stats::qbeta(0.75, a, b) - stats::qbeta(0.25, a, b)
      )
    },
    #' Logistic moments.
    logistic = {
      location <- par[1]
      scale <- par[2]
      list(
        mean   = location,
        median = location,
        sd     = scale * pi / sqrt(3),
        iqr    = stats::qlogis(0.75, location, scale) - stats::qlogis(0.25, location, scale)
      )
    },
    #' Pareto moments.
    pareto = {
      alpha <- par[1]
      xm <- par[2]
      list(
        mean   = if (alpha > 1) xm * alpha / (alpha - 1) else NA_real_,
        median = xm * 2^(1 / alpha),
        sd = if (alpha > 2) {
          xm * sqrt(alpha / ((alpha - 1)^2 * (alpha - 2)))
        } else {
          NA_real_
        },
        iqr = xm * 0.25^(-1 / alpha) - xm * 0.75^(-1 / alpha)
      )
    },
    #' Contaminated-Normal moments.
    contaminated = {
      p <- par[1]
      location <- par[2]
      sd1 <- par[3]
      sd2 <- par[4]
      list(
        mean   = location,
        median = location,
        sd     = sqrt(p * sd1^2 + (1 - p) * sd2^2),
        iqr    = contaminated_normal_iqr(p, location, sd1, sd2)
      )
    },
    stop("Unsupported distribution: '", dist, "'.", call. = FALSE)
  )
}

#' -----------------------------------------------------------------------------
#' Generate raw or theoretically standardized data from a supported distribution.
#'
#' The function accepts:
#' 1. a custom random-number generator;
#' 2. an empirical numeric vector; or
#' 3. the name of a supported distribution.
#'
#' When standardize = TRUE, the generated sample is centered using the requested
#' theoretical center and scaled using the theoretical SD, theoretical IQR, or
#' sample SD, in that order of preference.
#' -----------------------------------------------------------------------------
generate_data <- function(
    n,
    dist,
    par = NULL,
    center_by = c("median", "mean"),
    standardize = TRUE
) {
  
  #' Validate the requested sample size.
  n <- validate_positive_integer(n)
  
  #' Standardize and validate the centering method.
  center_by <- tolower(trimws(center_by[1L]))
  center_by <- match.arg(center_by, c("median", "mean"))
  
  #' ---------------------------------------------------------------------------
  #' 1. Generate data from the supplied input type.
  #' ---------------------------------------------------------------------------
  
  #' Custom generator: use sample-based reference values.
  if (is.function(dist)) {
    
    #' Generate the sample.
    samples <- as.numeric(dist(n))
    
    #' Estimate center and scale from the generated values.
    reference <- list(
      mean = mean(samples),
      median = stats::median(samples),
      sd = stats::sd(samples),
      iqr = stats::IQR(samples)
    )
    
    #' Record the source.
    distribution_label <- "custom_function"
    
    #' Empirical input: resample the supplied values with replacement.
  } else if (is.numeric(dist) && length(dist) > 0L) {
    
    #' Preserve the original empirical values.
    empirical_values <- as.numeric(dist)
    
    #' Draw the requested sample.
    samples <- sample(empirical_values, size = n, replace = TRUE)
    
    #' Use the full empirical vector for reference values.
    reference <- list(
      mean = mean(empirical_values),
      median = stats::median(empirical_values),
      sd = stats::sd(empirical_values),
      iqr = stats::IQR(empirical_values)
    )
    
    #' Record the source.
    distribution_label <- "empirical"
    
    #' Named distribution: use a built-in generator and theoretical reference values.
  } else if (is.character(dist) && length(dist) == 1L) {
    
    #' Standardize the distribution name.
    dist <- tolower(trimws(dist))
    
    #' Use supplied parameters or the framework defaults.
    par <- par %||% default_distribution_parameters(dist)
    
    #' Generate observations from the selected distribution.
    samples <- switch(
      dist,
      
      #' Normal.
      normal = stats::rnorm(n, mean = par[1], sd = par[2]),
      #' Student t.
      t = stats::rt(n, df = par[1]),
      #' Uniform.
      uniform = stats::runif(n, min = par[1], max = par[2]),
      #' Laplace.
      laplace = {
        require_namespace("LaplacesDemon", "Laplace data generation")
        LaplacesDemon::rlaplace(n, location = par[1], scale = par[2])
      },
      #' Cauchy.
      cauchy = stats::rcauchy(n, location = par[1], scale = par[2]),
      #' Chi-square.
      chi_square = stats::rchisq(n, df = par[1]),
      #' Gamma.
      gamma = stats::rgamma(n, shape = par[1], rate = par[2]),
      #' Exponential.
      exponential = stats::rexp(n, rate = par[1]),
      #' F.
      f = stats::rf(n, df1 = par[1], df2 = par[2]),
      #' Gumbel.
      gumbel = {
        require_namespace("evd", "Gumbel data generation")
        evd::rgumbel(n, loc = par[1], scale = par[2])
      },
      #' Negative-log Gamma transformation.
      log_gamma = {
        #' Shape 1 produces the standard Gumbel distribution.
        gamma_values <- stats::rgamma(n, shape = par[1], rate = par[2])
        #' Prevent log(0).
        -log(pmax(gamma_values, .Machine$double.xmin))
      },
      
      #' Weibull.
      weibull = stats::rweibull(n, shape = par[1], scale = par[2]),
      #' Lognormal.
      lognormal = stats::rlnorm(n, meanlog = par[1], sdlog = par[2]),
      #' Beta.
      beta = stats::rbeta(n, shape1 = par[1], shape2 = par[2]),
      #' Logistic.
      logistic = stats::rlogis(n, location = par[1], scale = par[2]),
      #' Pareto.
      pareto = {
        require_namespace("VGAM", "Pareto data generation")
        VGAM::rpareto(n, scale = par[2], shape = par[1])
      },
      #' Contaminated Normal mixture.
      contaminated = {
        #' Select the mixture component.
        component <- stats::rbinom(n, size = 1, prob = par[1])
        #' Assign the component-specific SD.
        component_sd <- ifelse(component == 1, par[3], par[4])
        #' Generate from the selected component.
        stats::rnorm(n, mean = par[2], sd = component_sd)
      },
      
      stop("Unsupported distribution: '", dist, "'.", call. = FALSE)
    )
    
    #' Retrieve theoretical center and scale values.
    reference <- distribution_reference_values(dist, par)
    #' Record the distribution name.
    distribution_label <- dist
    #' Reject unsupported input types.
  } else {
    stop(
      "dist must be a function, a nonempty numeric vector, or one character string.",
      call. = FALSE
    )
  }
  
  #' ---------------------------------------------------------------------------
  #' 2. Validate the generated sample.
  #' ---------------------------------------------------------------------------
  #' Require exactly n finite observations.
  if (length(samples) != n || any(!is.finite(samples))) {
    stop(
      "The generated sample must contain n finite numeric observations.",
      call. = FALSE
    )
  }
  
  #' Store the distribution label.
  attr(samples, "distribution") <- distribution_label
  #' Return raw observations when standardization is disabled.
  if (!isTRUE(standardize)) {
    attr(samples, "standardized") <- FALSE
    return(samples)
  }
  
  #' ---------------------------------------------------------------------------
  #' 3. Resolve the standardization center.
  #' ---------------------------------------------------------------------------
  #' Select the requested theoretical center.
  center_value <- reference[[center_by]]
  
  #' Fall back to the median when the requested center is undefined.
  if (!is.finite(center_value)) {
    warning("The requested center is undefined; using the median instead.",
      call. = FALSE
    )
    center_value <- reference$median
  }
  
  #' Stop when no finite center is available.
  if (!is.finite(center_value)) {
    stop("No finite center is available for standardization.",
      call. = FALSE
    )
  }
  
  #' ---------------------------------------------------------------------------
  #' 4. Resolve the standardization scale.
  #' ---------------------------------------------------------------------------
  #' First choice: theoretical SD.
  if (is.finite(reference$sd) && reference$sd > 0) {
    scale_factor <- reference$sd
    scale_method <- "standard_deviation"
    
    #' Second choice: theoretical IQR.
  } else if (is.finite(reference$iqr) && reference$iqr > 0) {
    scale_factor <- reference$iqr
    scale_method <- "interquartile_range"
    
    #' Final fallback: sample SD.
  } else {
    scale_factor <- stats::sd(samples)
    scale_method <- "sample_standard_deviation"
  }
  
  #' Require a valid positive scale.
  if (!is.finite(scale_factor) || scale_factor <= 0) {
    stop("No positive finite scale is available for standardization.",
      call. = FALSE
    )
  }
  
  #' ---------------------------------------------------------------------------
  #' 5. Standardize and annotate the sample.
  #' ---------------------------------------------------------------------------
  
  #' Center and scale the observations.
  standardized_sample <- (samples - center_value) / scale_factor
  
  #' Store reproducibility metadata.
  attr(standardized_sample, "distribution") <- distribution_label
  attr(standardized_sample, "center_by") <- center_by
  attr(standardized_sample, "center_value") <- center_value
  attr(standardized_sample, "theoretical_sd") <- reference$sd
  attr(standardized_sample, "theoretical_iqr") <- reference$iqr
  attr(standardized_sample, "scale_factor") <- scale_factor
  attr(standardized_sample, "scale_method") <- scale_method
  attr(standardized_sample, "standardized") <- TRUE
  
  #' Return the standardized sample.
  standardized_sample
}

#' =============================================================================
#' 3. NORMALITY TESTS, DIAGNOSTICS, AND CURVE INTEGRATION
#'
#' This section provides:
#' 1. classical and custom normality-test p-values;
#' 2. raw central-moment diagnostics; and
#' 3. numerical integration of ROC or related performance curves.
#' =============================================================================

#' -----------------------------------------------------------------------------
#' Apply a supported normality test and return one named p-value.
#'
#' The function accepts either:
#' 1. a supported character test code; or
#' 2. a user-supplied function that returns one numeric scalar.
#'
#' Nonfinite observations are removed before testing. The returned value is
#' named p.value so that all normality procedures use a consistent interface.
#' -----------------------------------------------------------------------------
generate_tests <- function(x, test, mu = 0, sigma = 1) {
  #' Validate the sample type.
  if (!is.numeric(x)) {
    stop("x must be a numeric vector.", call. = FALSE)
  }
  
  #' Retain finite observations only.
  x <- as.numeric(x[is.finite(x)])
  #' Require enough observations for normality testing.
  if (length(x) < 3L) {
    stop("x must contain at least three finite observations.", call. = FALSE)
  }
  
  #' Allow a user-defined normality test.
  if (is.function(test)) {
    #' Apply the custom test.
    result <- test(x)
    #' Require one numeric result.
    if (!is.numeric(result) || length(result) != 1L) {
      stop("A custom test must return one named numeric scalar.",
        call. = FALSE
      )
    }
    return(result)
  }
  
  #' Require one supported test code.
  if (!is.character(test) || length(test) != 1L) {
    stop("test must be one character string or a function.", call. = FALSE
    )
  }
  
  #' Standardize the test code.
  test <- toupper(trimws(test))
  #' Apply the selected normality procedure.
  switch(test,
    #' Shapiro-Wilk test.
    SW = c(p.value = stats::shapiro.test(x)$p.value),
    
    #' Shapiro-Francia test.
    SF = {require_namespace("nortest", "the Shapiro-Francia test")
      c(p.value = nortest::sf.test(x)$p.value)
    },
    
    #' Lilliefors test.
    LF = {require_namespace("nortest", "the Lilliefors test")
      c(p.value = nortest::lillie.test(x)$p.value)
    },
    
    #' Kolmogorov-Smirnov test against Normal(mu, sigma).
    KS = c(p.value = stats::ks.test(x, "pnorm", mean = mu, sd = sigma)$p.value),
    
    #' Anderson-Darling test with estimated Normal parameters.
    AD = {require_namespace("nortest", "the Anderson-Darling test")
      c(p.value = nortest::ad.test(x)$p.value)
    },
    
    #' Anderson-Darling test with specified Normal parameters.
    AD2 = {require_namespace("DescTools", "the specified-parameter Anderson-Darling test")
      c(p.value = DescTools::AndersonDarlingTest(x, null = "pnorm", mean = mu, sd = sigma)$p.value)
    },
    
    #' Cramer-von Mises test.
    CVM = {require_namespace("nortest", "the Cramer-von Mises test")
      c(p.value = nortest::cvm.test(x)$p.value)
    },
    
    #' Jarque-Bera test.
    JB = {require_namespace("tseries", "the Jarque-Bera test")
      c(p.value = tseries::jarque.bera.test(x)$p.value)
    },
    
    #' D'Agostino skewness test.
    DAG = {require_namespace("moments", "the D'Agostino skewness test")
      c(p.value = moments::agostino.test(x)$p.value)
    },
    
    #' Backward-compatible alias for D'Agostino skewness.
    DAP = {require_namespace("moments", "the D'Agostino skewness test")
      c(p.value = moments::agostino.test(x)$p.value)
    },
    
    #' Anscombe-Glynn kurtosis test.
    ANS = {require_namespace("moments", "the Anscombe-Glynn kurtosis test")
      c(p.value = moments::anscombe.test(x)$p.value)
    },
    
    #' Skewness z-test using the asymptotic Normal approximation.
    SKEW = {require_namespace("moments", "the skewness z-test")
      #' Sample skewness.
      sample_skewness <- moments::skewness(x)
      #' Standardized skewness statistic.
      z_statistic <- sample_skewness / sqrt(6 / length(x))
      #' Two-sided p-value.
      c(p.value = 2 * stats::pnorm(abs(z_statistic), lower.tail = FALSE))
    },
    
    #' Kurtosis z-test using the asymptotic Normal approximation.
    KURT = {require_namespace("moments", "the kurtosis z-test")
      #' Sample Pearson kurtosis.
      sample_kurtosis <- moments::kurtosis(x)
      #' Standardized excess-kurtosis statistic.
      z_statistic <- (sample_kurtosis - 3) / sqrt(24 / length(x))
      
      #' Two-sided p-value.
      c(p.value = 2 * stats::pnorm(abs(z_statistic), lower.tail = FALSE))
    },
    
    stop("Unknown test: ", test, ". Supported codes are: ",
      "SW, SF, LF, KS, AD, AD2, CVM, JB, DAG, DAP, ANS, SKEW, and KURT.",
      call. = FALSE
    )
  )
}


#' -----------------------------------------------------------------------------
#' Compute raw third or fourth central moments as descriptive diagnostics.
#'
#' These quantities are not converted to p-values and are not used as formal
#' tests in this function. They summarize sample asymmetry and tail magnitude
#' relative to the sample mean.
#' -----------------------------------------------------------------------------
generate_diagnostics <- function(x, diagnostic = c("MOM3", "MOM4")) {
  #' Validate the sample type.
  if (!is.numeric(x)) {
    stop("x must be a numeric vector.", call. = FALSE)
  }
  #' Retain finite observations only.
  x <- as.numeric(x[is.finite(x)])
  #' Require at least one finite value.
  if (length(x) == 0L) {
    stop("x must contain at least one finite observation.",call. = FALSE)
  }
  
  #' Standardize and validate the diagnostic code.
  diagnostic <- toupper(diagnostic[1L])
  diagnostic <- match.arg(diagnostic, c("MOM3", "MOM4"))
  
  #' Center observations around the sample mean.
  centered <- x - mean(x)
  
  #' Calculate the selected raw central moment.
  switch(
    diagnostic,
    #' Third central moment.
    MOM3 = c(moment3 = mean(centered^3)),
    #' Fourth central moment.
    MOM4 = c(moment4 = mean(centered^4))
  )
}


#' -----------------------------------------------------------------------------
#' Integrate a curve numerically using the trapezoidal rule.
#'
#' The function removes nonfinite coordinate pairs, sorts the curve by its
#' horizontal coordinate, optionally adds the conventional ROC endpoints, and
#' returns the numerical area. Partial curves may be normalized by their
#' observed horizontal range.
#' -----------------------------------------------------------------------------
compute_auc <- function(
    fpr,
    tpr,
    ensure_endpoints = TRUE,
    normalize = FALSE
) {
  
  #' Retain complete coordinate pairs.
  valid_points <- is.finite(fpr) & is.finite(tpr)
  
  #' Extract valid horizontal and vertical coordinates.
  x <- as.numeric(fpr[valid_points])
  y <- as.numeric(tpr[valid_points])
  
  #' Require at least two points for integration.
  if (length(x) < 2L) {
    return(NA_real_)
  }
  
  #' Sort the curve from left to right.
  sort_order <- order(x, y)
  x <- x[sort_order]
  y <- y[sort_order]
  
  #' Add the lower-left and upper-right ROC endpoints when requested.
  if (isTRUE(ensure_endpoints)) {
    
    #' Lower-left endpoint.
    if (x[1L] > 0 || y[1L] > 0) {
      x <- c(0, x)
      y <- c(0, y)
    }
    
    #' Upper-right endpoint.
    last_index <- length(x)
    if (x[last_index] < 1 || y[last_index] < 1) {
      x <- c(x, 1)
      y <- c(y, 1)
    }
  }
  
  #' Trapezoidal integration.
  auc <- sum(diff(x) * (head(y, -1L) + tail(y, -1L)) / 2)
  
  #' Normalize a partial curve by its observed horizontal width.
  if (isTRUE(normalize)) {
    x_range <- max(x) - min(x)
    auc <- if (is.finite(x_range) && x_range > 0) {
      auc / x_range
    } else {
      NA_real_
    }
  }
  
  #' Return the integrated area.
  auc
}


#' =============================================================================
#' 4. EFFECT-SIZE UTILITIES
#'
#' This section standardizes effect-size inputs used under the null and
#' alternative hypotheses. It supports scalar, vector, and explicit H0/H1
#' specifications while preserving names and dimensions.
#' =============================================================================

#' -----------------------------------------------------------------------------
#' Create a zero-valued null-effect object with the shape of the supplied effect.
#'
#' Numeric scalar or vector inputs are replaced with zeros of the same length.
#' Existing names are retained so downstream generators receive compatible
#' parameter structures.
#' -----------------------------------------------------------------------------
zero_like <- function(x) {
  #' Treat NULL as a scalar zero effect.
  if (is.null(x)) {
    return(0)
  }
  #' Require a numeric effect specification.
  if (!is.numeric(x)) {
    stop("effect_size must be numeric or list(H0 = ..., H1 = ...).",
      call. = FALSE
    )
  }
  #' Create matching zero values.
  zeros <- rep(0, length(x))
  #' Preserve parameter names.
  names(zeros) <- names(x)
  #' Return the null-effect object.
  zeros
}


#' -----------------------------------------------------------------------------
#' Split an effect-size specification into explicit H0 and H1 components.
#'
#' A named list may supply H0 and H1 directly. A numeric scalar or vector is
#' interpreted as the alternative effect, and a matching zero-valued null
#' effect is constructed automatically.
#' -----------------------------------------------------------------------------
split_effect_size <- function(effect_size) {
  #' Use explicitly supplied null and alternative effects.
  if (is.list(effect_size) && !is.data.frame(effect_size)) {
    #' Require both hypothesis labels.
    if (!all(c("H0", "H1") %in% names(effect_size))) {
      stop("A list effect_size must contain H0 and H1.", call. = FALSE
      )
    }
    #' Require numeric effects.
    if (!is.numeric(effect_size$H0) ||!is.numeric(effect_size$H1)) {
      stop("effect_size$H0 and effect_size$H1 must be numeric.",call. = FALSE)
    }
    #' Require matching parameter dimensions.
    if (length(effect_size$H0) != length(effect_size$H1)) {
      stop("effect_size$H0 and effect_size$H1 must have equal lengths.",
           call. = FALSE)
    }
    #' Return the explicit hypothesis effects.
    return(
      list(
        H0 = effect_size$H0,
        H1 = effect_size$H1
      )
    )
  }
  
  #' Treat a numeric input as the alternative effect.
  if (!is.numeric(effect_size)) {
    stop("effect_size must be numeric or list(H0 = ..., H1 = ...).",
      call. = FALSE
    )
  }
  #' Construct the matching null effect.
  list(
    H0 = zero_like(effect_size),
    H1 = effect_size
  )
}

#' =============================================================================
#' 5. ONE-SAMPLE DATA AND MAIN DOWNSTREAM TESTS
#'
#' This section defines the one-sample data generator and the inferential
#' procedures used in the one-sample demonstrations. The available procedures
#' include the ordinary t-test, sign test, randomized sign test, sign-flip
#' permutation test, and Box-Cox transformed t-test.
#' =============================================================================

#' -----------------------------------------------------------------------------
#' Generate a standardized one-sample dataset with the requested location shift.
#'
#' The selected distribution is centered and scaled by generate_data(). The
#' resulting observations are then multiplied by sd and shifted by effect_size.
#' Under the null hypothesis, effect_size is generally set to zero.
#' -----------------------------------------------------------------------------
onesample_data <- function(n = 10, effect_size = 0, sd = 1,dist = "exponential",
                           par = NULL,center_by = "median", ...) {
  #' Generate, scale, and shift the observations.
  effect_size + sd * generate_data(n = n, dist = dist, par = par, center_by = center_by)
}


#' -----------------------------------------------------------------------------
#' Return the named arguments required by the one-sample data generator.
#'
#' This helper translates the common simulation settings into the argument
#' structure expected by onesample_data().
#' -----------------------------------------------------------------------------
onesample_parameters <- function(n = 10, effect_size = 0.5, sd = 1, dist = "exponential",
                                 par = NULL, center_by = "median", ...) {
  #' Assemble the generator arguments.
  list(
    n = n,
    effect_size = effect_size,
    sd = sd,
    dist = dist,
    par = par,
    center_by = center_by
  )
}


#' -----------------------------------------------------------------------------
#' Return the supplied data unchanged for direct normality assessment.
#'
#' This identity helper allows raw observations to use the same normality-object
#' interface as residual-based ANOVA and regression analyses.
#' -----------------------------------------------------------------------------
raw_data <- function(data) {
  data
}

#' -----------------------------------------------------------------------------
#' Perform the ordinary one-sample t-test of a zero population mean.
#'
#' The function returns only the p-value required by the simulation framework.
#' -----------------------------------------------------------------------------
one_sample_t_test <- function(data) {
  #' Test H0: population mean = 0.
  test_result <- stats::t.test(data, mu = 0)
  #' Return the common downstream-test format.
  list(p.value = test_result$p.value)
}


#' -----------------------------------------------------------------------------
#' Perform an exact two-sided sign test for a specified population median.
#'
#' Observations equal to mu0 do not contribute a positive or negative sign and
#' are excluded. Under the null hypothesis, positive and negative signs each
#' have probability 0.5.
#' -----------------------------------------------------------------------------
sign_test <- function(data, mu0 = 0) {
  #' Center observations around the hypothesized median.
  centered_data <- as.numeric(data) - mu0
  #' Remove missing values and ties at the hypothesized median.
  centered_data <- centered_data[is.finite(centered_data) & centered_data != 0]
  #' No informative signs imply no evidence against H0.
  if (length(centered_data) == 0L) {
    return(list(p.value = 1))
  }
  
  #' Count observations above the hypothesized median.
  positive_signs <- sum(centered_data > 0)
  #' Calculate the exact two-sided binomial p-value.
  list(
    p.value = stats::binom.test(positive_signs, length(centered_data), p = 0.5)$p.value
  )
}


#' -----------------------------------------------------------------------------
#' Return a randomized exact sign-test p-value that is continuous under H0.
#'
#' The ordinary exact sign-test p-value is discrete. This version adds a
#' Uniform(0,1) fraction of the probability at the observed extremeness level,
#' producing a randomized p-value that is continuous under the null hypothesis.
#' -----------------------------------------------------------------------------
randomized_sign_pvalue <- function(data, mu0 = 0) {
  #' Center observations around the hypothesized median.
  centered_data <- as.numeric(data) - mu0
  #' Remove nonfinite values and ties.
  centered_data <- centered_data[is.finite(centered_data) & centered_data != 0]
  #' Number of informative signs.
  n_obs <- length(centered_data)
  
  #' No informative signs imply no rejection.
  if (n_obs == 0L) {
    return(list(p.value = 1))
  }
  
  #' Observed number of positive signs.
  observed_positive <- sum(centered_data > 0)
  #' Distance from the null-expected count.
  observed_deviation <- abs(observed_positive - n_obs / 2)
  #' Enumerate all possible positive-sign counts.
  possible_counts <- 0:n_obs
  possible_deviations <- abs(possible_counts - n_obs / 2)
  
  #' Null probabilities of the possible counts.
  count_probabilities <- stats::dbinom(
    possible_counts,
    size = n_obs,
    prob = 0.5
  )
  
  #' Probability of outcomes more extreme than observed.
  p_strict <- sum(
    count_probabilities[possible_deviations > observed_deviation]
  )
  
  #' Probability of outcomes equally extreme as observed.
  p_equal <- sum(
    count_probabilities[possible_deviations == observed_deviation]
  )
  
  #' Randomize within the probability mass at the observed boundary.
  list(
    p.value = p_strict + stats::runif(1) * p_equal
  )
}


#' -----------------------------------------------------------------------------
#' Perform a sign-flip permutation t-test for a specified population center.
#'
#' The procedure forms a studentized mean statistic and generates its null
#' distribution by randomly reversing the signs of centered observations.
#' This sign-flip argument requires symmetry of the population about mu0.
#' -----------------------------------------------------------------------------
one_sample_perm_test <- function(data, mu0 = 0, nresample = 1000) {
  
  #' Validate the number of sign-flip samples.
  nresample <- validate_positive_integer(nresample, "nresample")
  
  #' Retain finite observations.
  data <- as.numeric(data[is.finite(data)])
  n <- length(data)
  
  #' Require enough observations and positive sample variation.
  if (n < 3L || !is.finite(stats::sd(data)) || stats::sd(data) <= 0) {
    return(
      list(
        p.value = NA_real_,
        statistic = NA_real_
      )
    )
  }
  
  #' Observed studentized mean statistic.
  observed_statistic <- sqrt(n) * (mean(data) - mu0) / stats::sd(data)
  
  #' Generate the sign-flip null distribution.
  permuted_statistics <- replicate(nresample, {
    #' Assign an independent random sign to each centered observation.
    random_signs <- sample(c(-1, 1), n, replace = TRUE)
    
    #' Reflect observations around mu0.
    permuted_data <- mu0 + random_signs * (data - mu0)
    
    #' Recalculate the studentizing SD.
    permuted_sd <- stats::sd(permuted_data)
    #' Exclude degenerate resamples.
    if (!is.finite(permuted_sd) || permuted_sd <= 0) {
      return(NA_real_)
    }
    
    #' Permuted studentized mean statistic.
    sqrt(n) * (mean(permuted_data) - mu0) / permuted_sd
  })
  
  #' Return the Monte Carlo p-value and observed statistic.
  list(
    p.value = monte_carlo_p_value(permuted_statistics, observed_statistic),
    statistic = observed_statistic
  )
}


#' -----------------------------------------------------------------------------
#' Shift positive data, estimate the Box-Cox parameter, and transform the sample.
#'
#' Box-Cox transformations require strictly positive values. The sample and null
#' value are therefore shifted by the same constant when necessary. The
#' transformation parameter is selected by maximizing a manually evaluated
#' profile log-likelihood over a fixed lambda grid.
#'
#' The function returns NULL when a valid transformation cannot be obtained.
#' -----------------------------------------------------------------------------
boxcox_transform_data <- function(data, mu = 0) {
  
  #' Retain finite observations.
  data <- as.numeric(data)
  data <- data[is.finite(data)]
  
  #' Require enough observations for estimation.
  if (length(data) < 3L) {
    return(NULL)
  }
  
  #' A constant sample cannot support a t-test.
  if (!is.finite(stats::sd(data)) || stats::sd(data) == 0) {
    return(NULL)
  }
  
  #' Find the minimum across the data and null value.
  min_all <- min(c(data, mu), na.rm = TRUE)
  
  #' Shift both quantities when nonpositive values are present.
  shift <- 0
  
  if (min_all <= 0) {
    shift <- -min_all + 1e-6
  }
  
  #' Positive sample and corresponding shifted null value.
  x <- data + shift
  mu_shifted <- mu + shift
  
  #' Confirm that the transformation domain is valid.
  if (any(x <= 0) || mu_shifted <= 0) {
    return(NULL)
  }
  
  #' Apply the Box-Cox transformation for one lambda value.
  boxcox_transform <- function(z, lambda) {
    #' Use the logarithmic limit when lambda is approximately zero.
    if (abs(lambda) < 1e-8) {
      log(z)
    } else {
      (z^lambda - 1) / lambda
    }
  }
  
  #' Candidate transformation parameters.
  lambda_grid <- seq(-2, 2, by = 0.01)
  
  #' Evaluate the profile log-likelihood for each candidate lambda.
  loglik <- vapply(
    lambda_grid,
    function(lambda) {
      #' Transform the sample.
      y <- boxcox_transform(x, lambda)
      #' Reject invalid transformations.
      if (any(!is.finite(y))) {
        return(-Inf)
      }
      #' Residual variation around the transformed mean.
      y_bar <- mean(y)
      rss <- sum((y - y_bar)^2)
      
      #' Reject degenerate transformed samples.
      if (!is.finite(rss) || rss <= 0) {
        return(-Inf)
      }
      
      #' Box-Cox profile log-likelihood.
      n <- length(x)
      -n / 2 * log(rss / n) + (lambda - 1) * sum(log(x))
    },
    numeric(1)
  )
  
  #' Stop when no candidate transformation is valid.
  if (all(!is.finite(loglik))) {
    return(NULL)
  }
  
  #' Select the maximizing transformation parameter.
  lambda_hat <- lambda_grid[which.max(loglik)]
  
  #' Transform the sample and hypothesized center consistently.
  transformed_data <- boxcox_transform(x, lambda_hat)
  transformed_mu <- boxcox_transform(mu_shifted, lambda_hat)
  
  #' Require finite transformed variation.
  if (any(!is.finite(transformed_data)) ||
      !is.finite(stats::sd(transformed_data)) ||
      stats::sd(transformed_data) == 0) {
    return(NULL)
  }
  
  #' Return the transformed quantities and estimated parameters.
  list(
    transformed_data = transformed_data,
    transformed_mu = transformed_mu,
    lambda = lambda_hat,
    shift = shift
  )
}


#' -----------------------------------------------------------------------------
#' Perform a one-sample t-test after estimating a Box-Cox transformation.
#'
#' The observations and null value are transformed by boxcox_transform_data().
#' The function returns an unavailable p-value when the transformation or
#' transformed t-test cannot be completed.
#' -----------------------------------------------------------------------------
boxcox_t_test <- function(data, mu = 0) {
  
  #' Estimate and apply the transformation.
  bc <- boxcox_transform_data(data, mu)
  
  #' Return a missing result when transformation fails.
  if (is.null(bc)) {
    return(list(p.value = NA_real_))
  }
  
  #' Test the transformed mean against the transformed null value.
  test_result <- tryCatch(
    stats::t.test( bc$transformed_data, mu = bc$transformed_mu),
    error = function(error_condition) NULL
  )
  
  #' Handle numerical or testing failures.
  if (is.null(test_result)) {
    return(list(p.value = NA_real_))
  }
  
  #' Return inference and transformation details.
  list(
    p.value = test_result$p.value,
    statistic = unname(test_result$statistic),
    df = unname(test_result$parameter),
    lambda = bc$lambda,
    shift = bc$shift,
    transformed_mu = bc$transformed_mu,
    method = "Manual Box-Cox transformed one-sample t-test"
  )
}


#' =============================================================================
#' 6. TWO-SAMPLE DATA AND DOWNSTREAM TESTS
#'
#' This section generates two independent groups and defines the two-sample
#' procedures used by the framework. The available procedures are Welch's
#' t-test, the Mann-Whitney-Wilcoxon rank-sum test, and a label-permutation test
#' based on the Welch studentized statistic.
#' =============================================================================

#' -----------------------------------------------------------------------------
#' Generate two independent standardized samples separated by an effect size.
#'
#' Both samples use the same distribution family and parameter specification.
#' The second group is shifted by effect_size relative to the first group.
#' -----------------------------------------------------------------------------
two_sample_data <- function(
    n1 = 10,
    n2 = 10,
    mean1 = 0,
    effect_size = 0,
    sd1 = 1,
    sd2 = 1,
    dist = "exponential",
    par = NULL,
    center_by = "median",
    ...
) {
  
  #' Generate the reference group.
  group1 <- mean1 + sd1 * generate_data(n = n1, dist = dist, par = par, center_by = center_by)
  
  #' Generate the shifted comparison group.
  group2 <- mean1 + effect_size + sd2 * generate_data(n = n2, dist = dist, par = par, center_by = center_by)
  
  #' Combine the groups in long format.
  data.frame(
    group = factor(rep(c("x", "y"), c(n1, n2))),
    value = c(group1, group2)
  )
}


#' -----------------------------------------------------------------------------
#' Return the named arguments required by the two-sample data generator.
#'
#' A scalar n gives equal group sizes. A length-two vector supplies separate
#' sample sizes for the first and second groups.
#' -----------------------------------------------------------------------------
twosample_parameters <- function(
    n = 10,
    effect_size = 0.5,
    sd = 1,
    dist = "exponential",
    par = NULL,
    center_by = "median",
    ...
) {
  
  #' Use equal group sizes when one value is supplied.
  if (length(n) == 1L) {
    n1 <- n
    n2 <- n
    
    #' Use group-specific sizes when two values are supplied.
  } else if (length(n) == 2L) {
    n1 <- n[1L]
    n2 <- n[2L]
    
    #' Reject unsupported sample-size specifications.
  } else {
    stop("n must have length one or two for a two-sample analysis.", call. = FALSE)
  }
  
  #' Assemble the generator arguments.
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


#' -----------------------------------------------------------------------------
#' Perform Welch's two-sample t-test.
#'
#' The function extracts the two groups from the framework's long-format data
#' object and returns the p-value for equality of population means.
#' -----------------------------------------------------------------------------
twosample_t_test <- function(data) {
  
  #' Extract the two samples.
  x_data <- data$value[data$group == "x"]
  y_data <- data$value[data$group == "y"]
  
  #' Apply Welch's unequal-variance t-test.
  test_result <- stats::t.test(x = x_data, y = y_data)
  
  #' Return the common downstream-test format.
  list(p.value = test_result$p.value)
}


#' -----------------------------------------------------------------------------
#' Perform the two-sided Mann-Whitney-Wilcoxon rank-sum test.
#'
#' The procedure compares the two independent groups using their pooled ranks.
#' Its location interpretation requires compatible distributional shapes.
#' -----------------------------------------------------------------------------
Mann_whitney_U_test <- function(data) {
  
  #' Extract the two samples.
  x_data <- data$value[data$group == "x"]
  y_data <- data$value[data$group == "y"]
  
  #' Apply the two-sided rank-sum test.
  test_result <- stats::wilcox.test(x_data, y_data)
  
  #' Return the common downstream-test format.
  list(p.value = test_result$p.value)
}


#' -----------------------------------------------------------------------------
#' Perform a two-sample permutation test using the Welch t statistic.
#'
#' The observed statistic is compared with a null distribution obtained by
#' randomly reallocating the observed values to groups while preserving the
#' original group sizes.
#' -----------------------------------------------------------------------------
two_sample_perm_test <- function(data, nresample = 1000) {
  #' Validate the number of permutations.
  nresample <- validate_positive_integer(nresample, "nresample")
  #' Require the framework's long-format columns.
  if (!all(c("group", "value") %in% names(data))) {
    stop("data must contain group and value columns.", call. = FALSE)
  }
  
  #' Retain complete observations.
  data <- data[is.finite(data$value) & !is.na(data$group), , drop = FALSE]
  
  #' Extract the observed groups.
  x <- data$value[data$group == "x"]
  y <- data$value[data$group == "y"]
  
  #' Calculate a Welch-style studentized mean difference.
  welch_statistic <- function(x_values, y_values) {
    #' Standard error of the mean difference.
    denominator <- sqrt(stats::var(x_values) / length(x_values) + stats::var(y_values) / length(y_values)
    )
    
    #' Exclude degenerate group allocations.
    if (!is.finite(denominator) || denominator <= 0) {
      return(NA_real_)
    }
    
    #' Studentized mean difference.
    (mean(x_values) - mean(y_values)) / denominator
  }
  
  #' Observed test statistic.
  observed_statistic <- welch_statistic(x, y)
  
  #' Generate the label-permutation null distribution.
  permuted_statistics <- replicate(nresample, {
    #' Shuffle group labels while retaining the observed values.
    permuted_group <- sample(data$group, replace = FALSE)
    
    #' Reconstruct the permuted groups.
    x_permuted <- data$value[permuted_group == "x"]
    y_permuted <- data$value[permuted_group == "y"]
    
    #' Permuted Welch statistic.
    welch_statistic(x_permuted, y_permuted)
  })
  
  #' Return the Monte Carlo p-value and observed statistic.
  list(
    p.value = monte_carlo_p_value(permuted_statistics, observed_statistic),
    statistic = observed_statistic
  )
}


#' =============================================================================
#' 7. ONE-WAY ANOVA DATA AND DOWNSTREAM TESTS
#'
#' This section generates independent groups for one-way comparisons and
#' provides ordinary ANOVA, Kruskal-Wallis, and permutation-ANOVA procedures.
#' Normality assessments for these analyses are applied to model residuals.
#' =============================================================================

#' -----------------------------------------------------------------------------
#' Generate independent groups for a one-way ANOVA design.
#'
#' Each element of effect_size specifies the location shift for one group.
#' All groups use the same sample size, scale, distribution, and distribution
#' parameters.
#' -----------------------------------------------------------------------------
anova_gen_data <- function(
    n = 10,
    effect_size = c(0, 0, 0),
    sd = 1,
    dist = "normal",
    par = NULL,
    center_by = "median",
    ...
) {
  
  #' Require at least two group locations.
  if (!is.numeric(effect_size) || length(effect_size) < 2L) {
    stop("effect_size must contain at least two numeric group locations.", call. = FALSE)
  }
  
  #' Number and labels of groups.
  k <- length(effect_size)
  group_labels <- LETTERS[seq_len(k)]
  
  #' Generate each group using its specified location shift.
  values <- unlist(
    lapply(
      seq_along(effect_size),
      function(i) {
        effect_size[i] + sd * generate_data(n = n, dist = dist, par = par, center_by = center_by)
      }
    ),
    use.names = FALSE
  )
  
  #' Return the data in long format.
  data.frame(
    group = factor(rep(group_labels, each = n), levels = group_labels),
    value = values
  )
}


#' -----------------------------------------------------------------------------
#' Return the named arguments required by the one-way ANOVA data generator.
#' -----------------------------------------------------------------------------
anova_parameters <- function(
    n = 10,
    effect_size = c(0, 0, 0),
    sd = 1,
    dist = "exponential",
    par = NULL,
    center_by = "median",
    ...
) {
  
  #' Assemble the generator arguments.
  list(
    n = n,
    effect_size = effect_size,
    sd = sd,
    dist = dist,
    par = par,
    center_by = center_by
  )
}


#' -----------------------------------------------------------------------------
#' Fit a one-way ANOVA model and return its residuals for normality assessment.
#'
#' The adaptive procedure assesses the model errors through fitted residuals
#' rather than testing the pooled raw group observations.
#' -----------------------------------------------------------------------------
anova_residuals <- function(data) {
  #' Require the framework's long-format columns.
  if (!all(c("group", "value") %in% names(data))) {
    stop("data must contain group and value columns.", call. = FALSE)
  }
  #' Fit the additive group model and extract residuals.
  stats::residuals(stats::aov(value ~ group, data = data))
}


#' -----------------------------------------------------------------------------
#' Perform the ordinary one-way ANOVA F test.
#'
#' The null hypothesis states that all group means are equal.
#' -----------------------------------------------------------------------------
one_way_anova <- function(data) {
  #' Require the framework's long-format columns.
  if (!all(c("group", "value") %in% names(data))) {
    stop("data must contain group and value columns.", call. = FALSE)
  }
  
  #' Require at least two observed groups.
  if (length(unique(data$group)) < 2L) {
    stop("At least two groups are required for ANOVA.", call. = FALSE)
  }
  
  #' Fit the one-way ANOVA model.
  fitted_model <- stats::aov(value ~ group, data = data)
  #' Extract the omnibus F-test p-value.
  p_value <- summary(fitted_model)[[1L]][["Pr(>F)"]][1L]
  #' Return the common downstream-test format.
  list(p.value = p_value)
}


#' -----------------------------------------------------------------------------
#' Perform the Kruskal-Wallis test of equal group distributions.
#'
#' The procedure compares pooled ranks across groups. When group distributions
#' have comparable shapes, it is commonly interpreted as a comparison of group
#' locations.
#' -----------------------------------------------------------------------------
kruskal_wallis_test <- function(data) {
  
  #' Require the framework's long-format columns.
  if (!all(c("group", "value") %in% names(data))) {
    stop("data must contain group and value columns.", call. = FALSE)
  }
  
  #' Require at least two observed groups.
  if (length(unique(data$group)) < 2L) {
    stop("At least two groups are required for the Kruskal-Wallis test.",call. = FALSE)
  }
  
  #' Apply the rank-based omnibus test.
  test_result <- stats::kruskal.test(value ~ group, data = data)
  
  #' Return the common downstream-test format.
  list(p.value = test_result$p.value)
}


#' -----------------------------------------------------------------------------
#' Perform an approximate permutation one-way test using the coin package.
#'
#' The procedure tests the group effect by comparing the observed statistic
#' with an approximate permutation distribution based on nresample reallocations.
#' -----------------------------------------------------------------------------
permutation_anova <- function(data, nresample = 1000) {
  
  #' Require the permutation-testing package.
  require_namespace("coin", "permutation ANOVA")
  
  #' Validate the number of permutation samples.
  nresample <- validate_positive_integer(nresample, "nresample")
  
  #' Require the framework's long-format columns.
  if (!all(c("group", "value") %in% names(data))) {
    stop("data must contain group and value columns.", call. = FALSE)
  }
  
  #' Ensure that the grouping variable is categorical.
  data$group <- factor(data$group)
  
  #' Fit the approximate permutation test.
  test_result <- coin::oneway_test(value ~ group, data = data,
    distribution = coin::approximate(nresample = nresample)
  )
  
  #' Return the common downstream-test format.
  list(p.value = coin::pvalue(test_result))
}


#' =============================================================================
#' 8. REGRESSION DATA AND MAIN DOWNSTREAM TESTS
#'
#' This section generates data from a simple linear regression model and
#' provides ordinary least-squares, permutation, and rank-based slope tests.
#' Normality assessments are applied to fitted regression residuals.
#' =============================================================================


#' -----------------------------------------------------------------------------
#' Generate predictor and response data from a simple linear regression model.
#'
#' The response follows y = beta0 + beta1*x + error. The predictor and error
#' distributions may be specified independently. The legacy par argument is
#' retained as a fallback for both x_par and error_par.
#' -----------------------------------------------------------------------------
reg_data <- function(
    n = 10,
    beta0 = 0,
    beta1 = 0,
    x_dist = "exponential",
    error_sd = 1,
    dist = "normal",
    x_par = NULL,
    error_par = NULL,
    par = NULL,
    center_by = "median",
    ...
) {
  
  #' Use the legacy parameter input only when specific inputs are absent.
  if (!is.null(par)) {
    x_par <- x_par %||% par
    error_par <- error_par %||% par
  }
  
  #' Generate the predictor values.
  predictor <- generate_data(n = n, dist = x_dist, par = x_par, center_by = center_by)
  
  #' Generate and scale the model errors.
  error_term <- error_sd * generate_data( n = n, dist = dist, par = error_par, center_by = center_by)
  
  #' Construct the linear response.
  response <- beta0 + beta1 * predictor + error_term
  
  #' Return the regression data.
  data.frame(
    x = predictor,
    y = response
  )
}


#' -----------------------------------------------------------------------------
#' Return the named arguments required by the regression data generator.
#'
#' When beta1 is not supplied directly, effect_size is used as the slope under
#' the current simulation condition.
#' -----------------------------------------------------------------------------
reg_parameters <- function(
    n = 10,
    effect_size = 0,
    beta0 = 0,
    beta1 = NULL,
    x_dist = "exponential",
    error_sd = 1,
    dist = "normal",
    x_par = NULL,
    error_par = NULL,
    par = NULL,
    center_by = "median",
    ...
) {
  
  #' Use effect_size as the slope when beta1 is absent.
  beta1 <- beta1 %||% effect_size
  
  #' Assemble the generator arguments.
  list(
    n = n,
    beta0 = beta0,
    beta1 = beta1,
    x_dist = x_dist,
    error_sd = error_sd,
    dist = dist,
    x_par = x_par,
    error_par = error_par,
    par = par,
    center_by = center_by
  )
}


#' -----------------------------------------------------------------------------
#' Fit simple linear regression and return the residuals.
#'
#' The residuals are used as the normality-assessment object for regression
#' applications of the adaptive framework.
#' -----------------------------------------------------------------------------
reg_residuals <- function(data) {
  #' Fit the ordinary least-squares model and extract residuals.
  stats::residuals(stats::lm(y ~ x, data = data))
}


#' -----------------------------------------------------------------------------
#' Test whether the ordinary least-squares slope is zero.
#'
#' The function extracts the two-sided t-test p-value for the predictor
#' coefficient from the fitted simple linear regression model.
#' -----------------------------------------------------------------------------
simple_linear_reg <- function(data) {
  
  #' Fit the ordinary least-squares model.
  fitted_model <- stats::lm(y ~ x, data = data)
  
  #' Extract the slope p-value.
  p_value <- summary(fitted_model)$coefficients["x", "Pr(>|t|)"]
  
  #' Return the common downstream-test format.
  list(p.value = p_value)
}


#' -----------------------------------------------------------------------------
#' Perform a permutation test of the simple-regression slope.
#'
#' The observed least-squares slope is compared with slopes obtained after
#' randomly permuting the response values relative to the fixed predictor.
#' -----------------------------------------------------------------------------
perm_regression <- function(data, nresample = 10000) {
  #' Validate the number of permutations.
  nresample <- validate_positive_integer(nresample, "nresample")
  
  #' Retain complete predictor-response pairs.
  data <- data[is.finite(data$x) & is.finite(data$y), , drop = FALSE]
  
  #' Extract the analysis variables.
  x <- data$x
  y <- data$y
  
  #' Center the predictor once for repeated slope calculations.
  x_deviation <- x - mean(x)
  
  #' Fixed slope denominator.
  denominator <- sum(x_deviation^2)
  
  #' Require enough observations and predictor variation.
  if (length(x) < 3L || !is.finite(denominator) || denominator <= 0) {
    return(
      list(
        p.value = NA_real_,
        statistic = NA_real_
      )
    )
  }
  
  #' Observed least-squares slope.
  observed_slope <- sum(x_deviation * (y - mean(y)) ) / denominator
  
  #' Generate the permutation null distribution.
  permuted_slopes <- replicate(nresample, {
    #' Break the predictor-response association.
    permuted_y <- sample(y, replace = FALSE)
    #' Recalculate the slope.
    sum(x_deviation * (permuted_y - mean(permuted_y))) / denominator
  })
  
  #' Return the Monte Carlo p-value and observed slope.
  list(
    p.value = monte_carlo_p_value(permuted_slopes, observed_slope),
    statistic = observed_slope
  )
}


#' -----------------------------------------------------------------------------
#' Fit rank-based regression and return the slope p-value.
#'
#' The Rfit procedure estimates the regression relationship using rank-based
#' methods that are less sensitive to non-Normal errors and extreme responses.
#' -----------------------------------------------------------------------------
rank_regression <- function(data) {
  #' Require the rank-regression package.
  require_namespace("Rfit", "rank-based regression")
  
  #' Fit the rank-based regression model.
  fitted_model <- Rfit::rfit(y ~ x, data = data)
  
  #' Extract the coefficient table.
  coefficients <- summary(fitted_model)$coefficients
  
  #' Return the predictor p-value.
  list(
    p.value = coefficients["x", "p.value"]
  )
}

#' -----------------------------------------------------------------
#' (v) Permutation Regression
#' Performs permutation test for regression slope significance
#' 
perm_regression <- function(data, nresample = 10000) {
  data <- data[is.finite(data$x) & is.finite(data$y), , drop = FALSE]
  x <- data$x
  y <- data$y
  
  # Manual slope formula (algebraically identical to coef(lm(y ~ x))[["x"]])
  x_dev <- x - mean(x)
  original_slope <- sum(x_dev * (y - mean(y))) / sum(x_dev^2)
  
  perm_slopes <- replicate(nresample, {
    y_perm <- sample(y, replace = FALSE)
    sum(x_dev * (y_perm - mean(y_perm))) / sum(x_dev^2)
  })
  
  p_value <- (sum(abs(perm_slopes) >= abs(original_slope)) + 1) / (nresample + 1)
  return(list(p.value = p_value, statistic = original_slope))
}

#' -----------------------------------------------------------------
#' (vi) Rank-Based Regression
#' Performs rank-based (robust) regression using Rfit package
#' Requires Rfit package to be installed
#' 
rank_regression <- function(data) {
  # Fit rank-based regression
  result <- Rfit::rfit(y ~ x, data = data)
  # Extract p-value for slope coefficient 
  coefs <- summary(result)$coefficients
  p_value <- coefs["x", "p.value"]
  return(list(p.value = p_value))
}


#' -----------------------------------------------------------------
#' (vii) MM-Estimation Robust Regression (robustbase::lmrob)
#' Performs MM-estimation robust regression, robust to outliers and
#' heavy-tailed/contaminated error distributions
#' Requires robustbase package to be installed
#' 
mm_regression <- function(data) {
  # Fit MM-estimation robust regression
  result <- tryCatch(
    robustbase::lmrob(y ~ x, data = data),
    error = function(e) NULL,
    warning = function(w) suppressWarnings(robustbase::lmrob(y ~ x, data = data))
  )
  
  if (is.null(result) || !result$converged) {
    return(list(p.value = NA_real_))
  }
  
  # Extract p-value for slope coefficient using named indexing
  coefs <- summary(result)$coefficients
  p_value <- coefs["x", "Pr(>|t|)"]
  
  return(list(p.value = p_value))
}

#' -----------------------------------------------------------------
#' (viii) M-Estimation Robust Regression (MASS::rlm)
#' Performs M-estimation (Huber) robust regression, robust to
#' moderate outliers and non-normal error distributions
#' Requires MASS package to be installed
#' 
m_regression <- function(data) {
  # Fit M-estimation robust regression using Huber loss
  result <- tryCatch(
    MASS::rlm(y ~ x, data = data, method = "M"),
    error = function(e) NULL
  )
  
  if (is.null(result)) {
    return(list(p.value = NA_real_))
  }
  
  # rlm() does not return p-values directly; compute via t-approximation
  # using robust standard errors from the model summary
  s <- summary(result)
  coefs <- s$coefficients
  t_stat <- coefs["x", "t value"]
  df <- result$df.residual
  p_value <- 2 * pt(abs(t_stat), df = df, lower.tail = FALSE)
  
  return(list(p.value = p_value))
}


#' =============================================================================
#' 9. SAMPLE-FORMAT CONVERSION AND NORMALITY ASSESSMENT
#'
#' This section standardizes the sample layouts accepted by the framework and
#' provides a common interface for classical and custom normality assessments.
#' Supported inputs include numeric vectors, lists, matrices, wide data frames,
#' and long group-value data frames.
#' =============================================================================


#' -----------------------------------------------------------------------------
#' Convert supported sample-data layouts to a named list of numeric vectors.
#'
#' The returned object contains the extracted samples, their names, and the
#' number of samples. This common structure is used by custom normality methods
#' so that they can process different input layouts consistently.
#' -----------------------------------------------------------------------------
convert_to_sample_list <- function(data) {
  
  #' Classify each data-frame column.
  column_type <- function(x) {
    if (is.numeric(x)) {
      return("numeric")
    }
    
    if (is.character(x) || is.factor(x)) {
      return("character")
    }
    
    "other"
  }
  
  #' ---------------------------------------------------------------------------
  #' 1. Single numeric vector.
  #' ---------------------------------------------------------------------------
  #' Treat one numeric vector as one sample.
  if (is.numeric(data) && !is.matrix(data)) {
    return(
      list(
        samples = list(Sample1 = as.numeric(data)),
        sample_names = "Sample1",
        n_samples = 1L
      )
    )
  }
  #' ---------------------------------------------------------------------------
  #' 2. List of numeric vectors.
  #' ---------------------------------------------------------------------------
  #' Require every list element to be numeric.
  if (is.list(data) && !is.data.frame(data) && all(vapply(data, is.numeric, logical(1)))) {
    #' Convert each element to a plain numeric vector.
    samples <- lapply(data, as.numeric)
    #' Use supplied names when available.
    sample_names <- names(data)
    #' Create default names for an unnamed list.
    if (is.null(sample_names)) {
      sample_names <- paste0("Sample_", seq_along(samples))
    }
    #' Attach the resolved names to the sample list.
    names(samples) <- sample_names
    #' Return the standardized representation.
    return(
      list(samples = samples,
        sample_names = sample_names,
        n_samples = length(samples)
      )
    )
  }
  #' ---------------------------------------------------------------------------
  #' 3. Numeric matrix.
  #' ---------------------------------------------------------------------------
  #' Treat each matrix column as one sample.
  if (is.matrix(data) && is.numeric(data)) {
    #' Extract columns as numeric vectors.
    samples <- lapply(seq_len(ncol(data)), function(i) as.numeric(data[, i]))
    
    #' Use column names or create defaults.
    sample_names <- colnames(data) %||% paste0("Sample_", seq_along(samples))
    
    #' Attach sample names.
    names(samples) <- sample_names
    
    #' Return the standardized representation.
    return(
      list(samples = samples,
        sample_names = sample_names,
        n_samples = length(samples)
      )
    )
  }
  #' ---------------------------------------------------------------------------
  #' 4. Data frame.
  #' ---------------------------------------------------------------------------
  if (is.data.frame(data)) {
    #' Identify numeric, grouping, and unsupported columns.
    types <- vapply(data, column_type, character(1))
    numeric_columns <- types == "numeric"
    character_columns <- types == "character"
    
    #' Wide layout: each numeric column is one sample.
    if (all(numeric_columns)) {
      samples <- lapply(data, as.numeric)
      return(
        list(samples = samples,
          sample_names = names(samples),
          n_samples = length(samples)
        )
      )
    }
    
    #' Long layout: one grouping column and one numeric value column.
    if (sum(character_columns) == 1L && sum(numeric_columns) == 1L) {
      
      #' Locate the grouping and value columns.
      group_column <- names(data)[character_columns][1L]
      value_column <- names(data)[numeric_columns][1L]
      
      #' Split values into samples by group.
      samples <- lapply(split(data[[value_column]], data[[group_column]]), as.numeric)
      
      return(
        list(
          samples = samples,
          sample_names = names(samples),
          n_samples = length(samples)
        )
      )
    }
    
    #' Mixed layout: retain every numeric column as a sample.
    if (any(numeric_columns)) {
      samples <- lapply(data[numeric_columns], as.numeric)
      
      return(
        list(
          samples = samples,
          sample_names = names(samples),
          n_samples = length(samples)
        )
      )
    }
  }
  
  #' Reject unsupported layouts.
  stop("Unrecognized sample_data format. Supply a numeric vector, list of numeric ",
    "vectors, numeric matrix, numeric data frame, or group-value data frame.",
    call. = FALSE
  )
}


#' -----------------------------------------------------------------------------
#' Apply a classical normality test to one or more samples.
#'
#' The function accepts the same common sample layouts used by the framework.
#' A single p-value is calculated for each sample and returned as a named numeric
#' vector. For long-format data, group_col identifies the group variable and
#' value_col identifies the numeric observations.
#' -----------------------------------------------------------------------------
normality_test <- function(
    data,
    test = "SW",
    mu = 0,
    sigma = 1,
    group_col = 1,
    value_col = 2
) {
  
  #' Apply the selected test to one numeric sample.
  apply_test <- function(x) {
    #' Require numeric observations.
    if (!is.numeric(x)) {
      stop("Every sample must be numeric.", call. = FALSE)
    }
    #' Delegate to the common test dispatcher.
    generate_tests(x, test = test, mu = mu, sigma = sigma)
  }
  
  #' ---------------------------------------------------------------------------
  #' 1. Single numeric vector.
  #' ---------------------------------------------------------------------------
  if (is.numeric(data) && is.null(dim(data))) {
    #' Test the single sample.
    decision_values <- apply_test(data)
    #' Assign a consistent sample name.
    names(decision_values) <- "Sample1"
    #' ---------------------------------------------------------------------------
    #' 2. List of separate samples.
    #' ---------------------------------------------------------------------------
    
  } else if (is.list(data) && !is.data.frame(data)) {
    #' Apply the test to each list element.
    decision_values <- vapply(data, apply_test, numeric(1))
    #' Use supplied names or create defaults.
    names(decision_values) <- names(data) %||% paste0("Sample", seq_along(data))
    
    #' ---------------------------------------------------------------------------
    #' 3. Wide numeric data frame.
    #' ---------------------------------------------------------------------------
  } else if (is.data.frame(data) && all(vapply(data, is.numeric, logical(1)))) {
    #' Treat each numeric column as one sample.
    decision_values <- vapply(as.list(data), apply_test, numeric(1))
    #' Preserve column names.
    names(decision_values) <- names(data)
    
    #' ---------------------------------------------------------------------------
    #' 4. Long group-value data frame or matrix.
    #' ---------------------------------------------------------------------------
  } else if ((is.data.frame(data) || is.matrix(data)) && ncol(data) >= 2L) {
    
    #' Extract columns according to the input type.
    if (is.matrix(data)) {
      values <- data[, value_col]
      groups <- data[, group_col]
    } else {
      values <- data[[value_col]]
      groups <- data[[group_col]]
    }
    
    #' Convert the value column to numeric.
    values <- suppressWarnings(as.numeric(values))
    #' Reject a nonnumeric value column.
    if (all(!is.finite(values))) {
      stop("value_col must identify a numeric column.", call. = FALSE)
    }

    #' Split observations into group-specific samples.
    grouped_samples <- split(values, groups)
    #' Apply the test to each group.
    decision_values <- vapply(grouped_samples, apply_test, numeric(1))
    
    #' Reject unsupported layouts.
  } else {
    stop("Unsupported input format for normality_test().", call. = FALSE)
  }
  #' Return one named p-value per sample.
  decision_values
}


#' -----------------------------------------------------------------------------
#' Apply either a classical or custom normality assessment.
#'
#' Classical assessments use normality_test() and return p-values. Custom
#' assessments use a user-supplied function in config$fn and may return any
#' scalar score whose direction is understood by the routing procedure.
#'
#' Additional elements in config are passed to the custom function, except for
#' fn, label, and name, which are treated as framework metadata.
#' -----------------------------------------------------------------------------
assess_normality <- function(data, method = "classical", config = list()) {
  #' Standardize and validate the assessment method.
  method <- match.arg( tolower(trimws(method)), c("classical", "custom"))
  
  #' ---------------------------------------------------------------------------
  #' 1. Classical normality assessment.
  #' ---------------------------------------------------------------------------
  #' Delegate classical tests to the unified dispatcher.
  if (method == "classical") {
    return(
      normality_test(
        data,
        test = config$norm_test %||% "SW",
        mu = config$mu %||% 0,
        sigma = config$sigma %||% 1,
        group_col = config$group_col %||% 1,
        value_col = config$value_col %||% 2
      )
    )
  }
  
  #' ---------------------------------------------------------------------------
  #' 2. Custom normality assessment.
  #' ---------------------------------------------------------------------------
  #' Require a callable custom assessment function.
  if (is.null(config$fn) || !is.function(config$fn)) {
    stop("config$fn must be a function for method = 'custom'.", call. = FALSE)
  }
  
  #' Convert the input to a common sample-list structure.
  sample_information <- convert_to_sample_list(data)
  #' Store the custom function.
  custom_function <- config$fn
  #' Forward nonmetadata configuration values as function arguments.
  extra_arguments <- config[setdiff(names(config), c("fn", "label", "name"))]
  
  #' Apply the custom function separately to each sample.
  decision_values <- vapply(sample_information$samples,function(x) {
      #' Call the custom function with the sample and extra arguments.
      result <- do.call(custom_function,c(list(x), extra_arguments))

      #' Require one numeric score per sample.
      if (!is.numeric(result) || length(result) != 1L) {
        stop("config$fn must return one named numeric scalar.", call. = FALSE)
      }
      #' Remove any internal name before vapply stores the result.
      as.numeric(result)
    },
    numeric(1)
  )
  #' Restore the sample names.
  names(decision_values) <- sample_information$sample_names
  #' Return one custom decision score per sample.
  decision_values
}

#' =============================================================================
#' 10. NORMALITY ROC SIMULATION AND PLOTTING
#'
#' This section evaluates how well classical normality tests or custom
#' normality scores distinguish Normal samples from a selected non-Normal
#' alternative.
#'
#' The workflow:
#' 1. generates samples under Normality and non-Normality;
#' 2. calculates the requested test p-values or custom scores;
#' 3. converts stored values to ROC coordinates over a threshold grid;
#' 4. plots the ROC curves; and
#' 5. calculates the corresponding AUC values.
#' =============================================================================

#' -----------------------------------------------------------------------------
#' Simulate ROC curves for one or more custom normality-score functions.
#'
#' The custom function may return one score or a named vector of scores. For
#' multi-sample analyses, the function is applied separately to each sample and
#' the minimum score is retained because smaller scores indicate stronger
#' evidence against Normality.
#'
#' Samples are generated once under H0 and H1. The stored scores are then
#' evaluated over the complete threshold grid without regenerating data.
#' -----------------------------------------------------------------------------
fn_for_custom_test_roc_curve <- function(
    n = 10,
    threshold_grid = seq(0, 1, by = 0.025),
    H1_dist = "exponential",
    Nsim = 1000,
    gen_data = two_sample_data,
    get_parameters = twosample_parameters,
    fn_to_get_norm_obj = raw_data,
    custom_fn,
    custom_args = list(),
    ...
) {
  
  #' Require a callable custom scoring function.
  if (missing(custom_fn) || !is.function(custom_fn)) {
    stop("custom_fn must be supplied as a function.", call. = FALSE)
  }
  
  #' Standardize the threshold grid.
  threshold_grid <- sort(unique(as.numeric(threshold_grid)))
  #' Validate the simulation count.
  Nsim <- validate_positive_integer(Nsim, "Nsim")
  
  #' ---------------------------------------------------------------------------
  #' 1. Determine the custom-score structure.
  #' ---------------------------------------------------------------------------
  #' Generate one Normal probe dataset.
  probe_parameters <- get_parameters(n = n, dist = "normal", par = NULL, ...)
  
  probe_data <- do.call(gen_data, probe_parameters)
  
  #' Extract the sample vectors used for normality assessment.
  probe_samples <- convert_to_sample_list(fn_to_get_norm_obj(probe_data))$samples
  
  #' Apply the custom function to the first probe sample.
  probe_output <- do.call(custom_fn, c(list(probe_samples[[1L]]), custom_args))
  
  #' Require at least one numeric score.
  if (!is.numeric(probe_output) || length(probe_output) == 0L) {
    stop("custom_fn must return a nonempty numeric vector.", call. = FALSE)
  }
  
  #' Number and names of custom ROC curves.
  n_curves <- length(probe_output)
  curve_names <- names(probe_output) %||% paste0("custom_", seq_len(n_curves))
  
  #' ---------------------------------------------------------------------------
  #' 2. Define one custom-score evaluation.
  #' ---------------------------------------------------------------------------
  evaluate_custom_score <- function(data_object) {
    
    #' Convert the normality object to separate sample vectors.
    samples <- convert_to_sample_list(data_object)$samples
    
    #' Calculate all custom scores for each sample.
    score_matrix <- vapply(
      samples,
      function(sample_values) {
        #' Apply the custom scoring function.
        values <- as.numeric(do.call(custom_fn, c(list(sample_values),custom_args)))
        #' Require a consistent number of scores.
        if (length(values) != n_curves) {
          stop("custom_fn returned an inconsistent number of scores.", call. = FALSE)
        }
        
        values
      },
      numeric(n_curves)
    )
    
    #' Preserve a matrix structure for a single sample.
    if (is.null(dim(score_matrix))) {
      score_matrix <- matrix(score_matrix, nrow = n_curves)
    }
    
    #' Retain the least-Normal score across samples.
    combined_scores <- apply(score_matrix, MARGIN = 1L, FUN = min, na.rm = TRUE)
    #' Replace failed score combinations with missing values.
    combined_scores[!is.finite(combined_scores)] <- NA_real_
    #' Restore curve names.
    names(combined_scores) <- curve_names

    combined_scores
  }
  
  #' ---------------------------------------------------------------------------
  #' 3. Allocate score storage.
  #' ---------------------------------------------------------------------------
  #' Scores from Normal samples.
  scores_H0 <- matrix(NA_real_, nrow = n_curves, ncol = Nsim, dimnames = list(curve_names, NULL))
  
  #' Scores from the alternative distribution.
  scores_H1 <- matrix(NA_real_, nrow = n_curves, ncol = Nsim, dimnames = list(curve_names, NULL))
  
  #' Display simulation progress.
  cat("Computing ROC curves for custom normality scores\n")
  
  progress_bar <- utils::txtProgressBar(min = 0, max = 2 * Nsim, style = 3)
  
  on.exit(close(progress_bar),add = TRUE)
  
  #' ---------------------------------------------------------------------------
  #' 4. Simulate scores under Normality.
  #' ---------------------------------------------------------------------------
  for (simulation_index in seq_len(Nsim)) {

    #' Generate one Normal dataset.
    parameters_H0 <- get_parameters(n = n, dist = "normal", par = NULL, ...)
    
    data_H0 <- do.call(gen_data, parameters_H0)
    
    #' Extract the object used for normality assessment.
    normality_object_H0 <- fn_to_get_norm_obj(data_H0)
    
    #' Store the custom scores and retain failures as missing.
    scores_H0[, simulation_index] <- tryCatch(
      evaluate_custom_score(normality_object_H0),
      error = function(error_condition) {
        rep(NA_real_, n_curves)
      }
    )
    
    #' Update progress.
    utils::setTxtProgressBar(progress_bar, simulation_index)
  }
  
  #' ---------------------------------------------------------------------------
  #' 5. Simulate scores under the alternative.
  #' ---------------------------------------------------------------------------
  for (simulation_index in seq_len(Nsim)) {
    #' Generate one non-Normal dataset.
    parameters_H1 <- get_parameters(n = n, dist = H1_dist,...)
    
    data_H1 <- do.call(gen_data, parameters_H1)
    
    #' Extract the object used for normality assessment.
    normality_object_H1 <- fn_to_get_norm_obj(data_H1)
    
    #' Store the custom scores and retain failures as missing.
    scores_H1[, simulation_index] <- tryCatch(
      evaluate_custom_score(normality_object_H1),
      error = function(error_condition) {
        rep(NA_real_, n_curves)
      }
    )
    
    #' Update progress.
    utils::setTxtProgressBar(progress_bar, Nsim + simulation_index)
  }
  
  #' ---------------------------------------------------------------------------
  #' 6. Convert stored scores to ROC coordinates.
  #' ---------------------------------------------------------------------------
  #' Allocate false-positive rates.
  false_positive_rate <- matrix(NA_real_, nrow = n_curves, ncol = length(threshold_grid),
    dimnames = list(curve_names,paste0("t_", threshold_grid))
  )
  
  #' Use the same dimensions for true-positive rates.
  true_positive_rate <- false_positive_rate
  
  #' Calculate rejection rates over all thresholds.
  for (curve_index in seq_len(n_curves)) {
    
    #' False-positive rate under Normality.
    false_positive_rate[curve_index, ] <- vapply(
      threshold_grid,
      function(threshold) {
        mean(scores_H0[curve_index, ] <= threshold, na.rm = TRUE)
      },
      numeric(1)
    )
    
    #' True-positive rate under the alternative.
    true_positive_rate[curve_index, ] <- vapply(
      threshold_grid,
      function(threshold) {
        mean(scores_H1[curve_index, ] <= threshold, na.rm = TRUE)
      },
      numeric(1)
    )
  }
  
  #' Return ROC coordinates and the underlying scores.
  list(
    FPR = false_positive_rate,
    TPR = true_positive_rate,
    threshold = threshold_grid,
    raw = list(
      H0 = scores_H0,
      H1 = scores_H1
    )
  )
}


#' -----------------------------------------------------------------------------
#' Simulate ROC curves for one or more classical normality tests.
#'
#' Each simulated dataset is generated once and evaluated by every requested
#' test. For analyses containing multiple samples, the minimum p-value is used
#' because rejection by any sample indicates a departure from the joint
#' Normality requirement.
#'
#' The stored p-values are subsequently converted to ROC coordinates over the
#' supplied pretest-alpha grid.
#' -----------------------------------------------------------------------------
fn_for_norm_test_roc_curve <- function(
    n = 10,
    alpha_pretest = seq(0, 1, by = 0.025),
    H1_dist = "exponential",
    tests = c("SW", "SF", "LF", "JB"),
    Nsim = 1000,
    gen_data = two_sample_data,
    get_parameters = twosample_parameters,
    fn_to_get_norm_obj = raw_data,
    ...
) {
  
  #' Validate the simulation count.
  Nsim <- validate_positive_integer(Nsim, "Nsim")
  
  #' Standardize the threshold grid and test codes.
  alpha_pretest <- sort(unique(as.numeric(alpha_pretest)))
  tests <- unique(toupper(trimws(tests)))
  
  #' ---------------------------------------------------------------------------
  #' 1. Allocate p-value storage.
  #' ---------------------------------------------------------------------------
  #' P-values under Normality.
  pvalues_H0 <- matrix(NA_real_, nrow = length(tests), ncol = Nsim, 
                       dimnames = list(tests, NULL)
  )
  
  #' P-values under the alternative.
  pvalues_H1 <- pvalues_H0
  
  #' ---------------------------------------------------------------------------
  #' 2. Define the multi-sample p-value rule.
  #' ---------------------------------------------------------------------------
  minimum_p_value <- function(normality_object, test_name) {
    #' Apply the selected normality test.
    values <- tryCatch(normality_test(normality_object, test = test_name),
      error = function(error_condition) NA_real_
    )
    
    #' Retain finite p-values only.
    finite_values <- as.numeric(values)
    finite_values <- finite_values[is.finite(finite_values)]
    
    #' Use the minimum p-value across samples.
    if (length(finite_values) == 0L) {
      NA_real_
    } else {
      min(finite_values)
    }
  }
  
  #' Display simulation progress.
  cat("Computing ROC curves for classical normality tests\n")
  
  progress_bar <- utils::txtProgressBar(min = 0, max = Nsim, style = 3)
  
  on.exit(close(progress_bar),add = TRUE)
  
  #' ---------------------------------------------------------------------------
  #' 3. Simulate p-values under H0 and H1.
  #' ---------------------------------------------------------------------------
  for (simulation_index in seq_len(Nsim)) {
    #' Generate matched H0 and H1 parameter lists.
    parameters_H0 <- get_parameters(n = n, dist = "normal", par = NULL,...)
    parameters_H1 <- get_parameters(n = n, dist = H1_dist,...)
    #' Generate each dataset and extract its normality object.
    normality_object_H0 <- fn_to_get_norm_obj(
      do.call(gen_data, parameters_H0)
    )
    normality_object_H1 <- fn_to_get_norm_obj(do.call(gen_data, parameters_H1))
    #' Apply all requested tests to the same generated datasets.
    for (test_index in seq_along(tests)) {
  
      #' Current normality-test code.
      test_name <- tests[test_index]
      #' Store the H0 p-value.
      pvalues_H0[test_index, simulation_index] <- minimum_p_value(
        normality_object_H0, test_name
      )
      
      #' Store the H1 p-value.
      pvalues_H1[test_index, simulation_index] <- minimum_p_value(
        normality_object_H1, test_name
      )
    }
    
    #' Update progress.
    utils::setTxtProgressBar(progress_bar, simulation_index)
  }
  
  #' ---------------------------------------------------------------------------
  #' 4. Convert stored p-values to ROC coordinates.
  #' ---------------------------------------------------------------------------
  
  #' Allocate false-positive rates.
  false_positive_rate <- matrix(NA_real_, nrow = length(tests), ncol = length(alpha_pretest),
    dimnames = list( tests, paste0("alpha_", alpha_pretest))
  )
  #' Use the same dimensions for true-positive rates.
  true_positive_rate <- false_positive_rate
  #' Calculate rejection rates for each test and threshold.
  for (test_index in seq_along(tests)) {
    #' False-positive rates under Normality.
    false_positive_rate[test_index, ] <- vapply(
      alpha_pretest,
      function(alpha) {
        mean(pvalues_H0[test_index, ] <= alpha, na.rm = TRUE)
      },
      numeric(1)
    )
    
    #' True-positive rates under the alternative.
    true_positive_rate[test_index, ] <- vapply(
      alpha_pretest,
      function(alpha) {
        mean(pvalues_H1[test_index, ] <= alpha,na.rm = TRUE)
      },
      numeric(1)
    )
  }
  
  #' Return ROC coordinates and the underlying p-values.
  list(
    FPR = false_positive_rate,
    TPR = true_positive_rate,
    alpha = alpha_pretest,
    raw = list(
      H0 = pvalues_H0,
      H1 = pvalues_H1
    )
  )
}


#' -----------------------------------------------------------------------------
#' Plot one or more normality ROC curves and return their AUC values.
#'
#' Each element of curves must contain fpr, tpr, and label. An optional lty
#' value may be supplied to control the curve line type. The function removes
#' invalid coordinates, orders each curve, calculates its trapezoidal AUC, and
#' displays the AUC in the legend.
#' -----------------------------------------------------------------------------
plot_norm_roc_curve <- function(
    curves,
    title = NULL,
    dist_name = NULL,
    legend_title = "Normality Tests",
    legend_position = "bottomright",
    line_width = 3,
    point_size = 0.7
) {
  
  #' Require at least one curve specification.
  if (!is.list(curves) || length(curves) == 0L) {
    stop("curves must be a nonempty list of ROC-curve specifications.", call. = FALSE
    )
  }
  
  #' Construct an informative default title.
  if (is.null(title)) {
    title <- if (is.null(dist_name)) {
      "ROC Curves for Normality Tests"
    } else {
      paste("ROC Curves for Normality Tests | Alternative:", dist_name)
    }
  }
  
  #' Number of curves and plotting identifiers.
  n_curves <- length(curves)
  curve_colors <- seq_len(n_curves)
  curve_symbols <- seq(16, length.out = n_curves)
  
  #' Preserve the active graphics settings.
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  #' Set readable plot margins and text sizes.
  graphics::par(
    mar = c(5, 5, 4, 2) + 0.1,
    cex.lab = 1.2,
    cex.axis = 1.1,
    cex.main = 1.35,
    font.lab = 2,
    font.axis = 2
  )
  
  #' ---------------------------------------------------------------------------
  #' 1. Initialize the ROC plotting region.
  #' ---------------------------------------------------------------------------
  
  graphics::plot(
    NA,
    xlim = c(0, 1),
    ylim = c(0, 1),
    xlab = "False Positive Rate",
    ylab = "True Positive Rate",
    main = title
  )
  
  #' Add a reference grid.
  graphics::grid(col = "gray88", lty = "dotted")
  
  #' Add the random-classification reference line.
  graphics::abline(a = 0, b = 1, lty = 2, col = "gray60", lwd = 2)
  
  #' Allocate AUC storage.
  auc_values <- numeric(n_curves)
  
  #' ---------------------------------------------------------------------------
  #' 2. Draw each ROC curve.
  #' ---------------------------------------------------------------------------
  for (curve_index in seq_along(curves)) {
    #' Current curve specification.
    curve <- curves[[curve_index]]
    #' Retain complete ROC coordinates.
    valid_points <- is.finite(curve$fpr) & is.finite(curve$tpr)
    
    false_positive_rate <- as.numeric(curve$fpr[valid_points])
    true_positive_rate <- as.numeric(curve$tpr[valid_points])
    
    #' Order coordinates from left to right.
    sort_index <- order(false_positive_rate, true_positive_rate)
    
    false_positive_rate <- false_positive_rate[sort_index]
    true_positive_rate <- true_positive_rate[sort_index]
    
    #' Draw the ROC line.
    graphics::lines(
      false_positive_rate,
      true_positive_rate,
      col = curve_colors[curve_index],
      lwd = line_width,
      lty = curve$lty %||% 1
    )
    
    #' Mark evaluated thresholds.
    graphics::points(
      false_positive_rate,
      true_positive_rate,
      col = curve_colors[curve_index],
      pch = curve_symbols[curve_index],
      cex = point_size
    )
    
    #' Calculate the curve AUC.
    auc_values[curve_index] <- compute_auc(
      false_positive_rate,
      true_positive_rate,
      ensure_endpoints = TRUE
    )
  }
  
  #' ---------------------------------------------------------------------------
  #' 3. Construct the legend.
  #' ---------------------------------------------------------------------------
  
  #' Extract curve labels and line types.
  curve_labels <- vapply(curves, function(curve) curve$label,character(1))
  
  curve_line_types <- vapply(curves, function(curve) curve$lty %||% 1, numeric(1))
  
  #' Append AUC values to the labels.
  legend_text <- sprintf("%s (AUC = %.3f)", curve_labels, auc_values)
  
  #' Display the curve legend.
  graphics::legend(
    legend_position,
    legend = legend_text,
    col = curve_colors,
    pch = curve_symbols,
    lty = curve_line_types,
    lwd = line_width,
    pt.cex = 0.9,
    title = legend_title,
    cex = 0.9,
    bty = "o"
  )
  
  #' Return named AUC values without printing them.
  invisible(
    stats::setNames( auc_values,curve_labels)
  )
}


#' -----------------------------------------------------------------------------
#' Simulate and plot classical or custom normality ROC curves.
#'
#' Classical mode evaluates one or more normality tests in a common simulation
#' run. Custom mode evaluates the score function stored in norm_config$config$fn.
#' The first entry in distributions is treated as the non-Normal alternative.
#'
#' The function plots the selected curves and invisibly returns the complete ROC
#' object, including the stored H0 and H1 p-values or scores.
#' -----------------------------------------------------------------------------
plot_normality_roc_wrapper <- function(
    norm_config,
    single_n,
    threshold_grid,
    distributions,
    Nsim = 1000,
    norm_test = c("SW", "SF", "LF", "KS", "JB","SKEW", "KURT", "DAP", "AD", "CVM"),
    selected_tests = c("SW", "SF", "LF", "JB", "SKEW", "DAP", "AD", "CVM"),
    gen_data,
    get_parameters,
    fn_to_get_norm_obj,
    center_by = "median",
    ...
) {
  
  #' Standardize and validate the configured method.
  method <- match.arg(
    tolower(norm_config$method %||% "classical"),
    c("classical", "custom")
  )
  
  #' ---------------------------------------------------------------------------
  #' 1. Classical normality-test ROC analysis.
  #' ---------------------------------------------------------------------------
  if (method == "classical") {
    
    #' Read the primary configured test.
    configured_test <- norm_config$config$norm_test %||% "SW"
    #' Run the complete requested test set when applicable.
    tests_to_run <- if (configured_test %in% norm_test) {
      norm_test
    } else {
      configured_test
    }
    
    #' Simulate classical-test ROC coordinates.
    roc_object <- fn_for_norm_test_roc_curve(
      n = single_n,
      alpha_pretest = threshold_grid,
      H1_dist = distributions[1L],
      tests = tests_to_run,
      Nsim = Nsim,
      center_by = center_by,
      gen_data = gen_data,
      get_parameters = get_parameters,
      fn_to_get_norm_obj = fn_to_get_norm_obj,
      ...
    )
    
    #' Identify the tests produced by the simulation.
    computed_tests <- rownames(roc_object$FPR)
    #' Retain only the requested display curves.
    tests_to_plot <- intersect(selected_tests, computed_tests)
    
    #' Display all computed tests when no selection matches.
    if (length(tests_to_plot) == 0L) {
      tests_to_plot <- computed_tests
    }
    
    #' Convert matrix rows to plotting specifications.
    curves <- lapply(
      tests_to_plot,
      function(test_name) {
        list(
          fpr = as.numeric(roc_object$FPR[test_name, ]),
          tpr = as.numeric(roc_object$TPR[test_name, ]),
          label = test_name,
          lty = 1
        )
      }
    )
    
    #' Plot the selected classical curves.
    plot_norm_roc_curve(curves = curves, dist_name = distributions[1L])
    #' Return the complete simulation object.
    return(invisible(roc_object))
  }
  
  #' ---------------------------------------------------------------------------
  #' 2. Custom normality-score ROC analysis.
  #' ---------------------------------------------------------------------------
  
  #' Require a custom score function.
  if (is.null(norm_config$config$fn) ||!is.function(norm_config$config$fn)) {
    stop("norm_config$config$fn must be a function for method = 'custom'.", call. = FALSE)
  }
  
  #' Extract the custom scoring function.
  custom_function <- norm_config$config$fn
  
  #' Forward nonmetadata settings to the custom function.
  custom_arguments <- norm_config$config[setdiff(names(norm_config$config), c("fn", "label", "name"))]
  
  #' Simulate custom-score ROC coordinates.
  roc_object <- fn_for_custom_test_roc_curve(
    n = single_n,
    threshold_grid = threshold_grid,
    H1_dist = distributions[1L],
    Nsim = Nsim,
    gen_data = gen_data,
    get_parameters = get_parameters,
    fn_to_get_norm_obj = fn_to_get_norm_obj,
    custom_fn = custom_function,
    custom_args = custom_arguments,
    center_by = center_by,
    ...
  )
  
  #' Use the configured label for a single returned score.
  if (nrow(roc_object$FPR) == 1L) {
  
    #' Resolve the preferred display label.
    curve_label <- norm_config$config$label %||% norm_config$config$name %||% "Custom"

    #' Apply the same label to both ROC matrices.
    rownames(roc_object$FPR) <- curve_label
    rownames(roc_object$TPR) <- curve_label
  }
  
  #' Read the final custom curve names.
  curve_names <- rownames(roc_object$FPR)
  
  #' Convert matrix rows to plotting specifications.
  curves <- lapply(
    seq_along(curve_names),
    function(curve_index) {
      list(
        fpr = as.numeric(roc_object$FPR[curve_index, ]),
        tpr = as.numeric(roc_object$TPR[curve_index, ]),
        label = curve_names[curve_index],
        lty = 1
      )
    }
  )
  
  #' Plot the custom ROC curves.
  plot_norm_roc_curve(curves = curves,dist_name = distributions[1L]
  )
  
  #' Return the complete simulation object.
  invisible(roc_object)
}

#' =============================================================================
#' 11. CORE ADAPTIVE SIMULATION AND THRESHOLD SELECTION
#'
#' This section contains the main simulation engine for the adaptive procedure.
#' It generates downstream-test p-values and normality-assessment scores under
#' both hypotheses, then evaluates Type I error and power over a grid of
#' pretest thresholds without regenerating the simulated datasets.
#' =============================================================================

#' -----------------------------------------------------------------------------
#' Simulate downstream p-values and normality scores under H0 and H1.
#'
#' For each simulation replicate, the function:
#' 1. generates one dataset under H0 and one dataset under H1;
#' 2. applies both downstream tests to the same generated datasets;
#' 3. extracts the corresponding normality-assessment objects;
#' 4. calculates classical p-values or custom normality scores; and
#' 5. stores all results for later threshold evaluation.
#'
#' Failed downstream tests are stored as NA. Failed normality assessments are
#' also stored as NA so that the adaptive routing rule defaults to test 2.
#' -----------------------------------------------------------------------------
generate_pval <- function(
    Nsim,
    n,
    effect_size_H1,
    effect_size_H0 = zero_like(effect_size_H1),
    norm_config = list(method = "classical", config = list(norm_test = "SW")),
    dist = "normal",
    gen_data,
    get_parameters,
    fn_to_get_norm_obj,
    fn_for_ds_test_1,
    fn_for_ds_test_2,
    ...
) {
  
  #' Validate the simulation count.
  Nsim <- validate_positive_integer(Nsim, "Nsim")
  
  #' Standardize and validate the normality-assessment method.
  method <- match.arg(tolower(norm_config$method %||% ""),c("classical", "custom"))
  
  #' ---------------------------------------------------------------------------
  #' 1. Allocate simulation storage.
  #' ---------------------------------------------------------------------------
  #' Downstream test 1 p-values.
  pval_ds_test1_H0 <- rep(NA_real_, Nsim)
  pval_ds_test1_H1 <- rep(NA_real_, Nsim)
  
  #' Downstream test 2 p-values.
  pval_ds_test2_H0 <- rep(NA_real_, Nsim)
  pval_ds_test2_H1 <- rep(NA_real_, Nsim)
  
  #' Normality p-values or custom scores.
  norm_decision_H0 <- vector("list", Nsim)
  norm_decision_H1 <- vector("list", Nsim)
  
  #' ---------------------------------------------------------------------------
  #' 2. Define safe downstream p-value extraction.
  #' ---------------------------------------------------------------------------
  extract_p_value <- function(test_function, data_object) {
  
    #' Run the downstream test and trap failures.
    result <- tryCatch(test_function(data_object),
      error = function(error_condition) NULL
    )
    
    #' Return NA when no p-value is available.
    if (is.null(result) || is.null(result$p.value)) {
      return(NA_real_)
    }
    
    #' Retain the first finite p-value.
    p_value <- as.numeric(result$p.value)[1L]
    if (is.finite(p_value)) {
      p_value
    } else {
      NA_real_
    }
  }
  
  #' Display simulation progress.
  progress_bar <- utils::txtProgressBar(min = 0, max = Nsim, style = 3)
  
  on.exit(close(progress_bar), add = TRUE)
  
  #' ---------------------------------------------------------------------------
  #' 3. Simulate H0 and H1 results.
  #' ---------------------------------------------------------------------------
  for (simulation_index in seq_len(Nsim)) {
    
    #' Construct the H0 generator arguments.
    parameters_H0 <- get_parameters(n = n, effect_size = effect_size_H0,dist = dist,...)
    
    #' Construct the H1 generator arguments.
    parameters_H1 <- get_parameters(n = n, effect_size = effect_size_H1, dist = dist, ...)
    
    #' Generate one dataset under each hypothesis.
    data_H0 <- do.call(gen_data, parameters_H0)
    
    data_H1 <- do.call(gen_data, parameters_H1)
    
    #' Apply downstream test 1 to both datasets.
    pval_ds_test1_H0[simulation_index] <- extract_p_value(fn_for_ds_test_1, data_H0)
    
    pval_ds_test1_H1[simulation_index] <- extract_p_value(fn_for_ds_test_1, data_H1)
    
    #' Apply downstream test 2 to both datasets.
    pval_ds_test2_H0[simulation_index] <- extract_p_value(fn_for_ds_test_2, data_H0)
    
    pval_ds_test2_H1[simulation_index] <- extract_p_value(fn_for_ds_test_2, data_H1)
    
    #' Extract the objects used for normality assessment.
    normality_object_H0 <- fn_to_get_norm_obj(data_H0)
    normality_object_H1 <- fn_to_get_norm_obj(data_H1)
    
    #' Assess Normality under H0.
    norm_decision_H0[[simulation_index]] <- tryCatch(
      assess_normality(normality_object_H0, method = method, config = norm_config$config),
      error = function(error_condition) {
        c(score = NA_real_)
      }
    )
    
    #' Assess Normality under H1.
    norm_decision_H1[[simulation_index]] <- tryCatch(
      assess_normality(normality_object_H1, method = method,config = norm_config$config),
      error = function(error_condition) {
        c(score = NA_real_)
      }
    )
    
    #' Update progress.
    utils::setTxtProgressBar(progress_bar, simulation_index)
  }
  
  #' ---------------------------------------------------------------------------
  #' 4. Return all stored simulation quantities.
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


#' -----------------------------------------------------------------------------
#' Evaluate adaptive Type I error and power over a pretest-threshold grid.
#'
#' The function first calls generate_pval() once for each distribution. It then
#' reuses the stored downstream p-values and normality scores to evaluate every
#' threshold in threshold_grid.
#'
#' The first distribution is treated as the non-Normal alternative used for
#' Type I error assessment. The second distribution must be Normal and is used
#' to assess expected power behavior under the model assumptions.
#'
#' At each threshold, test 1 is selected when every required normality score
#' exceeds the threshold. Otherwise, the adaptive procedure uses test 2.
#' -----------------------------------------------------------------------------
perform_adaptive_analysis <- function(
    Nsim = 1000,
    n = 10,
    effect_size_H1 = 0.5,
    effect_size_H0 = zero_like(effect_size_H1),
    distributions = c("exponential", "normal"),
    norm_config = list(method = "classical", config = list(norm_test = "SW")),
    threshold_grid = seq(0, 1, by = 0.025),
    test_alpha = 0.05,
    gen_data = onesample_data,
    get_parameters = onesample_parameters,
    fn_to_get_norm_obj = raw_data,
    fn_for_ds_test_1 = one_sample_t_test,
    fn_for_ds_test_2 = sign_test,
    ...
) {
  
  #' Validate the simulation count.
  Nsim <- validate_positive_integer(Nsim, "Nsim")
  #' Standardize the threshold grid.
  threshold_grid <- sort(unique(as.numeric(threshold_grid)))
  
  #' Require at least two valid thresholds.
  if (length(threshold_grid) < 2L || any(!is.finite(threshold_grid))) {
    stop("threshold_grid must contain at least two finite values.", call. = FALSE)
  }
  
  #' Require a non-Normal distribution followed by Normal.
  if (length(distributions) < 2L) {
    stop("distributions must contain a non-normal distribution followed by normal.", call. = FALSE)
  }
  
  #' Confirm the Normal reference distribution.
  if (tolower(distributions[2L]) != "normal") {
    stop("distributions[2] must be 'normal'.", call. = FALSE)
  }
  
  #' ---------------------------------------------------------------------------
  #' 1. Allocate result containers.
  #' ---------------------------------------------------------------------------
  #' Raw simulation output for each distribution.
  all_pvalues <- list()
  #' Type I error summaries.
  typeI_rates <- list()
  #' Power summaries.
  power_rates <- list()
  
  #' ---------------------------------------------------------------------------
  #' 2. Define the normality-score routing rule.
  #' ---------------------------------------------------------------------------
  routing_score <- function(decision_values) {

    #' Retain finite sample-specific scores.
    finite_values <- as.numeric(decision_values)
    finite_values <- finite_values[is.finite(finite_values)]
    
    #' Use the minimum score across samples.
    if (length(finite_values) == 0L) {
      NA_real_
    } else {
      min(finite_values)
    }
  }
  
  #' Display distribution-level progress.
  progress_bar <- utils::txtProgressBar(min = 0, max = length(distributions), style = 3)
  
  on.exit(close(progress_bar),add = TRUE)
  
  #' ---------------------------------------------------------------------------
  #' 3. Generate simulation results for each distribution.
  #' ---------------------------------------------------------------------------
  
  for (distribution_index in seq_along(distributions)) {
    
    #' Current data-generating distribution.
    distribution <- distributions[distribution_index]
    
    #' Simulate all p-values and normality scores once.
    simulation_results <- generate_pval(
      Nsim = Nsim,
      n = n,
      effect_size_H1 = effect_size_H1,
      effect_size_H0 = effect_size_H0,
      norm_config = norm_config,
      dist = distribution,
      gen_data = gen_data,
      get_parameters = get_parameters,
      fn_to_get_norm_obj = fn_to_get_norm_obj,
      fn_for_ds_test_1 = fn_for_ds_test_1,
      fn_for_ds_test_2 = fn_for_ds_test_2,
      ...
    )
    
    #' Store the complete simulation output.
    all_pvalues[[distribution]] <- simulation_results
    
    #' Extract H0 downstream p-values.
    pvalue_H0_test1 <- simulation_results$pval_ds_test1_H0
    pvalue_H0_test2 <- simulation_results$pval_ds_test2_H0
    
    #' Extract H1 downstream p-values.
    pvalue_H1_test1 <- simulation_results$pval_ds_test1_H1
    pvalue_H1_test2 <- simulation_results$pval_ds_test2_H1
    
    #' Reduce H0 normality outputs to one routing score per replicate.
    score_H0 <- vapply(
      simulation_results$norm_decision_H0,
      routing_score,
      numeric(1)
    )
    
    #' Reduce H1 normality outputs to one routing score per replicate.
    score_H1 <- vapply(
      simulation_results$norm_decision_H1,
      routing_score,
      numeric(1)
    )
    
    #' -------------------------------------------------------------------------
    #' 4. Calculate fixed-test operating characteristics.
    #' -------------------------------------------------------------------------
    
    #' Type I error for each fixed downstream test.
    typeI_rates[[distribution]] <- list(
      test_1 = mean(pvalue_H0_test1 <= test_alpha, na.rm = TRUE),
      test_2 = mean(pvalue_H0_test2 <= test_alpha, na.rm = TRUE),
      adaptive = numeric(length(threshold_grid))
    )
    
    #' Power for each fixed downstream test.
    power_rates[[distribution]] <- list(
      test_1 = mean(pvalue_H1_test1 <= test_alpha, na.rm = TRUE),
      test_2 = mean(pvalue_H1_test2 <= test_alpha,na.rm = TRUE),
      adaptive = numeric(length(threshold_grid))
    )
    
    #' -------------------------------------------------------------------------
    #' 5. Evaluate every adaptive threshold.
    #' -------------------------------------------------------------------------
    
    for (threshold_index in seq_along(threshold_grid)) {
      
      #' Current routing threshold.
      threshold <- threshold_grid[threshold_index]
      
      #' Select test 1 when the H0 score exceeds the threshold.
      use_test1_H0 <- is.finite(score_H0) & score_H0 > threshold
      
      #' Select test 1 when the H1 score exceeds the threshold.
      use_test1_H1 <- is.finite(score_H1) & score_H1 > threshold
      
      #' Route H0 p-values between the two downstream tests.
      adaptive_H0 <- ifelse(use_test1_H0, pvalue_H0_test1, pvalue_H0_test2)
      
      #' Route H1 p-values between the two downstream tests.
      adaptive_H1 <- ifelse(use_test1_H1, pvalue_H1_test1, pvalue_H1_test2)
      
      #' Adaptive Type I error.
      typeI_rates[[distribution]]$adaptive[threshold_index] <- mean(
        adaptive_H0 <= test_alpha,
        na.rm = TRUE
      )
      
      #' Adaptive power.
      power_rates[[distribution]]$adaptive[threshold_index] <- mean(
        adaptive_H1 <= test_alpha,
        na.rm = TRUE
      )
    }
    
    #' Update distribution-level progress.
    utils::setTxtProgressBar(progress_bar, distribution_index)
  }
  
  #' ---------------------------------------------------------------------------
  #' 6. Return operating characteristics and reproducibility settings.
  #' ---------------------------------------------------------------------------
  
  list(
    typeI_rates = typeI_rates,
    power_rates = power_rates,
    all_pvalues = all_pvalues,
    param_grid = threshold_grid,
    param_name = "threshold",
    settings = list(
      Nsim = Nsim,
      n = n,
      effect_size_H1 = effect_size_H1,
      effect_size_H0 = effect_size_H0,
      distributions = distributions,
      norm_config = norm_config,
      threshold_grid = threshold_grid,
      test_alpha = test_alpha,
      gen_data = deparse(substitute(gen_data)),
      get_parameters = deparse(substitute(get_parameters)),
      fn_to_get_norm_obj = deparse(substitute(fn_to_get_norm_obj)),
      fn_for_ds_test_1 = deparse(substitute(fn_for_ds_test_1)),
      fn_for_ds_test_2 = deparse(substitute(fn_for_ds_test_2))
    )
  )
}

#' -----------------------------------------------------------------------------
#' Compute power and Type I error metrics for the adaptive procedure.
#'
#' The function compares the adaptive procedure with the two fixed downstream
#' tests under:
#' 1. the non-Normal distribution stored in the first list element; and
#' 2. the Normal distribution stored in the second list element.
#'
#' Expected power loss measures the power sacrificed relative to test 1 under
#' Normality. Expected power gain measures the improvement relative to test 1
#' under non-Normality. Type I error inflation measures the difference between
#' the adaptive Type I error rate and the nominal test level.
#' -----------------------------------------------------------------------------
compute_roc_metrics <- function(typeI_rates, power_rates, test_alpha) {
  
  #' ---------------------------------------------------------------------------
  #' 1. Validate the distribution ordering.
  #' ---------------------------------------------------------------------------
  
  #' Read the distribution names from the power results.
  distribution_names <- names(power_rates)
  
  #' Confirm that the Normal distribution is stored second.
  if (!is.null(distribution_names) && length(distribution_names) >= 2L) {
    
    if (tolower(distribution_names[2L]) != "normal") {
      stop("Distributions must be ordered as non-normal followed by normal; ",
        "the second element is '", distribution_names[2L], "'.", call. = FALSE
      )
    }
    
    #' Warn when the ordering cannot be verified.
  } else {
    warning("Cannot verify distribution order: power_rates is not a named list.",
      call. = FALSE)
  }
  
  #' ---------------------------------------------------------------------------
  #' 2. Extract operating characteristics by distribution.
  #' ---------------------------------------------------------------------------
  
  #' First distribution: non-Normal setting.
  non_normal <- list(
    power_test1 = power_rates[[1L]]$test_1,
    power_test2 = power_rates[[1L]]$test_2,
    power_adaptive = power_rates[[1L]]$adaptive,
    error_adaptive = typeI_rates[[1L]]$adaptive
  )
  
  #' Second distribution: Normal setting.
  normal <- list(
    power_test1 = power_rates[[2L]]$test_1,
    power_test2 = power_rates[[2L]]$test_2,
    power_adaptive = power_rates[[2L]]$adaptive,
    error_adaptive = typeI_rates[[2L]]$adaptive
  )
  
  #' ---------------------------------------------------------------------------
  #' 3. Compute adaptive-procedure metrics.
  #' ---------------------------------------------------------------------------
  
  #' Expected power loss under Normality.
  EPL <- normal$power_test1 - normal$power_adaptive
  
  #' Expected power gain under non-Normality.
  EPG <- non_normal$power_adaptive - non_normal$power_test1
  
  #' Type I error inflation under Normality.
  ETIE_normal <- normal$error_adaptive - test_alpha
  
  #' Type I error inflation under non-Normality.
  ETIE_non_normal <- non_normal$error_adaptive - test_alpha
  
  #' ---------------------------------------------------------------------------
  #' 4. Compute fixed-test benchmark differences.
  #' ---------------------------------------------------------------------------
  #' Power gain from always using test 2 under non-Normality.
  benchmark_power_gain <- non_normal$power_test2 - non_normal$power_test1
  
  #' Power loss from always using test 2 under Normality.
  benchmark_power_loss <- normal$power_test1 - normal$power_test2
  
  #' ---------------------------------------------------------------------------
  #' 5. Return all metric vectors.
  #' ---------------------------------------------------------------------------
  
  list(
    Expected_power_loss = EPL,
    Expected_power_gain = EPG,
    Expected_type1_error_inflation_normal = ETIE_normal,
    Expected_type1_error_inflation_non_normal = ETIE_non_normal,
    benchmark_power_gain = benchmark_power_gain,
    benchmark_power_loss = benchmark_power_loss
  )
}

#' -----------------------------------------------------------------------------
#' Select the threshold that prioritizes Type I error control and power.
#'
#' The procedure uses a two-stage rule:
#' 1. retain thresholds whose worst-case positive Type I error inflation is no
#'    greater than tol_pos, then select among them using power performance; or
#' 2. when no threshold satisfies the Type I error tolerance, minimize the
#'    worst-case positive inflation first and then apply the power criterion.
#'
#' Within a candidate set, thresholds with power loss no greater than loss_tol
#' are preferred. Among those thresholds, the largest non-Normal power gain is
#' selected. Remaining ties are resolved by choosing the smallest threshold.
#'
#' Nonfinite metric values are removed before selection. The function also
#' records whether threshold-grid endpoints 0 and 1 are present because those
#' endpoints help diagnose why the controlled candidate set may be empty.
#' -----------------------------------------------------------------------------
select_optimal_parameter <- function(
    distributions = c("exponential", "normal"),
    param_grid,
    metrics,
    alpha = 0.05,
    tol_pos = 0.01,
    loss_tol = 0.01,
    param_name = "threshold",
    effect_size = NULL
) {
  
  #' ---------------------------------------------------------------------------
  #' 1. Extract and validate metric vectors.
  #' ---------------------------------------------------------------------------

  #' Convert stored metrics to plain numeric vectors.
  power_loss_normal <- as.numeric(metrics$Expected_power_loss)
  power_gain_non_normal <- as.numeric(metrics$Expected_power_gain)
  type1_infl_normal <- as.numeric(metrics$Expected_type1_error_inflation_normal)
  type1_infl_non_normal <- as.numeric(metrics$Expected_type1_error_inflation_non_normal)
  
  #' Record the length of every threshold-dependent quantity.
  input_lengths <- c(
    param_grid = length(param_grid),
    power_loss_normal = length(power_loss_normal),
    power_gain_non_normal = length(power_gain_non_normal),
    type1_infl_normal = length(type1_infl_normal),
    type1_infl_non_normal = length(type1_infl_non_normal)
  )
  
  #' Require one metric value for every threshold.
  if (length(unique(input_lengths)) != 1L) {
    stop("param_grid and all metric vectors must have identical lengths. ",
      paste(names(input_lengths),input_lengths, sep = " = ", collapse = ", "),
      call. = FALSE
    )
  }
  
  #' Preserve positions from the original threshold grid.
  original_index <- seq_along(param_grid)
  
  #' ---------------------------------------------------------------------------
  #' 2. Remove thresholds with undefined metrics.
  #' ---------------------------------------------------------------------------
  
  #' Retain only thresholds with finite parameter and metric values.
  valid <- is.finite(param_grid) & is.finite(power_loss_normal) &
    is.finite(power_gain_non_normal) & is.finite(type1_infl_normal) &
    is.finite(type1_infl_non_normal)
  
  #' Stop when no candidate can be evaluated.
  if (!any(valid)) {
    stop("No valid candidate thresholds remain after removing non-finite values.",
      call. = FALSE)
  }
  
  #' Record the number of excluded candidates.
  n_removed_nonfinite <- sum(!valid)
  
  #' Restrict all quantities to valid thresholds.
  param_values <- param_grid[valid]
  original_index <- original_index[valid]
  power_loss <- power_loss_normal[valid]
  power_gain <- power_gain_non_normal[valid]
  type1_infl_normal <- type1_infl_normal[valid]
  type1_infl_non_normal <- type1_infl_non_normal[valid]
  
  #' ---------------------------------------------------------------------------
  #' 3. Summarize positive Type I error inflation.
  #' ---------------------------------------------------------------------------
  
  #' Numerical tolerance for floating-point comparisons.
  eps <- .Machine$double.eps^0.5
  #' Ignore conservative deviations below the nominal level.
  pos_infl_normal <- pmax(type1_infl_normal, 0)
  pos_infl_non_normal <- pmax(type1_infl_non_normal, 0)
  
  #' Use the larger positive inflation across the two scenarios.
  worst_pos_infl <- pmax(pos_infl_normal, pos_infl_non_normal)
  
  #' Check whether the routing endpoints are represented.
  has_zero_endpoint <- any(abs(param_values - 0) <= eps)
  has_one_endpoint <- any(abs(param_values - 1) <= eps)
  
  has_both_endpoints <- has_zero_endpoint && has_one_endpoint
  
  #' ---------------------------------------------------------------------------
  #' 4. Define the power-based rule within a candidate set.
  #' ---------------------------------------------------------------------------
  
  select_best <- function(candidates) {
    #' Require at least one candidate threshold.
    if (length(candidates) == 0L) {
      stop("select_best() received an empty candidate set.", call. = FALSE)
    }
    
    #' Extract power losses for the candidate thresholds.
    loss_values <- power_loss[candidates]
    
    #' Prefer candidates satisfying the permitted power-loss tolerance.
    good_loss <- candidates[loss_values <= loss_tol + eps]
    
    if (length(good_loss) > 0L) {
      #' Maximize non-Normal power gain among low-loss candidates.
      max_gain <- max(power_gain[good_loss],na.rm = TRUE)
      best <- good_loss[abs(power_gain[good_loss] - max_gain) <= eps * (1 + abs(max_gain))]
      loss_controlled <- TRUE
      
    } else {
      
      #' No candidate satisfies loss_tol, so minimize power loss first.
      min_loss <- min(loss_values, na.rm = TRUE)
      min_loss_group <- candidates[abs(loss_values - min_loss) <= eps * (1 + abs(min_loss))]
      
      #' Maximize power gain among the minimum-loss candidates.
      max_gain <- max(power_gain[min_loss_group], na.rm = TRUE)
      best <- min_loss_group[abs(power_gain[min_loss_group] - max_gain) <= eps * (1 + abs(max_gain))]
      
      loss_controlled <- FALSE
    }
    
    #' Resolve remaining ties by selecting the smallest threshold.
    selected <- if (length(best) > 1L) {
      best[which.min(param_values[best])]
    } else {
      best
    }
    
    #' Retain whether the selected threshold met loss_tol.
    attr(selected, "loss_controlled") <- loss_controlled
    
    selected
  }
  
  #' ---------------------------------------------------------------------------
  #' 5. Apply the two-stage threshold-selection rule.
  #' ---------------------------------------------------------------------------
  
  #' Stage 1 candidates satisfy the Type I error inflation tolerance.
  controlled_index <- which(worst_pos_infl <= tol_pos + eps)
  
  if (length(controlled_index) > 0L) {
    
    #' Select by power performance within the controlled set.
    candidate_index <- controlled_index
    optimal_index <- select_best(candidate_index)
    
    #' Record successful Type I error control.
    selection_stage <- 1L
    inflation_controlled <- TRUE
    
    selection_note <- sprintf(
      paste0("Stage 1: worst-case positive Type I error inflation ", "within tolerance (<= %.3g)"
      ),
      tol_pos
    )
    
  } else {
    
    #' Stage 2 minimizes worst-case inflation when no threshold is controlled.
    minimum_inflation <- min(worst_pos_infl, na.rm = TRUE)
    candidate_index <- which(abs(worst_pos_infl - minimum_inflation) <= eps * (1 + abs(minimum_inflation)))
    
    #' Apply the same power rule within the minimum-inflation set.
    optimal_index <- select_best(candidate_index)
    
    #' Record use of the fallback rule.
    selection_stage <- 2L
    inflation_controlled <- FALSE
    
    #' Explain whether missing grid endpoints may have contributed.
    selection_note <- if (has_both_endpoints) {
      sprintf(paste0("Stage 2 fallback: endpoints were present but no threshold met ",
          "tol_pos; minimized worst-case positive inflation (min = %.3g)"
          ), 
          minimum_inflation)
    } else {
      sprintf(
        paste0("Stage 2 fallback: no threshold met tol_pos, possibly because the ",
          "grid omitted 0 or 1; minimized worst-case positive inflation ", "(min = %.3g)"
        ),
        minimum_inflation
      )
    }
  }
  
  #' Read the power-loss status stored by select_best().
  loss_controlled <- isTRUE(attr(optimal_index, "loss_controlled")
  )
  
  #' ---------------------------------------------------------------------------
  #' 6. Return the selected threshold and complete diagnostics.
  #' ---------------------------------------------------------------------------
  
  list(
    param_star = param_values[optimal_index],
    filtered_index = optimal_index,
    optimal_index = optimal_index,
    original_index = original_index[optimal_index],
    power_loss = power_loss[optimal_index],
    power_gain = power_gain[optimal_index],
    inflation_normal = type1_infl_normal[optimal_index],
    inflation_non_normal = type1_infl_non_normal[optimal_index],
    worst_case_inflation = worst_pos_infl[optimal_index],
    inflation_controlled = inflation_controlled,
    loss_controlled = loss_controlled,
    feasible = inflation_controlled,
    selection_stage = selection_stage,
    selection_note = selection_note,
    n_removed_nonfinite = n_removed_nonfinite,
    threshold_grid_has_zero_endpoint = has_zero_endpoint,
    threshold_grid_has_one_endpoint = has_one_endpoint,
    threshold_grid_has_both_endpoints = has_both_endpoints,
    param_name = param_name,
    effect_size = effect_size,
    
    #' Retain all filtered values used in the selection.
    all_metrics = list(
      param_values = param_values,
      original_index = original_index,
      power_loss = power_loss,
      power_gain = power_gain,
      type1_infl_normal = type1_infl_normal,
      type1_infl_non_normal = type1_infl_non_normal,
      pos_infl_normal = pos_infl_normal,
      pos_infl_non_normal = pos_infl_non_normal,
      worst_case_positive_inflation = worst_pos_infl,
      controlled_idx = controlled_index,
      candidate_idx = candidate_index,
      distributions = distributions,
      alpha = alpha,
      tol_pos = tol_pos,
      loss_tol = loss_tol,
      effect_size = effect_size,
      threshold_grid_has_zero_endpoint = has_zero_endpoint,
      threshold_grid_has_one_endpoint = has_one_endpoint,
      threshold_grid_has_both_endpoints = has_both_endpoints,
      param_name = param_name
    )
  )
}


#' -----------------------------------------------------------------------------
#' Plot the four threshold trade-off metrics and the selected threshold.
#'
#' The function creates a four-panel display for:
#' 1. expected power loss under Normality;
#' 2. expected power gain under non-Normality;
#' 3. Type I error inflation under Normality; and
#' 4. Type I error inflation under non-Normality.
#'
#' Each panel uses a data-driven vertical scale. A dashed vertical line marks the
#' selected threshold, and horizontal zero-reference lines are added to the Type I
#' error inflation panels. The function also prints the selected threshold and
#' its operating characteristics in an outer summary line.
#' -----------------------------------------------------------------------------
plot_tradeoff_results_generic <- function(
    optimal_result,
    outer_title = "Power and Type I Error Trade-off for Normality Pretesting",
    text_size_main = 1.0,
    text_size_labels = 1.0,
    line_width = 3
) {
  
  #' ---------------------------------------------------------------------------
  #' 1. Extract and validate the stored metrics.
  #' ---------------------------------------------------------------------------
  
  #' Retrieve all threshold-dependent quantities.
  metrics <- optimal_result$all_metrics
  
  #' Require the complete metric object produced by select_optimal_parameter().
  if (is.null(metrics)) {
    stop("optimal_result must contain all_metrics.", call. = FALSE)
  }
  
  #' Read the threshold grid and labels.
  param_values <- metrics$param_values
  param_name <- metrics$param_name %||% "threshold"
  distributions <- metrics$distributions %||% c("non-normal", "normal")
  
  #' ---------------------------------------------------------------------------
  #' 2. Define the four plot panels.
  #' ---------------------------------------------------------------------------
  
  #' Store panel-specific values, labels, colors, and reference lines.
  panel_definitions <- list(
    list(
        values = metrics$power_loss, 
         color = "blue", 
         ylab = "Expected Power Loss",
        title = sprintf("Expected Power Loss: %s", distributions[2L]),
      reference = NULL
    ),
    list(
      values = metrics$power_gain,
      color = "red",
      ylab = "Expected Power Gain",
      title = sprintf("Expected Power Gain: %s", distributions[1L]),
      reference = NULL
    ),
    list(
      values = metrics$type1_infl_normal,
      color = "orange",
      ylab = "Type I Error Inflation",
      title = sprintf("Type I Error Inflation: %s", distributions[2L]),
      reference = 0
    ),
    list(
      values = metrics$type1_infl_non_normal,
      color = "green4",
      ylab = "Type I Error Inflation",
      title = sprintf("Type I Error Inflation: %s", distributions[1L]),
      reference = 0
    )
  )
  
  #' Resolve a readable horizontal-axis label.
  xlab <- switch(
    param_name,
    threshold = "Pretest Threshold",
    pretest_threshold = expression(alpha[pre]),
    decision_threshold = "Decision Threshold",
    param_name
  )
  
  #' ---------------------------------------------------------------------------
  #' 3. Configure the four-panel graphics layout.
  #' ---------------------------------------------------------------------------
  
  #' Preserve the active graphics settings.
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  
  #' Use a two-by-two panel layout with space for outer annotations.
  graphics::par(
    mfrow = c(2, 2),
    oma = c(0.8, 0.2, 3.8, 0.2),
    mar = c(4, 3.8, 1.8, 0.8),
    mgp = c(2.2, 0.7, 0)
  )
  
  #' ---------------------------------------------------------------------------
  #' 4. Draw each trade-off metric.
  #' ---------------------------------------------------------------------------
  
  for (panel in panel_definitions) {
    
    #' Convert the current metric to a plain numeric vector.
    values <- as.numeric(panel$values)
    #' Determine the initial vertical range.
    y_limits <- range(values, na.rm = TRUE)
    
    #' Use a fallback range when all values are nonfinite.
    if (!all(is.finite(y_limits))) {
      y_limits <- c(-0.01, 0.01)
      
      #' Add visible padding when all finite values are identical.
    } else if (diff(y_limits) == 0) {
      padding <- max(abs(y_limits[1L]) * 0.05, 0.01)
      
      y_limits <- y_limits + c(-padding, padding)
      
      #' Add proportional padding for a nonconstant series.
    } else {
      padding <- 0.04 * diff(y_limits)
      y_limits <- y_limits + c(-padding, padding)
    }
    
    #' Plot the metric over the threshold grid.
    graphics::plot(
      param_values,
      values,
      type = "l",
      col = panel$color,
      lwd = line_width,
      ylim = y_limits,
      xlab = xlab,
      ylab = panel$ylab,
      main = panel$title,
      cex.lab = text_size_labels,
      cex.main = text_size_main,
      font.lab = 2
    )
    
    #' Add a light reference grid.
    graphics::grid(col = "gray88", lty = "dotted")
    
    #' Mark the selected threshold.
    graphics::abline(v = optimal_result$param_star, lty = 2, col = "darkred", lwd = 1.5)
    
    #' Add the zero-inflation reference line when applicable.
    if (!is.null(panel$reference)) {
      graphics::abline(h = panel$reference, lty = 3, col = "gray60")
    }
  }
  
  #' ---------------------------------------------------------------------------
  #' 5. Add the overall title and threshold summary.
  #' ---------------------------------------------------------------------------
  #' Add the common title above all four panels.
  graphics::mtext(outer_title, side = 3, outer = TRUE, line = 1.8, cex = 1.45, font = 2)
  
  #' Summarize the selected threshold and its main operating characteristics.
  summary_line <- sprintf(
    paste0("Selected threshold = %.4f | EPG = %.4f | EPL = %.4f | ",
      "Inflation = (%.4f, %.4f) | tol_pos = %.4f | loss_tol = %.4f"
    ),
    optimal_result$param_star,
    optimal_result$power_gain,
    optimal_result$power_loss,
    optimal_result$inflation_normal,
    optimal_result$inflation_non_normal,
    metrics$tol_pos,
    metrics$loss_tol
  )
  
  #' Display the summary beneath the overall title.
  graphics::mtext(summary_line, side = 3, outer = TRUE, line = 0.35, cex = 0.90, col = "#7aa2ff")
  
  #' Return the selection result without printing it.
  invisible(optimal_result)
}


#' -----------------------------------------------------------------------------
#' Construct power-versus-Type-I-error curves from stored simulation p-values.
#'
#' For each distribution, the function simulates downstream-test p-values once,
#' applies the fixed pretest threshold to route each replicate to test 1 or
#' test 2, and then evaluates power and Type I error over sig_levels.
#'
#' The resulting data frame contains one row for each distribution, method, and
#' significance level. The simulation settings are attached as an attribute for
#' plotting, verification, and reproducibility.
#' -----------------------------------------------------------------------------
power_vs_error_roc_data <- function(
    N = 1e3,
    n = 10,
    distributions = c("exponential", "normal"),
    norm_config = list(method = "classical", config = list(norm_test = "SW")),
    effect_size_H1,
    effect_size_H0 = zero_like(effect_size_H1),
    pretest_threshold = 0.05,
    optimal_result = NULL,
    sig_levels = seq(0, 1, by = 0.05),
    gen_data = anova_gen_data,
    get_parameters = anova_parameters,
    fn_to_get_norm_obj = anova_residuals,
    fn_for_ds_test_1 = one_way_anova,
    fn_for_ds_test_2 = kruskal_wallis_test,
    ...
) {
  
  #' ---------------------------------------------------------------------------
  #' 1. Validate the configuration and threshold.
  #' ---------------------------------------------------------------------------
  #' Require an explicit normality-assessment method.
  if (is.null(norm_config$method)) {
    stop("norm_config must have a $method slot: 'classical' or 'custom'.", call. = FALSE)
  }
  
  #' Use the selected optimal threshold when one is supplied.
  if (!is.null(optimal_result) && !is.null(optimal_result$param_star)) {
    pretest_threshold <- optimal_result$param_star
  }
  
  #' ---------------------------------------------------------------------------
  #' 2. Initialize result storage and progress reporting.
  #' ---------------------------------------------------------------------------
  #' Store all curve coordinates in long format.
  roc_results <- data.frame()
  #' Count all distribution-by-alpha evaluations.
  total_operations <- length(distributions) * length(sig_levels)
  
  progress_bar <- utils::txtProgressBar(min = 0, max = total_operations, style = 3)
  
  on.exit(close(progress_bar), add = TRUE)
  
  operation_counter <- 0L
  
  #' ---------------------------------------------------------------------------
  #' 3. Simulate and evaluate each distribution.
  #' ---------------------------------------------------------------------------
  
  for (distribution in distributions) {
    #' Report the current distribution.
    cat("Processing distribution:", distribution, "\n")
    
    #' Generate downstream p-values and normality scores once.
    simulation_results <- generate_pval(
      Nsim = N,
      n = n,
      effect_size_H1 = effect_size_H1,
      effect_size_H0 = effect_size_H0,
      norm_config = norm_config,
      dist = distribution,
      gen_data = gen_data,
      get_parameters = get_parameters,
      fn_to_get_norm_obj = fn_to_get_norm_obj,
      fn_for_ds_test_1 = fn_for_ds_test_1,
      fn_for_ds_test_2 = fn_for_ds_test_2,
      ...
    )
    
    #' -------------------------------------------------------------------------
    #' 4. Apply the fixed pretest routing rule.
    #' -------------------------------------------------------------------------
    #' Route H0 replicates to test 1 when Normality is supported.
    use_test1_H0 <- vapply(
      simulation_results$norm_decision_H0,
      use_parametric_test,
      logical(1),
      threshold = pretest_threshold
    )
    
    #' Route H1 replicates using the same rule.
    use_test1_H1 <- vapply(
      simulation_results$norm_decision_H1,
      use_parametric_test,
      logical(1),
      threshold = pretest_threshold
    )
    
    #' Construct adaptive H0 p-values.
    adaptive_pvalues_H0 <- ifelse(
      use_test1_H0,
      simulation_results$pval_ds_test1_H0,
      simulation_results$pval_ds_test2_H0
    )
    
    #' Construct adaptive H1 p-values.
    adaptive_pvalues_H1 <- ifelse(
      use_test1_H1,
      simulation_results$pval_ds_test1_H1,
      simulation_results$pval_ds_test2_H1
    )
    
    #' -------------------------------------------------------------------------
    #' 5. Sweep the downstream significance level.
    #' -------------------------------------------------------------------------
    
    for (alpha in sig_levels) {
      
      #' Add operating characteristics for both fixed tests and the adaptive test.
      roc_results <- rbind(
        roc_results,
        
        data.frame(
                Distribution = distribution,
                Method = "test_1",
                Alpha = alpha,
                Power = mean(simulation_results$pval_ds_test1_H1 <= alpha,na.rm = TRUE),
          TypeI_error = mean(simulation_results$pval_ds_test1_H0 <= alpha, na.rm = TRUE),
          stringsAsFactors = FALSE),
        
        data.frame(
          Distribution = distribution,
          Method = "test_2",
          Alpha = alpha,
          Power = mean(simulation_results$pval_ds_test2_H1 <= alpha, na.rm = TRUE),
          TypeI_error = mean(simulation_results$pval_ds_test2_H0 <= alpha,na.rm = TRUE),
          stringsAsFactors = FALSE),
        
        data.frame(
          Distribution = distribution,
          Method = "adaptive",
          Alpha = alpha,
          Power = mean(adaptive_pvalues_H1 <= alpha, na.rm = TRUE),
          TypeI_error = mean(adaptive_pvalues_H0 <= alpha, na.rm = TRUE),
          stringsAsFactors = FALSE)
      )
      
      #' Update progress.
      operation_counter <- operation_counter + 1L
      
      utils::setTxtProgressBar(progress_bar, operation_counter)
    }
  }
  
  #' ---------------------------------------------------------------------------
  #' 6. Attach settings and return the curve data.
  #' ---------------------------------------------------------------------------
  #' Retain the settings used to construct the curves.
  attr(roc_results, "settings") <- list(
    N = N,
    n = n,
    distributions = distributions,
    norm_config = norm_config,
    pretest_threshold = pretest_threshold,
    sig_levels = sig_levels,
    effect_size_H1 = effect_size_H1,
    effect_size_H0 = effect_size_H0,
    timestamp = Sys.time()
  )
  
  roc_results
}


#' -----------------------------------------------------------------------------
#' Plot full and zoomed power-versus-Type-I-error curves.
#'
#' For each distribution, the function produces:
#' 1. a full plot over the complete Type I error range;
#' 2. a zoomed plot over the range from zero to zoom_xlim; and
#' 3. a legend panel reporting the operating point nearest nominal_alpha.
#'
#' Curves are shown for test 1, test 2, and the adaptive procedure. The subtitle
#' reports the normality-assessment approach and pretest threshold used when the
#' curve data were generated.
#' -----------------------------------------------------------------------------
power_vs_error_roc_plot <- function(
    roc_results,
    nominal_alpha = 0.05,
    optimal_result = NULL,
    zoom_xlim = 0.10,
    method_labels = NULL
) {
  
  #' ---------------------------------------------------------------------------
  #' 1. Recover the stored simulation settings.
  #' ---------------------------------------------------------------------------
  #' Read settings attached by power_vs_error_roc_data().
  settings <- attr(roc_results, "settings")
  
  #' Use fallback values when the settings attribute is unavailable.
  if (is.null(settings)) {
    warning("roc_results is missing a settings attribute. Using defaults.", call. = FALSE
    )
    
    settings <- list(
      norm_config = list(method = "classical"),
      pretest_threshold = 0.05
    )
  }
  
  #' Read the threshold used to construct the adaptive p-values.
  pretest_threshold <- settings$pretest_threshold
  
  #' Use the supplied optimal threshold for display when available.
  if (!is.null(optimal_result) && !is.null(optimal_result$param_star)) {
    pretest_threshold <- optimal_result$param_star
  }
  
  #' Resolve a readable normality-assessment label.
  method_label <- if (!is.null(settings$norm_config$method)) {
    switch(
      tolower(settings$norm_config$method),
      classical = "Classical",
      custom = "Custom",
      settings$norm_config$method
    )
  } else {
    "Unknown"
  }
  
  #' ---------------------------------------------------------------------------
  #' 2. Configure the multi-panel layout.
  #' ---------------------------------------------------------------------------
  
  #' Identify the displayed distributions.
  distributions <- unique(roc_results$Distribution)
  n_distributions <- length(distributions)
  
  #' Use one column per distribution and three plot rows.
  layout_matrix <- matrix(seq_len(3L * n_distributions), nrow = 3L, 
                          ncol = n_distributions, byrow = TRUE)
  
  #' Preserve the active graphics settings.
  old_par <- graphics::par(no.readonly = TRUE)
  
  on.exit(
    {
      graphics::layout(1)
      graphics::par(old_par)
    },
    add = TRUE
  )
  
  #' Allocate larger rows to the full and zoomed plots.
  graphics::layout(layout_matrix, heights = c(1, 1, 0.28))
  
  #' Set margins and readable text sizes.
  graphics::par(
    mar = c(4, 4, 3.0, 1),
    oma = c(1.2, 0, 4.2, 0),
    mgp = c(2.2, 0.7, 0),
    cex.axis = 1.0,
    cex.lab = 1.2,
    font.lab = 2,
    font.axis = 2
  )
  
  #' ---------------------------------------------------------------------------
  #' 3. Define method colors and labels.
  #' ---------------------------------------------------------------------------
  
  #' Preserve consistent colors across all panels.
  test_colors <- c(test_1 = "red", test_2 = "blue", adaptive = "green")
  
  #' Resolve user-facing names for the internal method codes.
  method_names <- resolve_method_labels(
    names(test_colors),
    method_labels
  )
  
  #' Store nominal-alpha operating points for the legend row.
  legend_coordinates <- list()
  
  #' ---------------------------------------------------------------------------
  #' 4. Define the shared panel-drawing helper.
  #' ---------------------------------------------------------------------------
  
  plot_roc_panel <- function(
    distribution_data,
    x_limits,
    main_title,
    store_coordinates = FALSE
  ) {
    
    #' Initialize the plotting region.
    graphics::plot(
      NA,
      xlim = x_limits,
      ylim = c(0, 1),
      xlab = "P(Type I Error)",
      ylab = "Power",
      main = main_title,
      cex.axis = 1.1,
      cex.lab = 1.2
    )
    
    #' Add the equal-power-and-error reference line.
    graphics::abline(a = 0, b = 1, col = "gray80", lty = 2)
    #' Mark the nominal downstream significance level.
    graphics::abline(v = nominal_alpha, col = "red", lty = 2)
    
    #' Store operating points when requested.
    coordinates <- list()
    
    #' Draw one curve for each available method.
    for (method in names(test_colors)) {
      
      #' Skip methods absent from the current distribution.
      if (!method %in% distribution_data$Method) {
        next
      }
      
      #' Extract and order the current method data.
      method_data <- distribution_data[distribution_data$Method == method, ,drop = FALSE]
      method_data <- method_data[order(method_data$Alpha), ,drop = FALSE]
      
      #' Draw the power-versus-Type-I-error curve.
      graphics::lines(method_data$TypeI_error,
        method_data$Power,
        col = test_colors[method],
        lwd = 3
      )
      
      #' Locate the point nearest nominal_alpha.
      nominal_index <- which.min(abs(method_data$Alpha - nominal_alpha))
      
      #' Mark the nominal-alpha operating point.
      graphics::points(
        method_data$TypeI_error[nominal_index],
        method_data$Power[nominal_index],
        col = test_colors[method],
        pch = 16,
        cex = 2.0
      )
      
      #' Retain the coordinates for the legend panel.
      if (isTRUE(store_coordinates)) {
        coordinates[[method]] <- c(
          method_data$TypeI_error[nominal_index],
          method_data$Power[nominal_index]
        )
      }
    }
    
    if (isTRUE(store_coordinates)) {
      coordinates
    } else {
      invisible(NULL)
    }
  }
  
  #' ---------------------------------------------------------------------------
  #' 5. Draw the full-view panels.
  #' ---------------------------------------------------------------------------
  
  for (distribution in distributions) {
    
    #' Extract results for the current distribution.
    distribution_data <- roc_results[roc_results$Distribution == distribution, , drop = FALSE]
    
    #' Draw the complete Type I error range.
    coordinates <- plot_roc_panel(
      distribution_data = distribution_data,
      x_limits = c(0, 1),
      main_title = paste(distribution, "- Full View"),
      store_coordinates = TRUE
    )
    
    #' Save the nominal-alpha operating points.
    legend_coordinates[[distribution]] <- coordinates
  }
  
  #' ---------------------------------------------------------------------------
  #' 6. Draw the zoomed panels.
  #' ---------------------------------------------------------------------------
  
  for (distribution in distributions) {
    
    #' Extract results for the current distribution.
    distribution_data <- roc_results[roc_results$Distribution == distribution, , drop = FALSE]
    
    #' Restrict the horizontal range to the Type I error region of interest.
    plot_roc_panel(
      distribution_data = distribution_data,
      x_limits = c(0, zoom_xlim),
      main_title = paste(distribution, "- Zoomed View")
    )
  }
  
  #' ---------------------------------------------------------------------------
  #' 7. Draw one legend panel per distribution.
  #' ---------------------------------------------------------------------------
  
  #' Use minimal margins for the legend row.
  graphics::par(mar = c(0.2, 0.2, 0.2, 0.2))
  
  for (distribution in distributions) {
    
    #' Initialize an empty legend panel.
    graphics::plot.new()

    legend_text <- character(0)
    legend_colors <- character(0)
    
    #' Build one summary entry for each available method.
    for (method in names(test_colors)) {
      
      #' Skip methods without stored coordinates.
      if (!method %in% names(legend_coordinates[[distribution]])) {
        next
      }
      
      #' Retrieve Type I error and power at nominal_alpha.
      coordinates <- legend_coordinates[[distribution]][[method]]
      
      legend_text <- c(
        legend_text,
        sprintf("%s (alpha = %.3f): TIE = %.3f, Power = %.3f",
          method_names[method], nominal_alpha, coordinates[1L],
          coordinates[2L]
        )
      )
      
      legend_colors <- c(legend_colors, test_colors[method])
    }
    
    #' Display the distribution-specific operating-point summary.
    graphics::legend(
      "center",
      legend = legend_text,
      col = legend_colors,
      lwd = 4,
      pch = 16,
      pt.cex = 1.6,
      bty = "b",
      cex = 1.15,
      text.font = 2,
      title = paste(distribution, "Distribution"),
      title.adj = 0.5,
      title.cex = 1.2
    )
  }
  
  #' ---------------------------------------------------------------------------
  #' 8. Add common titles.
  #' ---------------------------------------------------------------------------
  
  #' Add the overall figure title.
  graphics::mtext(
    sprintf("ROC-Like Curves: Power vs Type I Error Rate (%s Approach)", method_label),
    side = 3, outer = TRUE, cex = 1.3, line = 2.4, font = 2
  )
  
  #' Report the nominal level and routing threshold.
  graphics::mtext(sprintf("Nominal alpha = %.3f | Pretest threshold = %.4f",
      nominal_alpha, pretest_threshold), side = 3, outer = TRUE, cex = 1.0, line = 0.65
  )
  
  #' Return the plotted data without printing them.
  invisible(roc_results)
}


#' -----------------------------------------------------------------------------
#' Apply the selected adaptive downstream test to one observed dataset.
#'
#' The function first evaluates the object used for normality assessment. It
#' then applies test 1 when all finite normality decision values exceed the
#' selected pretest threshold; otherwise, it applies test 2.
#'
#' The function returns only the p-value from the selected downstream test.
#' -----------------------------------------------------------------------------
run_adaptive_test_dual <- function(
    data,
    norm_config = list(method = "classical", config = list(norm_test = "SW")),
    pretest_threshold,
    fn_to_get_norm_obj = raw_data,
    fn_for_ds_test_1 = twosample_t_test,
    fn_for_ds_test_2 = Mann_whitney_U_test,
    ...
) {
  
  #' Require the routing threshold.
  if (missing(pretest_threshold) || is.null(pretest_threshold)) {
    stop("pretest_threshold must be provided.", call. = FALSE)
  }
  
  #' Require the normality-assessment method.
  if (is.null(norm_config$method)) {
    stop("norm_config must have a $method slot: 'classical' or 'custom'.", call. = FALSE)
  }
  
  #' Extract the sample or residual object used for normality assessment.
  normality_object <- fn_to_get_norm_obj(data)
  
  #' Compute classical p-values or custom normality scores.
  normality_decision <- assess_normality(
    normality_object,
    method = norm_config$method,
    config = norm_config$config
  )
  
  #' Select test 1 when the normality criterion is satisfied.
  use_test_1 <- use_parametric_test(
    normality_decision,
    threshold = pretest_threshold
  )
  
  #' Apply the selected downstream procedure.
  downstream_result <- if (use_test_1) {
    fn_for_ds_test_1(data)
  } else {
    fn_for_ds_test_2(data)
  }
  
  #' Return the selected test p-value.
  downstream_result$p.value
}

#' =============================================================================
#' 12. POWER AND TYPE I ERROR ACROSS SAMPLE SIZE AND EFFECT SIZE
#'
#' This section evaluates the Type I error and power of the fixed and adaptive
#' downstream procedures across sample sizes. The simulated p-values are stored
#' and reused to construct operating-characteristic summaries and AUC tables.
#' =============================================================================


#' -----------------------------------------------------------------------------
#' Estimate power and Type I error across sample sizes.
#'
#' For each distribution and sample size, the function generates the downstream
#' p-values and normality decisions once. It then evaluates the requested fixed
#' tests and adaptive procedure from those stored values.
#'
#' The adaptive procedure may use one common threshold, the threshold contained
#' in optimal_result, or named sample-size-specific values in threshold_by_n.
#' -----------------------------------------------------------------------------
run_ds_simulation <- function(
    sample_sizes = c(10, 20, 30, 40, 50),
    distributions = c("exponential", "normal"),
    N = 1000,
    alpha = 0.05,
    effect_size_H1 = 0.5,
    effect_size_H0 = zero_like(effect_size_H1),
    norm_config = list(method = "classical", config = list(norm_test = "SW")),
    pretest_threshold = 0.05,
    optimal_result = NULL,
    threshold_by_n = NULL,
    ds_test_methods = c("test_1", "test_2", "adaptive"),
    gen_data = onesample_data,
    get_parameters = onesample_parameters,
    fn_to_get_norm_obj = raw_data,
    fn_for_ds_test_1 = one_sample_t_test,
    fn_for_ds_test_2 = sign_test,
    ...
) {
  
  #' ---------------------------------------------------------------------------
  #' 1. Validate the simulation settings.
  #' ---------------------------------------------------------------------------
  
  #' Validate the number of simulation replicates.
  N <- validate_positive_integer(N, "N")
  
  #' Require at least one numeric sample size.
  if (!is.numeric(sample_sizes) || length(sample_sizes) == 0L) {
    stop("sample_sizes must contain positive integers.", call. = FALSE)
  }
  
  #' Validate each sample size.
  sample_sizes <- vapply(
    sample_sizes,
    validate_positive_integer,
    integer(1),
    name = "each sample size"
  )
  
  #' Remove duplicates and sort the sample sizes.
  sample_sizes <- sort(unique(sample_sizes))
  
  #' Validate the requested downstream methods.
  valid_methods <- c("test_1", "test_2", "adaptive")
  
  ds_test_methods <- match.arg(ds_test_methods, choices = valid_methods, several.ok = TRUE)
  
  #' Preserve the standard method order.
  ds_test_methods <- valid_methods[
    valid_methods %in% ds_test_methods
  ]
  
  #' Check whether named thresholds are supplied by sample size.
  use_per_n_thresholds <- !is.null(threshold_by_n) &&
    is.numeric(threshold_by_n) && !is.null(names(threshold_by_n))
  
  #' ---------------------------------------------------------------------------
  #' 2. Determine the source of the adaptive threshold.
  #' ---------------------------------------------------------------------------
  
  #' A single optimal threshold overrides the fallback scalar.
  if (!use_per_n_thresholds &&
      !is.null(optimal_result$param_star)) {
    
    pretest_threshold <- as.numeric(
      optimal_result$param_star)
    
    message(
      sprintf("Using optimal pretest threshold = %.4f", pretest_threshold)
    )
    
    #' Use named thresholds for individual sample sizes.
  } else if (use_per_n_thresholds) {
    
    message("Using sample-size-specific pretest thresholds.")
    
    #' Retain the supplied scalar threshold as a fallback.
  } else {
    
    warning(
      sprintf("No optimal_result supplied; using pretest_threshold = %.4f.", pretest_threshold),
      call. = FALSE
    )
  }
  
  #' Return the threshold assigned to one sample size.
  resolve_threshold <- function(current_n) {
    
    #' Use the common threshold when per-n thresholds are unavailable.
    if (!use_per_n_thresholds) {
      return(pretest_threshold)
    }
    
    #' Retrieve the named threshold for the current sample size.
    threshold <- threshold_by_n[as.character(current_n)]
    
    #' Fall back to the common threshold when no valid match exists.
    if (length(threshold) == 0L || !is.finite(threshold)) {
      
      warning(
        sprintf(
          paste0("No threshold_by_n value was found for n = %d; ","using %.4f."),
          current_n,
          pretest_threshold
        ),
        call. = FALSE
      )
      return(pretest_threshold)
    }
    as.numeric(threshold)
  }
  
  #' ---------------------------------------------------------------------------
  #' 3. Initialize the result containers.
  #' ---------------------------------------------------------------------------
  #' Store complete simulation output by distribution.
  simulation_results <- list()
  #' Store compact data frames for plots and summaries.
  plot_data <- list()
  
  #' ---------------------------------------------------------------------------
  #' 4. Process each distribution.
  #' ---------------------------------------------------------------------------
  
  for (distribution in distributions) {
    #' Display the current distribution.
    cat("\n=== Distribution:", distribution,"===\n")
    
    #' Allocate distribution-specific storage.
    distribution_results <- list(
      power = setNames(vector("list", length(ds_test_methods)), ds_test_methods),
      type1 = setNames(vector("list", length(ds_test_methods)), ds_test_methods),
      pvals = list(H0 = setNames(vector("list", length(ds_test_methods)), ds_test_methods),
        H1 = setNames(vector("list", length(ds_test_methods)), ds_test_methods)),
      thresholds_used = numeric(length(sample_sizes)),
      timing = list()
    )
    
    #' Allocate one result position per sample size.
    for (method in ds_test_methods) {
      
      distribution_results$power[[method]] <- numeric(length(sample_sizes))
      
      distribution_results$type1[[method]] <- numeric(length(sample_sizes))
      
      distribution_results$pvals$H0[[method]] <- vector("list", length(sample_sizes))
      
      distribution_results$pvals$H1[[method]] <- vector("list", length(sample_sizes))
    }
    
    #' Display progress across sample sizes.
    progress_bar <- utils::txtProgressBar(min = 0, max = length(sample_sizes), style = 3)
    
    #' -------------------------------------------------------------------------
    #' 5. Process each sample size.
    #' -------------------------------------------------------------------------
    
    for (sample_index in seq_along(sample_sizes)) {
      
      #' Current sample size.
      current_n <- sample_sizes[sample_index]
      #' Threshold used at the current sample size.
      current_threshold <- resolve_threshold(current_n)
      #' Record the selected threshold.
      distribution_results$thresholds_used[sample_index] <- current_threshold
      
      #' Start the computation timer.
      start_time <- Sys.time()
      
      #' Generate every p-value and normality score once for this n.
      stored_values <- generate_pval(
        Nsim = N,
        n = current_n,
        effect_size_H1 = effect_size_H1,
        effect_size_H0 = effect_size_H0,
        norm_config = norm_config,
        dist = distribution,
        gen_data = gen_data,
        get_parameters = get_parameters,
        fn_to_get_norm_obj = fn_to_get_norm_obj,
        fn_for_ds_test_1 = fn_for_ds_test_1,
        fn_for_ds_test_2 = fn_for_ds_test_2,
        ...
      )
      
      #' Apply the adaptive routing rule under H0.
      use_test1_H0 <- vapply(
        stored_values$norm_decision_H0,
        use_parametric_test,
        logical(1),
        threshold = current_threshold
      )
      
      #' Apply the adaptive routing rule under H1.
      use_test1_H1 <- vapply(
        stored_values$norm_decision_H1,
        use_parametric_test,
        logical(1),
        threshold = current_threshold
      )
      
      #' -----------------------------------------------------------------------
      #' 6. Construct fixed and adaptive p-value vectors.
      #' -----------------------------------------------------------------------
    
      #' P-values under H0.
      pvalues_H0 <- list(
        test_1 = stored_values$pval_ds_test1_H0,
        test_2 = stored_values$pval_ds_test2_H0,
        adaptive = ifelse(
          use_test1_H0,
          stored_values$pval_ds_test1_H0,
          stored_values$pval_ds_test2_H0
        )
      )
      
      #' P-values under H1.
      pvalues_H1 <- list(
        test_1 = stored_values$pval_ds_test1_H1,
        test_2 = stored_values$pval_ds_test2_H1,
        adaptive = ifelse(
          use_test1_H1,
          stored_values$pval_ds_test1_H1,
          stored_values$pval_ds_test2_H1
        )
      )
      
      #' -----------------------------------------------------------------------
      #' 7. Compute and store operating characteristics.
      #' -----------------------------------------------------------------------
      for (method in ds_test_methods) {
        
        #' Estimate Type I error.
        distribution_results$type1[[method]][sample_index] <- mean(
          pvalues_H0[[method]] <= alpha,
          na.rm = TRUE
        )
        
        #' Estimate power.
        distribution_results$power[[method]][sample_index] <- mean(
          pvalues_H1[[method]] <= alpha,
          na.rm = TRUE
        )
        
        #' Store the complete H0 p-value vector.
        distribution_results$pvals$H0[[method]][[sample_index]] <-
          pvalues_H0[[method]]
        
        #' Store the complete H1 p-value vector.
        distribution_results$pvals$H1[[method]][[sample_index]] <-
          pvalues_H1[[method]]
      }
      
      #' Record elapsed time for the current sample size.
      distribution_results$timing[[paste0("n=", current_n)]] <- as.numeric(
        difftime(
          Sys.time(),
          start_time,
          units = "secs"
        )
      )
      
      #' Update progress.
      utils::setTxtProgressBar(progress_bar, sample_index)
    }
    
    #' Close the sample-size progress bar.
    close(progress_bar)
    
    #' Store the evaluated sample sizes.
    distribution_results$sample_sizes <- sample_sizes
    #' Save the complete distribution results.
    simulation_results[[distribution]] <- distribution_results
    
    #' -------------------------------------------------------------------------
    #' 8. Create compact plotting tables.
    #' -------------------------------------------------------------------------
    
    #' Initialize the power table.
    power_data <- data.frame(n = sample_sizes, Distribution = distribution, stringsAsFactors = FALSE)
    
    #' Initialize the Type I error table.
    type1_data <- power_data
    
    #' Add one column for each requested method.
    for (method in ds_test_methods) {
      
      power_data[[method]] <- distribution_results$power[[method]]
      
      type1_data[[method]] <- distribution_results$type1[[method]]
    }
    
    #' Arrange the power-table columns.
    power_data <- power_data[, c("n", ds_test_methods, "Distribution"), drop = FALSE]
    
    #' Arrange the Type I error-table columns.
    type1_data <- type1_data[, c("n", ds_test_methods, "Distribution"), drop = FALSE]
    
    #' Store data for the current distribution.
    plot_data[[distribution]] <- list(
      power = power_data,
      type1 = type1_data
    )
    
    cat("Completed:", distribution, "\n")
  }
  
  #' ---------------------------------------------------------------------------
  #' 9. Combine results across distributions.
  #' ---------------------------------------------------------------------------
  
  #' Stack the power tables.
  combined_power <- do.call(rbind, lapply(plot_data, `[[`, "power"))
  
  #' Stack the Type I error tables.
  combined_type1 <- do.call(rbind, lapply(plot_data, `[[`, "type1"))
  
  #' Remove inherited row names.
  rownames(combined_power) <- NULL
  rownames(combined_type1) <- NULL
  
  #' ---------------------------------------------------------------------------
  #' 10. Compute normalized sample-size AUC summaries.
  #' ---------------------------------------------------------------------------
  
  #' Summarize each sample-size curve by its normalized trapezoidal area.
  auc_tables <- setNames(
    lapply(
      distributions,
      function(distribution) {
        
        #' Retrieve plotting data for the current distribution.
        distribution_data <- plot_data[[distribution]]
        
        data.frame(
          Method = ds_test_methods,
          
          AUC_Power = vapply(
            ds_test_methods,
            function(method) {
              compute_auc(sample_sizes, distribution_data$power[[method]], 
                          ensure_endpoints = FALSE, normalize = TRUE)
            },
            numeric(1)
          ),
          
          AUC_TypeI = vapply(
            ds_test_methods,
            function(method) {
              compute_auc(
                sample_sizes, distribution_data$type1[[method]], 
                ensure_endpoints = FALSE, normalize = TRUE)
            },
            numeric(1)
          ),
    
          row.names = NULL
        )
      }
    ),
    distributions
  )
  
  #' ---------------------------------------------------------------------------
  #' 11. Assemble the final result object.
  #' ---------------------------------------------------------------------------
  
  results <- list(
    raw = simulation_results,
    plot_data = list(
      combined_power = combined_power,
      combined_type1 = combined_type1
    ),
    auc_tables = auc_tables,
    summary = list(
      sample_sizes = sample_sizes,
      distributions = distributions,
      N = N,
      alpha = alpha,
      effect_size_H1 = effect_size_H1,
      effect_size_H0 = effect_size_H0,
      norm_config = norm_config,
      pretest_threshold = pretest_threshold,
      threshold_by_n = threshold_by_n,
      ds_test_methods = ds_test_methods,
      timestamp = Sys.time()
    )
  )
  
  #' Assign a class for downstream methods.
  class(results) <- "ds_simulation_results"
  
  results
}

#' -----------------------------------------------------------------------------
#' Plot power and Type I error by sample size.
#'
#' The function creates one power panel and one Type I error panel for each
#' distribution. It uses readable method labels, a shared legend, and a dynamic
#' Type I error axis that keeps the nominal significance level visible when the
#' observed error rates are conservative.
#'
#' The selected adaptive threshold and normality-assessment approach are shown
#' in the outer title when that information is available.
#' -----------------------------------------------------------------------------
plot_power_type1_results <- function(
    combined_power,
    combined_type1,
    ds_test_methods = c("test_1", "test_2", "adaptive"),
    distributions = c("exponential", "normal"),
    sample_sizes,
    test_alpha = 0.05,
    optimal_result = NULL,
    norm_config = NULL,
    method_labels = NULL,
    type1_ylim = NULL,
    line_width = 3,
    point_size = 1.1
) {
  
  #' ---------------------------------------------------------------------------
  #' 1. Resolve display labels and plotting identifiers.
  #' ---------------------------------------------------------------------------
  
  #' Read the selected adaptive threshold when available.
  threshold <- optimal_result$param_star %||% NULL
  
  #' Create the threshold subtitle.
  threshold_label <- if (is.null(threshold)) {
    "default threshold"
  } else {
    sprintf("threshold = %.4f", threshold)
  }
  
  #' Create a readable label for the normality-assessment approach.
  approach_label <- switch(
    tolower(norm_config$method %||% "custom"),
    classical = "Classical Approach",
    custom = "Custom Approach",
    norm_config$method %||% "Custom Approach"
  )
  
  #' Resolve user-facing names for the requested methods.
  display_labels <- unname(
    resolve_method_labels(
      ds_test_methods,
      method_labels
    )
  )
  
  #' Assign consistent colors and plotting symbols.
  colors <- c("red", "blue", "green")[seq_along(ds_test_methods)]
  shapes <- c(19, 17, 15)[seq_along(ds_test_methods)]
  
  #' ---------------------------------------------------------------------------
  #' 2. Configure the multi-panel layout.
  #' ---------------------------------------------------------------------------
  
  #' Use one power and one Type I error panel per distribution.
  n_distributions <- length(distributions)
  n_panels <- 2L * n_distributions
  #' Reserve the final panel for the shared legend.
  legend_panel <- n_panels + 1L
  #' Preserve the current graphics settings.
  old_par <- par(no.readonly = TRUE)
  
  on.exit(
    {
      layout(1)
      par(old_par)
    },
    add = TRUE
  )
  
  #' Arrange power panels, Type I error panels, and the legend row.
  layout(
    matrix(c(seq_len(n_panels), rep(legend_panel, n_distributions)),
      nrow = 3, byrow = TRUE), heights = c(1, 1, 0.25)
  )
  
  #' Set readable margins and text sizes.
  par(
    mar = c(3.6, 3.8, 2.8, 0.8),
    oma = c(0, 0, 3.0, 0),
    mgp = c(2.2, 0.60, 0),
    tcl = -0.2,
    cex.axis = 1.0,
    cex.lab = 1.2,
    font.lab = 2,
    font.axis = 2
  )
  
  #' Add a small horizontal margin around the sample-size values.
  x_range <- range(sample_sizes, finite = TRUE)
  
  x_padding <- if (diff(x_range) > 0) {
    0.03 * diff(x_range)
  } else {
    1
  }
  
  x_limits <- x_range + c(-x_padding, x_padding)
  
  #' ---------------------------------------------------------------------------
  #' 3. Plot power by sample size.
  #' ---------------------------------------------------------------------------
  
  for (distribution in distributions) {
    
    #' Extract power results for the current distribution.
    distribution_power <- combined_power[combined_power$Distribution == distribution, , drop = FALSE]
    
    #' Initialize the power panel.
    plot(
      NA,
      xlim = x_limits,
      ylim = c(0, 1.01),
      xlab = "Sample Size",
      ylab = "Power",
      main = paste("Power -", distribution),
      cex.main = 1.3,
      font.main = 2
    )
    
    #' Add a light reference grid.
    grid(col = "gray88", lty = "dotted")
    
    #' Draw one power curve for each requested method.
    for (j in seq_along(ds_test_methods)) {
      
      method <- ds_test_methods[j]
      
      lines(
        distribution_power$n,
        distribution_power[[method]],
        type = "b",
        col = colors[j],
        pch = shapes[j],
        lwd = line_width,
        cex = point_size
      )
    }
  }
  
  #' ---------------------------------------------------------------------------
  #' 4. Determine the Type I error vertical scale.
  #' ---------------------------------------------------------------------------
  #' Keep the nominal alpha line visible when all estimates are conservative.
  if (is.null(type1_ylim)) {
    
    #' Find the largest observed Type I error across methods.
    observed_maximum <- suppressWarnings(
      max(combined_type1[, ds_test_methods, drop = FALSE], na.rm = TRUE)
    )
    
    #' Use the nominal level when no finite estimate is available.
    if (!is.finite(observed_maximum)) {
      observed_maximum <- test_alpha
    }
    
    #' Add headroom above both the nominal level and observed estimates.
    type1_upper <- max(test_alpha * 1.25, observed_maximum * 1.10, 0.075)
    type1_ylim <- c(0, type1_upper)
  }
  
  #' ---------------------------------------------------------------------------
  #' 5. Plot Type I error by sample size.
  #' ---------------------------------------------------------------------------
  
  for (distribution in distributions) {
    
    #' Extract Type I error results for the current distribution.
    distribution_type1 <- combined_type1[combined_type1$Distribution == distribution, , drop = FALSE]
    
    #' Initialize the Type I error panel.
    plot(
      NA,
      xlim = x_limits,
      ylim = type1_ylim,
      xlab = "Sample Size",
      ylab = "Type I Error Rate",
      main = paste("Type I Error -", distribution),
      cex.main = 1.3,
      font.main = 2
    )
    
    #' Add a light reference grid.
    grid(col = "gray88", lty = "dotted")
    
    #' Mark the nominal Type I error level.
    abline(h = test_alpha, col = "gray50", lty = 2, lwd = 1.4)
    
    #' Draw one Type I error curve for each requested method.
    for (j in seq_along(ds_test_methods)) {
      
      method <- ds_test_methods[j]
      
      lines(
        distribution_type1$n,
        distribution_type1[[method]],
        type = "b",
        col = colors[j],
        pch = shapes[j],
        lwd = line_width,
        cex = point_size
      )
    }
  }
  
  #' ---------------------------------------------------------------------------
  #' 6. Add the shared legend and outer title.
  #' ---------------------------------------------------------------------------
  
  #' Use minimal margins for the legend row.
  par(mar = c(0.1, 0.1, 0.1, 0.1))
  #' Initialize the legend panel.
  plot.new()
  
  #' Draw the shared method legend.
  legend(
    "center",
    legend = display_labels,
    col = colors,
    pch = shapes,
    lwd = line_width,
    horiz = TRUE,
    bty = "n",
    cex = 1.25,
    title = "Test Methods"
  )
  
  #' Display the assessment approach and threshold.
  mtext(
    sprintf("Test Methods Comparison (%s) | %s", approach_label, threshold_label),
    side = 3, outer = TRUE, line = 0.5, cex = 1.35, font = 2
  )
  
  #' Return the plotted data without printing them.
  invisible(
    list(
      power = combined_power,
      type1 = combined_type1
    )
  )
}


#' -----------------------------------------------------------------------------
#' Estimate power across an effect-size grid.
#'
#' For each distribution and effect size, the function generates the downstream
#' p-values and normality decisions once. It then estimates power for test 1,
#' test 2, and the adaptive procedure from the stored H1 p-values.
#'
#' A scalar effect is inserted into the final component when the data generator
#' requires a vector-valued effect, such as a one-way ANOVA group-location
#' pattern. The function returns the formatted effect-size results produced by
#' format_power_effect_results().
#' -----------------------------------------------------------------------------
perform_ds_power_by_effect <- function(
    fixed_n = 10,
    effect_sizes = seq(0.1, 1.0, by = 0.1),
    distributions = c("exponential", "normal"),
    Nsim = 1000,
    norm_config = list(method = "classical", config = list(norm_test = "SW")
    ),
    pretest_threshold = 0.05,
    optimal_result = NULL,
    gen_data,
    get_parameters,
    fn_to_get_norm_obj,
    fn_for_ds_test_1,
    fn_for_ds_test_2,
    ds_test_methods = c("test_1", "test_2", "adaptive"),
    test_alpha = 0.05,
    effect_template = NULL,
    verbose = TRUE,
    ...
) {
  
  #' ---------------------------------------------------------------------------
  #' 1. Validate the simulation settings.
  #' ---------------------------------------------------------------------------
  #' Validate the fixed sample size and simulation count.
  fixed_n <- validate_positive_integer(fixed_n, "fixed_n")
  Nsim <- validate_positive_integer(Nsim, "Nsim")
  
  #' Convert the effect-size grid to numeric values.
  effect_sizes <- as.numeric(effect_sizes)
  
  #' Require at least one finite effect size.
  if (length(effect_sizes) == 0L || any(!is.finite(effect_sizes))) {
    stop("effect_sizes must contain finite numeric values.", call. = FALSE)
  }
  
  #' Validate the requested downstream methods.
  valid_methods <- c("test_1", "test_2", "adaptive")
  
  ds_test_methods <- match.arg(ds_test_methods, choices = valid_methods, several.ok = TRUE)
  
  #' Preserve the standard method order.
  ds_test_methods <- valid_methods[valid_methods %in% ds_test_methods]
  
  #' ---------------------------------------------------------------------------
  #' 2. Resolve the adaptive threshold.
  #' ---------------------------------------------------------------------------
  
  #' Use the selected threshold when a trade-off result is available.
  if (!is.null(optimal_result$param_star)) {
    pretest_threshold <- as.numeric(optimal_result$param_star)
    
    message(
      sprintf("Using optimal pretest threshold = %.4f", pretest_threshold)
    )
  } else {
    warning(
      sprintf("No optimal_result supplied; using pretest_threshold = %.4f.",pretest_threshold),
      call. = FALSE
    )
  }
  
  #' ---------------------------------------------------------------------------
  #' 3. Determine the effect-size structure.
  #' ---------------------------------------------------------------------------
  
  #' Infer the effect shape when no template is supplied.
  if (is.null(effect_template)) {
    base_parameters <- get_parameters(fixed_n, dist = distributions[1L],...)
    effect_template <- base_parameters$effect_size
  }
  
  #' Number of elements in the effect specification.
  effect_length <- length(effect_template)
  #' Store results separately for each distribution.
  results_by_distribution <- list()
  
  #' ---------------------------------------------------------------------------
  #' 4. Process each distribution.
  #' ---------------------------------------------------------------------------
  for (distribution in distributions) {
  
    #' Display the current distribution when requested.
    if (isTRUE(verbose)) {
      cat("\n=== Distribution:", distribution, "===\n")
    }
    
    #' Allocate power estimates for each method.
    power_by_method <- setNames(
      lapply(ds_test_methods, function(method) numeric(length(effect_sizes))),
      ds_test_methods
    )
    
    #' Allocate H1 p-value storage for verification.
    pvalues_H1 <- setNames(
      lapply(ds_test_methods, function(method) vector("list", length(effect_sizes))),
      ds_test_methods
    )
    
    #' Store elapsed time for each effect size.
    timing <- list()
    #' Display progress across effect sizes.
    progress_bar <- utils::txtProgressBar(min = 0, max = length(effect_sizes), style = 3)
    
    #' -------------------------------------------------------------------------
    #' 5. Process each effect size.
    #' -------------------------------------------------------------------------
    
    for (effect_index in seq_along(effect_sizes)) {
      
      #' Current scalar effect value.
      current_value <- effect_sizes[effect_index]
      #' Start the computation timer.
      start_time <- Sys.time()
      
      #' Place a scalar effect in the final component of a vector-valued effect.
      current_effect <- if (effect_length > 1L) {
        effect_vector <- rep(0, effect_length)
        effect_vector[effect_length] <- current_value
        effect_vector
      } else {
        current_value
      }
      
      #' Generate all downstream p-values and normality decisions once.
      simulation_values <- generate_pval(
        Nsim = Nsim,
        n = fixed_n,
        effect_size_H1 = current_effect,
        effect_size_H0 = zero_like(current_effect),
        norm_config = norm_config,
        dist = distribution,
        gen_data = gen_data,
        get_parameters = get_parameters,
        fn_to_get_norm_obj = fn_to_get_norm_obj,
        fn_for_ds_test_1 = fn_for_ds_test_1,
        fn_for_ds_test_2 = fn_for_ds_test_2,
        ...
      )
      
      #' Apply the adaptive routing rule under H1.
      use_test1_H1 <- vapply(
        simulation_values$norm_decision_H1,
        use_parametric_test,
        logical(1),
        threshold = pretest_threshold
      )
      
      #' -----------------------------------------------------------------------
      #' 6. Construct the H1 p-value vectors.
      #' -----------------------------------------------------------------------
      method_pvalues <- list(
        test_1 = simulation_values$pval_ds_test1_H1,
        test_2 = simulation_values$pval_ds_test2_H1,
        adaptive = ifelse(
          use_test1_H1,
          simulation_values$pval_ds_test1_H1,
          simulation_values$pval_ds_test2_H1
        )
      )
      
      #' -----------------------------------------------------------------------
      #' 7. Estimate power and retain the p-values.
      #' -----------------------------------------------------------------------
      for (method in ds_test_methods) {
      
        #' Estimate power at the selected significance level.
        power_by_method[[method]][effect_index] <- mean(
          method_pvalues[[method]] <= test_alpha,
          na.rm = TRUE
        )
        
        #' Store the complete H1 p-value vector.
        pvalues_H1[[method]][[effect_index]] <- method_pvalues[[method]]
      }
      
      #' Record elapsed time for the current effect size.
      timing[[paste0("effect=", current_value)]] <- as.numeric(
        difftime(Sys.time(), start_time, units = "secs"))
      
      #' Update progress.
      utils::setTxtProgressBar(progress_bar, effect_index)
    }
    
    #' Close the effect-size progress bar.
    close(progress_bar)
    
    #' -------------------------------------------------------------------------
    #' 8. Store the distribution-specific results.
    #' -------------------------------------------------------------------------
    results_by_distribution[[distribution]] <- list(
      power = power_by_method,
      pvalues = list(H1 = pvalues_H1),
      timing = timing,
      effect_sizes = effect_sizes,
      pretest_threshold = pretest_threshold,
      norm_config = norm_config,
      fixed_n = fixed_n,
      distribution = distribution,
      test_alpha = test_alpha,
      effect_length = effect_length,
      eff_len = effect_length
    )
  }
  
  #' ---------------------------------------------------------------------------
  #' 9. Format and return the combined effect-size results.
  #' ---------------------------------------------------------------------------
  format_power_effect_results(
    results_by_dist = results_by_distribution,
    distributions = distributions,
    methods = ds_test_methods
  )
}

#' -----------------------------------------------------------------------------
#' Convert effect-size simulation results to a plot-ready combined data frame.
#'
#' The function combines the distribution-specific power results returned by
#' perform_ds_power_by_effect() into one long-format data frame. It also retains
#' the raw results and shared simulation settings for later plotting and
#' verification.
#' -----------------------------------------------------------------------------
format_power_effect_results <- function(
    results_by_dist,
    distributions,
    methods
) {
  
  #' Store one power table for each distribution.
  plot_data <- list()
  
  #' Combine method-specific power estimates within each distribution.
  for (dist in distributions) {

    #' Retrieve the current distribution results.
    dist_results <- results_by_dist[[dist]]
    
    #' Initialize the distribution-specific power table.
    power_df <- data.frame(
      effect_value = dist_results$effect_sizes,
      Distribution = dist,
      n = dist_results$fixed_n,
      stringsAsFactors = FALSE
    )
    
    #' Add one power column for each requested method.
    for (method in methods) {
      power_df[[method]] <- dist_results$power[[method]]
    }
    
    #' Store the completed table.
    plot_data[[dist]] <- power_df
  }
  
  #' Stack all distribution-specific power tables.
  combined_power <- do.call(rbind, plot_data)
  #' Remove inherited row names.
  rownames(combined_power) <- NULL
  
  #' Read shared settings from the first distribution.
  first <- results_by_dist[[distributions[1]]]
  
  #' Return raw results, plotting data, and common settings.
  list(raw_results = results_by_dist,
    plot_data = list(power_by_effect = combined_power),
    parameters = list(
      distributions = distributions,
      methods = methods,
      norm_config = first$norm_config,
      pretest_threshold = first$pretest_threshold,
      fixed_n = first$fixed_n,
      effect_sizes = first$effect_sizes,
      test_alpha = first$test_alpha,
      eff_len = first$eff_len,
      timestamp = Sys.time()
    )
  )
}


#' -----------------------------------------------------------------------------
#' Plot power across an effect-size grid for every distribution.
#'
#' The function creates one power panel for each distribution and a shared legend
#' below the panels. Curves are drawn for the fixed downstream tests and the
#' adaptive procedure. When the effect specification is vector-valued, the
#' horizontal axis is labeled as the final group effect size.
#'
#' A vertical reference line may be added through reference_effect, and the
#' selected adaptive threshold is displayed in the outer subtitle when available.
#' -----------------------------------------------------------------------------
plot_power_by_effect_size <- function(
    combined_power,
    fixed_n,
    optimal_result = NULL,
    norm_config = NULL,
    eff_len = 1,
    method_labels = NULL,
    reference_effect = NULL,
    line_width = 3,
    point_size = 1.1
) {
  
  #' ---------------------------------------------------------------------------
  #' 1. Validate the plotting data.
  #' ---------------------------------------------------------------------------

  #' Require the columns used to define the horizontal axis and panels.
  required_columns <- c("effect_value", "Distribution")
  
  if (!all(required_columns %in% names(combined_power))) {
    stop("combined_power must contain effect_value and Distribution.", call. = FALSE)
  }
  
  #' Identify the available downstream methods.
  methods <- intersect(c("test_1", "test_2", "adaptive"), names(combined_power)
  )
  
  #' Identify the displayed distributions.
  distributions <- unique(as.character(combined_power$Distribution))
  
  #' Resolve readable method labels.
  display_labels <- unname(resolve_method_labels(methods, method_labels))
  
  #' Assign consistent plotting colors and symbols.
  colors <- c("red", "blue", "green")[seq_along(methods)]
  shapes <- c(19, 17, 15)[seq_along(methods)]
  
  #' Read the selected threshold when available.
  threshold <- optimal_result$param_star %||% NA_real_
  
  #' Resolve the normality-assessment approach label.
  approach <- switch(
    tolower(norm_config$method %||% "custom"),
    classical = "Classical Approach",
    custom = "Custom Approach",
    norm_config$method %||% "Custom Approach"
  )
  
  #' ---------------------------------------------------------------------------
  #' 2. Configure the multi-panel layout.
  #' ---------------------------------------------------------------------------
  
  #' Preserve the active graphics settings.
  old_par <- par(no.readonly = TRUE)
  
  on.exit(
    {
      layout(1)
      par(old_par)
    },
    add = TRUE
  )
  
  #' Use one plot panel per distribution and one shared legend panel.
  n_distributions <- length(distributions)
  
  layout(matrix(
      c(seq_len(n_distributions), rep(n_distributions + 1L, n_distributions)),
      nrow = 2, byrow = TRUE), heights = c(1, 0.22))
  
  #' Set readable margins and text sizes.
  par(
    mar = c(4.2, 4.2, 3, 1),
    oma = c(0.8, 0, 3.8, 0),
    mgp = c(2.2, 0.7, 0),
    cex.axis = 1.0,
    cex.lab = 1.2,
    font.lab = 2,
    font.axis = 2
  )
  
  #' Determine the shared effect-size range.
  effect_limits <- range(
    combined_power$effect_value,
    finite = TRUE
  )
  
  #' ---------------------------------------------------------------------------
  #' 3. Plot power against effect size.
  #' ---------------------------------------------------------------------------
  for (distribution in distributions) {
  
    #' Extract results for the current distribution.
    distribution_power <- combined_power[combined_power$Distribution == distribution, , drop = FALSE]
    
    #' Initialize the distribution-specific panel.
    plot(
      NA,
      xlim = effect_limits,
      ylim = c(0, 1),
      xlab = if (eff_len > 1L) {
        "Last Group Effect Size"
      } else {
        "Effect Size"
      },
      ylab = "Power",
      main = paste(distribution, "(n =", fixed_n, ")")
    )
    
    #' Add a light reference grid.
    grid(col = "gray88", lty = "dotted")
    
    #' Add the optional reference effect.
    if (!is.null(reference_effect) && is.finite(reference_effect)) {
      abline(v = reference_effect, col = "gray45", lty = 2, lwd = 1.4)
    }
    
    #' Draw one power curve for each available method.
    for (j in seq_along(methods)) {
      
      method <- methods[j]
  
      lines(
        distribution_power$effect_value,
        distribution_power[[method]],
        type = "b",
        col = colors[j],
        pch = shapes[j],
        lwd = line_width,
        cex = point_size
      )
    }
  }
  
  #' ---------------------------------------------------------------------------
  #' 4. Add the shared legend.
  #' ---------------------------------------------------------------------------
  #' Use minimal margins for the legend panel.
  par(mar = c(0.2, 0.2, 0.2, 0.2))
  #' Initialize the legend panel.
  plot.new()
  #' Draw the shared method legend.
  legend(
    "center",
    title = "Test Methods",
    legend = display_labels,
    col = colors,
    pch = shapes,
    lwd = line_width,
    horiz = TRUE,
    bty = "n",
    cex = 1.2
  )
  
  #' ---------------------------------------------------------------------------
  #' 5. Add the overall title and threshold subtitle.
  #' ---------------------------------------------------------------------------
  #' Resolve the threshold subtitle.
  subtitle <- if (is.finite(threshold)) {
    sprintf("Threshold = %.4f", threshold)
  } else {
    "Stored simulation threshold"
  }
  
  #' Add the main outer title.
  mtext(paste0("Power by Effect Size (", approach,")"),
    side = 3, outer = TRUE, line = 2.0, cex = 1.35, font = 2)
  #' Add the threshold subtitle.
  mtext(subtitle, side = 3, outer = TRUE, line = 0.6, cex = 0.95)
  #' Return the plotted data without printing them.
  invisible(combined_power)
}

#' =============================================================================
#' 13. MAIN ANALYSIS ORCHESTRATOR
#'
#' This section coordinates the complete normality-pretesting workflow. It
#' validates the requested phases, runs the selected analyses, saves all plots,
#' and stores the resulting objects and reproducibility settings in one
#' structured results object.
#' =============================================================================

#' -----------------------------------------------------------------------------
#' Run selected analysis phases and save one structured results object.
#'
#' The available phases are:
#' 1. normality-test ROC analysis;
#' 2. threshold trade-off analysis;
#' 3. power-versus-Type-I-error curves;
#' 4. power and Type I error across sample sizes; and
#' 5. power across effect sizes.
#'
#' Phase 2 may select one common threshold or separate thresholds by sample
#' size. When Phase 2 is skipped, pretest_threshold must be supplied for any
#' later phase that requires adaptive routing.
#' -----------------------------------------------------------------------------
run_simulation <- function(Nsim           = 1e3,
                           N_tradeoff     = 1e3,
                           test_type      = "one_sample_t_vs_sign_test",
                           distributions  = c("exponential", "normal"),
                           
                           norm_config    = list(method = "classical", config = list(norm_test = "SW")),
                           threshold_grid = seq(0, 1, by = 0.025),
                           normality_roc_grid = NULL,
                           
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
                           method_labels   = NULL,
                           type1_ylim      = NULL,
                           reference_effect = NULL,
                           
                           norm_test      = c("SW","SF","LF","KS","JB","SKEW", "KURT","DAP","AD","CVM"),
                           selected_tests = c("SW","SF","LF","JB","DAP","AD"),
                           
                           phases            = 1:5,
                           pretest_threshold = NULL,
                           per_n_thresholds  = FALSE,
                           output_dir        = NULL,
                           save_results      = TRUE,
                           save_compression  = "gzip",
                           
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
  if (6L %in% phases) {
    warning("Phase 6 is no longer a separate analytical phase; results are saved ",
      "during finalization. The value 6 has been ignored.",
      call. = FALSE
    )
    phases <- setdiff(phases, 6L)
  }
  
  valid_phases <- 1:5
  bad_phases <- setdiff(phases, valid_phases)
  
  if (length(bad_phases) > 0L) {
    stop("Invalid phase(s): ", paste(bad_phases, collapse = ", "),
      ". Valid analytical phases are 1 through 5.",
      call. = FALSE
    )
  }
  
  if (length(phases) == 0L) {
    stop("At least one analytical phase must be selected.", call. = FALSE)
  }
  
  #' --------------------------------------------------------------------------
  #' Normality ROC grid
  #' --------------------------------------------------------------------------
  #' Use a separate grid for the normality ROC curve.
  #' If not supplied, use threshold_grid to preserve the old behavior.
  if (is.null(normality_roc_grid)) {
    normality_roc_grid <- threshold_grid
  }
  
  normality_roc_grid <- sort(unique(as.numeric(normality_roc_grid)))
  
  if (length(normality_roc_grid) < 2L || any(!is.finite(normality_roc_grid))) {
    stop("normality_roc_grid must contain at least two finite numeric values.",
      call. = FALSE
    )
  }
  
  threshold_grid <- sort(unique(as.numeric(threshold_grid)))
  
  if (length(threshold_grid) < 2L || any(!is.finite(threshold_grid))) {
    stop("threshold_grid must contain at least two finite numeric values.",
      call. = FALSE
    )
  }
  
  run_phase <- function(p) p %in% phases
  
  #' Phases 3, 4, 5 all need param_star.
  #' It comes from Phase 2 unless the user bypasses Phase 2 with pretest_threshold.
  needs_threshold <- any(c(3, 4, 5) %in% phases)
  skip_phase2     <- !run_phase(2)
  
  if (needs_threshold && skip_phase2) {
    if (is.null(pretest_threshold) || !is.numeric(pretest_threshold)) {
      stop("Phase 2 (trade-off analysis) is skipped but Phases ",
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
      message(sprintf("Phase 2 skipped. Using user-supplied pretest_threshold = %.4f.",
        pretest_threshold
      ))
    }
  }
  
  #' --------------------------------------------------------------------------
  #' STEP 3: TIMING / INIT
  #' --------------------------------------------------------------------------
  start_time     <- Sys.time()
  method_display <- switch(method, classical = "Classical", custom = "Custom")
  
  cat("\n", strrep("=", 60), "\n", sep = "")
  cat("STARTING SIMULATION WORKFLOW (", method_display, " Approach)\n", sep = "")
  cat("Test type    :", test_type, "\n")
  cat("Distributions:", paste(distributions, collapse = " vs "), "\n")
  cat("Phases       :", paste(sort(phases), collapse = ", "), "\n")
  cat(strrep("=", 60), "\n\n")
  
  eff       <- split_effect_size(effect_size)
  effect_H0 <- eff$H0
  effect_H1 <- eff$H1
  
  #' If single_n is not supplied, use the first sample size.
  if (is.null(single_n)) {
    single_n <- sample_sizes[1]
    message(sprintf("single_n was not supplied. Using single_n = %d.", single_n))
  }
  
  if (isTRUE(per_n_thresholds) && !single_n %in% sample_sizes) {
    stop("single_n must be included in sample_sizes when per_n_thresholds = TRUE.",
      call. = FALSE
    )
  }
  
  #' Initialize threshold-related objects.
  #' These are filled in Phase 2 or when Phase 2 is skipped.
  optimal_result <- NULL
  param_star     <- NULL
  threshold_by_n <- NULL
  
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
      normality_roc_grid = normality_roc_grid,
      effect_sizes_plot = effect_sizes_plot,
      sig_levels        = sig_levels,
      ds_test_methods   = ds_test_methods,
      method_labels     = method_labels,
      type1_ylim        = type1_ylim,
      reference_effect  = reference_effect,
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
      output_dir        = output_dir,
      save_results      = save_results,
      save_compression  = save_compression,
      timestamp_start   = start_time
    ),
    objects = list(),
    files   = list()
  )
  
  output_parent <- output_dir %||% file.path(getwd(), "results")
  results$root_dir <- file.path(output_parent, test_type)
  dir.create(results$root_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (!dir.exists(results$root_dir)) {
    stop("Could not create output directory: ", results$root_dir, call. = FALSE)
  }
  
  results$root_dir <- normalizePath(results$root_dir, winslash = "/", mustWork = TRUE)
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
      selection_note = sprintf("User-supplied pretest_threshold = %.4f (Phase 2 skipped)", threshold),
      feasible              = NA,
      inflation_controlled  = NA,
      loss_controlled       = NA,
      selection_stage       = NA_integer_,
      power_loss            = NA,
      power_gain            = NA,
      inflation_normal      = NA,
      inflation_non_normal  = NA,
      worst_case_inflation  = NA,
      effect_size           = effect_H1,
      all_metrics  = list(
        param_values  = threshold,
        distributions = distributions,
        effect_size   = effect_H1,
        param_name    = "threshold"
      )
    )
  }
  
  #' ==========================================================================
  #' PHASE 1: NORMALITY PRETEST ROC ANALYSIS
  #' ==========================================================================
  if (run_phase(1)) {
    cat("\n1. GENERATING NORMALITY TEST ROC CURVES\n")
    
    results$files$norm_roc <- file.path(results$root_dir, paste0(test_type, "_normality_roc.pdf")
    )
    phase1_out <- tryCatch(
      save_pdf_plot(
        file = results$files$norm_roc,
        width = 9,
        height = 7,
        plot_function = function() {
          plot_normality_roc_wrapper(
            norm_config = norm_config,
            single_n = single_n,
            threshold_grid = normality_roc_grid,
            distributions = distributions,
            Nsim = Nsim,
            norm_test = norm_test,
            selected_tests = selected_tests,
            gen_data = gen_data,
            get_parameters = get_parameters,
            fn_to_get_norm_obj = fn_to_get_norm_obj,
            center_by = center_by,
            ...
          )
        }
      ),
      error = function(e) {
        message("  Phase 1 plotting error: ", conditionMessage(e))
        NULL
      }
    )
    
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
    #' Each sample size receives an independent trade-off simulation.
    #' A sequential loop keeps the progress bar and multi-page PDF reliable.
    #' Replicate-level data generation remains inside generate_pval().
    #' ------------------------------------------------------------------
    if (per_n_thresholds) {
      cat("\n2. PERFORMING TRADE-OFF ANALYSIS (per sample size)\n")
      
      #' Capture additional generator arguments once for the per-n loop.
      extra_args <- list(...)
      
      progress_bar <- utils::txtProgressBar(min = 0, max = length(sample_sizes), style = 3)
      
      results$files$tradeoff <- file.path(results$root_dir, paste0(test_type, "_tradeoff_analysis.pdf")
      )
      
      grDevices::pdf(results$files$tradeoff, width = 10, height = 8)
      
      #' Run one independent trade-off analysis for each sample size.
      per_n_results <- tryCatch(
        lapply(seq_along(sample_sizes), function(sample_index) {
          n_i <- sample_sizes[sample_index]
          
          analysis_i <- do.call(
            perform_adaptive_analysis,
            c(
              list(
                Nsim = N_tradeoff,
                n = n_i,
                effect_size_H1 = effect_H1,
                effect_size_H0 = effect_H0,
                distributions = distributions,
                norm_config = norm_config,
                threshold_grid = threshold_grid,
                test_alpha = test_alpha,
                gen_data = gen_data,
                get_parameters = get_parameters,
                fn_to_get_norm_obj = fn_to_get_norm_obj,
                fn_for_ds_test_1 = fn_for_ds_test_1,
                fn_for_ds_test_2 = fn_for_ds_test_2,
                center_by = center_by
              ),
              extra_args
            )
          )
          
          metrics_i <- compute_roc_metrics(
            typeI_rates = analysis_i$typeI_rates,
            power_rates = analysis_i$power_rates,
            test_alpha = test_alpha
          )
          
          optimal_i <- select_optimal_parameter(
            distributions = distributions,
            param_grid = analysis_i$param_grid,
            metrics = metrics_i,
            alpha = test_alpha,
            tol_pos = tol_pos,
            loss_tol = loss_tol,
            param_name = analysis_i$param_name,
            effect_size = effect_H1
          )
          
          plot_tradeoff_results_generic(
            optimal_result = optimal_i,
            outer_title = sprintf("Trade-off (%s Approach) | n = %d | threshold = %.4f",
              method_display,
              n_i,
              optimal_i$param_star
            )
          )
          
          utils::setTxtProgressBar(progress_bar, sample_index)
          
          list(
            n_i = n_i,
            analysis = analysis_i,
            metrics = metrics_i,
            optimal = optimal_i
          )
        }),
        finally = {
          close(progress_bar)
          grDevices::dev.off()
        }
      )
      
      #' Keep the result objects in ascending sample-size order.
      per_n_results <- per_n_results[order(sapply(per_n_results, `[[`, "n_i"))]
      
      #' Collect optimal thresholds and result objects into named vectors/lists.
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
        cat(sprintf("   n = %d  ->  threshold = %.4f\n", sample_sizes[i], threshold_by_n[i]))
      }
      
      #' For Phases 3 and 5, which operate at single_n, extract that entry.
      optimal_result <- optimal_result_by_n[[as.character(single_n)]]
      param_star     <- threshold_by_n[as.character(single_n)]
      
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
        param_name    = analysis_ds_tests$param_name,
        effect_size   = effect_H1
      )
      results$objects$optimal_result <- optimal_result
      param_star <- optimal_result$param_star
      results$objects$param_star <- param_star
      threshold_by_n <- NULL
      
      results$files$tradeoff <- file.path(results$root_dir, paste0(test_type, "_tradeoff_analysis.pdf")
      )
      save_pdf_plot(
        file = results$files$tradeoff,
        width = 10,
        height = 8,
        plot_function = function() {
          plot_tradeoff_results_generic(
            optimal_result = optimal_result,
            outer_title = sprintf("Power and Type I Error Trade-off (%s Approach)",
              method_display
            )
          )
        }
      )
      cat("  Saved:", results$files$tradeoff, "\n")
      cat("  Optimal threshold =", round(param_star, 4), "\n")
    }
    
  } else {
    cat("\n2. TRADE-OFF ANALYSIS  [skipped]\n")
    threshold_by_n <- NULL
    
    if (needs_threshold) {
      #' Accept either a named vector or a scalar for pretest_threshold
      if (is.numeric(pretest_threshold) && !is.null(names(pretest_threshold))) {
        threshold_by_n <- pretest_threshold
        single_threshold <- pretest_threshold[as.character(single_n)]
        
        if (length(single_threshold) == 0L || !is.finite(single_threshold)) {
          warning("No named pretest threshold was supplied for single_n; using the first value.",
            call. = FALSE
          )
          single_threshold <- pretest_threshold[[1L]]
        }
        
        param_star <- as.numeric(single_threshold)
        optimal_result <- make_bypass_optimal(param_star)
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
    
    roc_data <- power_vs_error_roc_data(
      N                  = Nsim,
      n                  = single_n,
      distributions      = distributions,
      norm_config        = norm_config,
      effect_size_H1     = effect_H1,
      effect_size_H0     = effect_H0,
      optimal_result     = optimal_result,
      sig_levels         = sig_levels,
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
    save_pdf_plot(
      file = results$files$power_error_roc,
      width = 10,
      height = 10,
      plot_function = function() {
        power_vs_error_roc_plot(
          roc_results = roc_data,
          nominal_alpha = test_alpha,
          optimal_result = optimal_result,
          method_labels = method_labels
        )
      }
    )
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
      threshold_by_n     = threshold_by_n,
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
    save_pdf_plot(
      file = results$files$power_type1,
      width = 10,
      height = 8,
      plot_function = function() {
        plot_power_type1_results(
          combined_power = sim_output$plot_data$combined_power,
          combined_type1 = sim_output$plot_data$combined_type1,
          ds_test_methods = ds_test_methods,
          distributions = distributions,
          sample_sizes = sample_sizes,
          test_alpha = test_alpha,
          optimal_result = optimal_result,
          norm_config = norm_config,
          method_labels = method_labels,
          type1_ylim = type1_ylim
        )
      }
    )
    cat("  Saved:", results$files$power_type1, "\n")
    
    #' Print AUC summary to console
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
    
    base_paras      <- get_parameters(single_n, dist = distributions[1], center_by = center_by, ...)
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
    save_pdf_plot(
      file = results$files$power_by_effect,
      width = 10,
      height = 7,
      plot_function = function() {
        plot_power_by_effect_size(
          combined_power = sim_power_list$plot_data$power_by_effect,
          fixed_n = single_n,
          optimal_result = optimal_result,
          norm_config = norm_config,
          eff_len = eff_len,
          method_labels = method_labels,
          reference_effect = reference_effect
        )
      }
    )
    cat("  Saved:", results$files$power_by_effect, "\n")
    
  } else {
    cat("\n5. POWER VS EFFECT SIZE  [skipped]\n")
  }
  
  #' ==========================================================================
  #' FINALIZE
  #' ==========================================================================
  end_time   <- Sys.time()
  total_time <- difftime(end_time, start_time, units = "mins")
  
  results$params$timestamp_end         <- end_time
  results$params$total_runtime_minutes <- as.numeric(total_time)
  
  #' ==========================================================================
  #' FINALIZATION: SAVE THE STRUCTURED RESULTS OBJECT
  #' ==========================================================================
  rdata_results <- file.path(
    results$root_dir,
    paste0(test_type, "_results.RData")
  )
  
  #' Store the path before saving so it is present after the file is reloaded.
  results$files$results_rdata <- normalizePath(rdata_results, winslash = "/", mustWork = FALSE)
  
  if (isTRUE(save_results)) {
    cat("\nSAVING RESULTS\n")
    
    valid_compression <- identical(save_compression, FALSE) ||
      identical(save_compression, TRUE) || save_compression %in% c("gzip", "bzip2", "xz")
    
    if (!valid_compression) {
      stop("save_compression must be FALSE, TRUE, 'gzip', 'bzip2', or 'xz'.",
        call. = FALSE
      )
    }
    
    save(
      results,
      file = rdata_results,
      compress = save_compression
    )
    
    cat("  Saved results:", rdata_results, "\n")
  } else {
    cat("\nRESULT SAVING SKIPPED\n")
  }
  
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
#' 14. OPTIONAL CUSTOM AND ML NORMALITY INTEGRATION
#'
#' This section provides utilities for combining classical tests, locating and
#' loading external ML resources, applying sample-size-specific trained models,
#' and constructing custom normality configurations for ROC analysis and
#' adaptive routing.
#' =============================================================================

#' -----------------------------------------------------------------------------
#' Combine Shapiro-Wilk and Anderson-Darling p-values using Fisher's method.
#'
#' The two p-values are transformed into a chi-square statistic with four
#' degrees of freedom and returned through the common named p-value interface.
#' -----------------------------------------------------------------------------
fisher_combined <- function(x) {
  #' Require the package used for the Anderson-Darling test.
  require_namespace("nortest", "the Fisher-combined SW and AD test")
  
  #' Compute the two component p-values.
  p_sw <- stats::shapiro.test(x)$p.value
  p_ad <- nortest::ad.test(x)$p.value
  
  #' Form Fisher's combined chi-square statistic.
  statistic <- -2 * sum(log(c(p_sw, p_ad)))
  
  #' Return the combined upper-tail p-value.
  c(p.value = stats::pchisq(statistic, df = 4, lower.tail = FALSE))
}

#' -----------------------------------------------------------------------------
#' Return the directory containing the currently sourced script when available.
#'
#' The function searches active call frames for a source-file path. When no
#' sourced script can be identified, it returns the current working directory.
#' -----------------------------------------------------------------------------
get_script_dir <- function() {
  source_files <- vapply(
    sys.frames(),
    function(frame) frame$ofile %||% "",
    character(1)
  )
  
  source_files <- source_files[nzchar(source_files)]
  
  if (length(source_files) > 0L) {
    #' Return the sourced script directory.
    return(dirname(normalizePath(tail(source_files, 1L), mustWork = TRUE)))
  }
  
  #' Fall back to the current working directory.
  normalizePath(getwd(), mustWork = TRUE)
}

#' -----------------------------------------------------------------------------
#' Return the first existing file from a vector of candidate paths.
#'
#' Candidate paths are expanded, deduplicated, and checked in their supplied
#' order. NULL is returned when none of the files exists.
#' -----------------------------------------------------------------------------
first_existing_file <- function(paths) {
  #' Expand, remove blanks, and deduplicate candidate paths.
  paths <- unique(path.expand(paths[nzchar(paths)]))
  
  #' Retain files that currently exist.
  existing_paths <- paths[file.exists(paths)]
  
  if (length(existing_paths) == 0L) {
    return(NULL)
  }
  
  #' Return the first valid candidate in normalized form.
  normalizePath(existing_paths[1L], mustWork = TRUE)
}

#' -----------------------------------------------------------------------------
#' Resolve the ML framework and trained-model files without changing directories.
#'
#' Explicit arguments, R options, environment variables, expected project
#' locations, and interactive file selection are checked in that order. The
#' function returns normalized paths to the framework script and trained-model
#' RData file.
#' -----------------------------------------------------------------------------
resolve_ml_files <- function(framework_file = NULL,
                             model_file     = NULL,
                             project_dir    = NULL,
                             allow_choose   = interactive()) {
  #' Resolve the project directory from the argument, option, or environment.
  project_dir <- project_dir %||% getOption(
    "ml.normality.project_dir",
    Sys.getenv("ML_NORMALITY_PROJECT_DIR", unset = "")
  )
  
  #' Resolve an explicitly configured ML framework path.
  framework_file <- framework_file %||% getOption(
    "ml.normality.framework_file",
    Sys.getenv("ML_FRAMEWORK_FILE", unset = "")
  )
  
  #' Resolve an explicitly configured trained-model path.
  model_file <- model_file %||% getOption(
    "ml.normality.model_file",
    Sys.getenv("ML_MODEL_FILE", unset = "")
  )
  
  #' Build the directories searched for ML resources.
  search_roots <- unique(c(if (nzchar(project_dir)) project_dir else character(0), get_script_dir(), getwd()))
  
  #' Retain existing search directories only.
  search_roots <- path.expand(search_roots[dir.exists(path.expand(search_roots))])
  
  #' Construct expected framework-file locations.
  framework_candidates <- c(
    framework_file,
    unlist(
      lapply(search_roots, function(root) {
        file.path(
          root,
          c("ML_framework_func.R", "ML_framework.R")
        )
      }),
      use.names = FALSE
    )
  )
  
  #' Select the first framework file found.
  resolved_framework <- first_existing_file(framework_candidates)
  
  #' Prompt the user when automatic framework discovery fails.
  if (is.null(resolved_framework) && isTRUE(allow_choose)) {
    message("Select the ML framework function file.")
    resolved_framework <- normalizePath(file.choose(), mustWork = TRUE)
  }
  
  if (is.null(resolved_framework)) {
    stop(
      "The ML framework file could not be found. Supply framework_file or set ",
      "options(ml.normality.framework_file = '<path>').",
      call. = FALSE
    )
  }
  
  #' Search for trained models relative to the resolved framework.
  framework_dir <- dirname(resolved_framework)
  model_candidates <- c(
    model_file,
    file.path(
      framework_dir,
      "ml_oob_workflow_full_user_features",
      "stage2_final_refit",
      "trained_models.RData"
    ),
    file.path(framework_dir, "stage2_final_refit", "trained_models.RData"),
    file.path(framework_dir, "trained_models.RData")
  )
  
  #' Select the first expected trained-model file found.
  resolved_model <- first_existing_file(model_candidates)
  
  #' Search recursively only when the expected model paths are absent.
  if (is.null(resolved_model)) {
    discovered_models <- list.files(
      framework_dir,
      pattern = "^trained_models\\.RData$",
      recursive = TRUE,
      full.names = TRUE
    )
    
    #' Prefer models saved in the final-refit stage.
    preferred_models <- discovered_models[
      grepl("stage2_final_refit", discovered_models, fixed = TRUE)
    ]
    
    if (length(preferred_models) == 1L) {
      resolved_model <- normalizePath(preferred_models, mustWork = TRUE)
    } else if (length(discovered_models) == 1L) {
      resolved_model <- normalizePath(discovered_models, mustWork = TRUE)
    } else if (length(preferred_models) > 1L) {
      stop(
        "Multiple final-refit model files were found. Supply model_file explicitly.",
        call. = FALSE
      )
    }
  }
  
  #' Prompt the user when automatic model discovery fails.
  if (is.null(resolved_model) && isTRUE(allow_choose)) {
    message("Select the final trained_models.RData file.")
    resolved_model <- normalizePath(file.choose(), mustWork = TRUE)
  }
  
  if (is.null(resolved_model)) {
    stop(
      "The trained-model file could not be found. Supply model_file or set ",
      "options(ml.normality.model_file = '<path>').",
      call. = FALSE
    )
  }
  
  #' Return both verified file paths.
  list(framework_file = resolved_framework, model_file = resolved_model)
}

#' -----------------------------------------------------------------------------
#' Source ML functions and load trained models only when an ML analysis is requested.
#'
#' The framework script is sourced into an isolated environment, required
#' functions are verified, and trained_models is loaded from the selected RData
#' file. When requested sample sizes are unavailable, each size is mapped to the
#' closest trained sample size instead of stopping the analysis.
#' -----------------------------------------------------------------------------
load_ml_resources <- function(
    framework_file = NULL,
    model_file = NULL,
    project_dir = NULL,
    required_sample_sizes = NULL,
    required_functions = c(
      "calculate_features",
      "prepare_prediction_features",
      "predict_prob_normal_fast"
    )
) {
  #' Locate the framework script and trained-model file.
  files <- resolve_ml_files(
    framework_file = framework_file,
    model_file = model_file,
    project_dir = project_dir
  )
  
  #' Isolate sourced ML functions from the global workspace.
  ml_env <- new.env(parent = globalenv())
  
  #' Source the ML framework into the isolated environment.
  sys.source(files$framework_file, envir = ml_env, chdir = TRUE)
  
  #' Identify required ML functions missing from the sourced script.
  missing_functions <- required_functions[
    !vapply(
      required_functions,
      exists,
      logical(1),
      envir = ml_env,
      mode = "function",
      inherits = FALSE
    )
  ]
  
  if (length(missing_functions) > 0L) {
    stop(
      "The ML framework is missing function(s): ",
      paste(missing_functions, collapse = ", "),
      call. = FALSE
    )
  }
  
  #' Load trained objects into a separate temporary environment.
  model_environment <- new.env(parent = emptyenv())
  loaded_objects <- load(files$model_file, envir = model_environment)
  
  if (!"trained_models" %in% loaded_objects) {
    stop(
      "The model file does not contain trained_models.",
      call. = FALSE
    )
  }
  
  #' Extract the saved sample-size-specific model bundles.
  trained_models <- model_environment$trained_models
  
  if (!is.list(trained_models) || length(trained_models) == 0L) {
    stop(
      "trained_models must be a non-empty list of sample-size-specific bundles.",
      call. = FALSE
    )
  }
  
  #' Read numeric sample sizes from the trained-model list names.
  available_keys <- names(trained_models)
  available_sizes <- suppressWarnings(as.integer(available_keys))
  valid_sizes <- is.finite(available_sizes)
  
  if (!any(valid_sizes)) {
    stop(
      "The trained-model object does not contain numeric sample-size names.",
      call. = FALSE
    )
  }
  
  available_keys <- available_keys[valid_sizes]
  available_sizes <- available_sizes[valid_sizes]
  trained_models <- trained_models[available_keys]
  
  #' Sort model bundles by their trained sample size.
  size_order <- order(available_sizes)
  available_sizes <- available_sizes[size_order]
  available_keys <- available_keys[size_order]
  trained_models <- trained_models[available_keys]
  
  #' Map requested sizes to their closest available trained sizes.
  if (!is.null(required_sample_sizes)) {
    required_sample_sizes <- sort(
      unique(as.integer(required_sample_sizes))
    )
    
    if (length(required_sample_sizes) == 0L ||
        any(!is.finite(required_sample_sizes)) ||
        any(required_sample_sizes <= 0L)) {
      stop(
        "required_sample_sizes must contain positive integers.",
        call. = FALSE
      )
    }
    
    sample_size_map <- data.frame(
      requested_n = required_sample_sizes,
      selected_n = vapply(
        required_sample_sizes,
        function(current_n) {
          size_distance <- abs(available_sizes - current_n)
          closest_sizes <- available_sizes[
            size_distance == min(size_distance)
          ]
          
          #' Use the smaller size when two trained sizes are equally close.
          min(closest_sizes)
        },
        integer(1)
      ),
      stringsAsFactors = FALSE
    )
  } else {
    sample_size_map <- data.frame(
      requested_n = available_sizes,
      selected_n = available_sizes,
      stringsAsFactors = FALSE
    )
  }
  
  sample_size_map$exact_match <-
    sample_size_map$requested_n == sample_size_map$selected_n
  
  #' Check model availability only in bundles that may be used.
  selected_keys <- unique(as.character(sample_size_map$selected_n))
  
  available_models <- Reduce(
    intersect,
    lapply(
      trained_models[selected_keys],
      function(bundle) names(bundle$models)
    )
  )
  
  if (length(available_models) == 0L) {
    stop(
      "No common ML model type is available across the selected model bundles.",
      call. = FALSE
    )
  }
  
  #' Return the sourced functions, models, size mapping, and resolved paths.
  list(
    ml_env = ml_env,
    trained_models = trained_models,
    available_models = available_models,
    available_sample_sizes = available_sizes,
    sample_size_map = sample_size_map,
    framework_file = files$framework_file,
    model_file = files$model_file
  )
}

#' -----------------------------------------------------------------------------
#' Retrieve the exact or closest available trained ML model bundle.
#'
#' The target size is the length of the normality-assessment object unless
#' custom_sample_size is supplied. When no exact model exists, the bundle trained
#' at the closest available sample size is selected. If two sizes are equally
#' close, the smaller trained size is used.
#' -----------------------------------------------------------------------------
get_ml_bundle <- function(
    x,
    trained_models,
    custom_sample_size = NULL
) {
  #' Use the requested size or infer it from the assessed numeric vector.
  target_n <- custom_sample_size %||% length(x)
  target_n <- as.integer(target_n)
  
  if (length(target_n) != 1L ||
      !is.finite(target_n) ||
      target_n <= 0L) {
    stop(
      "The target ML sample size must be one positive integer.",
      call. = FALSE
    )
  }
  
  if (!is.list(trained_models) || length(trained_models) == 0L) {
    stop(
      "trained_models must be a non-empty list.",
      call. = FALSE
    )
  }
  
  #' Read valid numeric sample sizes from the model-bundle names.
  available_keys <- names(trained_models)
  available_sizes <- suppressWarnings(as.integer(available_keys))
  valid_sizes <- is.finite(available_sizes)
  
  if (!any(valid_sizes)) {
    stop(
      "The trained-model object does not contain numeric sample-size names.",
      call. = FALSE
    )
  }
  
  available_keys <- available_keys[valid_sizes]
  available_sizes <- available_sizes[valid_sizes]
  
  #' Identify the closest trained sample size.
  size_distance <- abs(available_sizes - target_n)
  closest_sizes <- available_sizes[
    size_distance == min(size_distance)
  ]
  
  #' Resolve equal-distance ties using the smaller trained size.
  selected_n <- min(closest_sizes)
  selected_key <- available_keys[
    match(selected_n, available_sizes)
  ]
  
  list(
    bundle = trained_models[[selected_key]],
    target_n = target_n,
    selected_n = selected_n,
    exact_match = selected_n == target_n
  )
}

#' -----------------------------------------------------------------------------
#' Return P(Normal) from one fitted ML model using its saved preprocessing.
#'
#' The function selects the exact or closest trained sample-size bundle and then
#' delegates feature preparation and probability prediction to the ML framework.
#' -----------------------------------------------------------------------------
predict_ml_prob_normal <- function(
    x,
    trained_models,
    model_name,
    ml_env,
    custom_sample_size = NULL
) {
  if (missing(ml_env) || !is.environment(ml_env)) {
    stop(
      "ml_env must be supplied from load_ml_resources().",
      call. = FALSE
    )
  }
  
  #' Select the exact or closest trained model bundle.
  bundle_information <- get_ml_bundle(
    x = x,
    trained_models = trained_models,
    custom_sample_size = custom_sample_size
  )
  
  model_bundle <- bundle_information$bundle
  
  if (!model_name %in% names(model_bundle$models)) {
    stop(
      "Model '",
      model_name,
      "' is unavailable in the bundle trained for n = ",
      bundle_information$selected_n,
      ".",
      call. = FALSE
    )
  }
  
  #' Apply the saved preprocessing and return P(Normal).
  probability_normal <- ml_env$predict_prob_normal_fast(
    x = x,
    model_bundle = model_bundle,
    model_name = model_name
  )
  
  #' Retain the size-selection details without changing the numeric return value.
  attr(probability_normal, "target_n") <- bundle_information$target_n
  attr(probability_normal, "selected_n") <- bundle_information$selected_n
  attr(probability_normal, "exact_match") <- bundle_information$exact_match
  
  probability_normal
}

#' -----------------------------------------------------------------------------
#' Classify one or more samples using exact or closest trained ML model bundles.
#'
#' Each sample is converted to the required feature representation and evaluated
#' by the requested models. The output records the observed sample size, the
#' trained model size used, model-specific probabilities and classes, and an
#' optional majority-vote class.
#' -----------------------------------------------------------------------------
classify_sample_prob <- function(
    sample_data,
    trained_models,
    ml_env,
    model_name = "RF",
    use_majority_vote = TRUE,
    decision_threshold = 0.50,
    custom_sample_size = NULL
) {
  #' Convert all supported inputs to a named list of samples.
  sample_information <- convert_to_sample_list(sample_data)
  
  #' Validate the requested model codes.
  valid_models <- c("LR", "RF", "ANN", "GBM", "SVM", "KNN")
  invalid_models <- setdiff(model_name, valid_models)
  
  if (length(invalid_models) > 0L) {
    stop(
      "Invalid model(s): ",
      paste(invalid_models, collapse = ", "),
      call. = FALSE
    )
  }
  
  #' Classify each sample separately.
  results <- lapply(
    seq_along(sample_information$samples),
    function(i) {
      #' Extract the current sample and display name.
      sample_values <- sample_information$samples[[i]]
      sample_name <- sample_information$sample_names[i]
      
      #' Select the exact or closest trained model bundle.
      bundle_information <- get_ml_bundle(
        x = sample_values,
        trained_models = trained_models,
        custom_sample_size = custom_sample_size
      )
      
      model_bundle <- bundle_information$bundle
      
      #' Verify that all requested models exist in the selected bundle.
      missing_models <- setdiff(
        model_name,
        names(model_bundle$models)
      )
      
      if (length(missing_models) > 0L) {
        stop(
          "The following model(s) are unavailable for sample '",
          sample_name,
          "' in the bundle trained for n = ",
          bundle_information$selected_n,
          ": ",
          paste(missing_models, collapse = ", "),
          ".",
          call. = FALSE
        )
      }
      
      available_models <- model_name
      
      #' Apply the preprocessing saved with the selected model bundle.
      prepared_features <- ml_env$prepare_prediction_features(
        x = sample_values,
        model_bundle = model_bundle
      )
      
      #' Initialize one output row and record the size mapping.
      output <- data.frame(
        Sample_Name = sample_name,
        Sample_Size = length(sample_values),
        ML_Target_Size = bundle_information$target_n,
        ML_Model_Size = bundle_information$selected_n,
        Exact_Size_Match = bundle_information$exact_match,
        stringsAsFactors = FALSE
      )
      
      #' Store model-specific non-Normal probabilities for voting.
      non_normal_probabilities <- numeric(
        length(available_models)
      )
      names(non_normal_probabilities) <- available_models
      
      #' Generate probabilities and classes from each model.
      for (current_model in available_models) {
        predicted_probabilities <- predict(
          model_bundle$models[[current_model]],
          newdata = prepared_features,
          type = "prob"
        )
        
        #' Extract class probabilities from the prediction output.
        probability_normal <- as.numeric(
          predicted_probabilities[, "Normal"]
        )
        probability_non_normal <- as.numeric(
          predicted_probabilities[, "Non_Normal"]
        )
        
        #' Store probabilities and the threshold-based class.
        output[[paste0(current_model, "_Prob_Normal")]] <-
          probability_normal
        
        output[[paste0(current_model, "_Prob_Non_Normal")]] <-
          probability_non_normal
        
        output[[paste0(current_model, "_Class")]] <- if (
          probability_non_normal >= decision_threshold
        ) {
          "Non_Normal"
        } else {
          "Normal"
        }
        
        non_normal_probabilities[current_model] <-
          probability_non_normal
      }
      
      #' Add an ensemble decision when several models are requested.
      if (isTRUE(use_majority_vote) &&
          length(available_models) > 1L) {
        mean_non_normal <- mean(
          non_normal_probabilities
        )
        
        output$Majority_Vote_Class <- if (
          mean_non_normal >= decision_threshold
        ) {
          "Non_Normal"
        } else {
          "Normal"
        }
        
        output$Proportion_Non_Normal <- mean_non_normal
      }
      
      output
    }
  )
  
  #' Combine sample-level rows into one result table.
  output <- do.call(rbind, results)
  rownames(output) <- NULL
  
  #' Assign a class for downstream methods.
  class(output) <- c(
    "ml_normality_result",
    "data.frame"
  )
  
  output
}

#' -----------------------------------------------------------------------------
#' Return one P(Normal) score for each requested ML model for ROC analysis.
#'
#' A named probability is returned for each requested model. When requested, an
#' additional MajorityVote score is calculated as the mean P(Normal).
#' -----------------------------------------------------------------------------
ml_fn_roc <- function(x,
                      trained_models,
                      ml_env,
                      model_name,
                      use_majority_vote = FALSE,
                      custom_sample_size = NULL,
                      ...) {
  #' Obtain one P(Normal) value from each requested model.
  probabilities <- vapply(
    model_name,
    function(current_model) {
      predict_ml_prob_normal(
        x = x,
        trained_models = trained_models,
        model_name = current_model,
        ml_env = ml_env,
        custom_sample_size = custom_sample_size
      )
    },
    numeric(1)
  )
  
  #' Preserve model names for separate ROC curves.
  names(probabilities) <- model_name
  
  #' Add the mean model probability as an optional ensemble score.
  if (isTRUE(use_majority_vote) && length(probabilities) > 1L) {
    probabilities <- c(probabilities, MajorityVote = mean(probabilities))
  }
  
  probabilities
}

#' -----------------------------------------------------------------------------
#' Return one scalar P(Normal) score for adaptive ML pretest routing.
#'
#' Adaptive routing requires exactly one selected model. The returned probability
#' is named ML_Prob_Normal and is compared directly with the routing threshold.
#' -----------------------------------------------------------------------------
ml_fn_decision <- function(x,
                           trained_models,
                           ml_env,
                           model_name,
                           custom_sample_size = NULL,
                           ...) {
  if (length(model_name) != 1L) {
    stop("ml_fn_decision() requires exactly one model name.", call. = FALSE)
  }
  
  #' Obtain P(Normal) from the selected routing model.
  probability_normal <- predict_ml_prob_normal(
    x = x,
    trained_models = trained_models,
    model_name = model_name,
    ml_env = ml_env,
    custom_sample_size = custom_sample_size
  )
  
  #' Return the named scalar expected by the adaptive framework.
  setNames(as.numeric(probability_normal), "ML_Prob_Normal")
}

#' -----------------------------------------------------------------------------
#' Create custom ML configurations for ROC analysis and adaptive routing.
#'
#' The ROC configuration may contain several models and an optional majority
#' vote. The decision configuration contains exactly one chosen model for
#' threshold-based routing in the adaptive procedure.
#' -----------------------------------------------------------------------------
make_ml_normality_configs <- function(ml_resources,
                                      chosen_model = "RF",
                                      models_for_roc = NULL,
                                      use_majority_vote = TRUE) {
  #' Use all shared models unless a smaller ROC set is requested.
  models_for_roc <- models_for_roc %||% ml_resources$available_models
  
  #' Identify requested ROC models that are unavailable.
  unavailable_roc_models <- setdiff(models_for_roc, ml_resources$available_models)
  
  if (length(unavailable_roc_models) > 0L) {
    stop(
      "Unavailable ROC model(s): ",
      paste(unavailable_roc_models, collapse = ", "),
      call. = FALSE
    )
  }
  
  #' Require the routing model to be available across the selected model bundles.
  if (!chosen_model %in% ml_resources$available_models) {
    stop("The chosen ML model is unavailable: ", chosen_model, call. = FALSE)
  }
  
  #' Build separate configurations for ROC comparison and adaptive routing.
  list(
    roc = list(
      method = "custom",
      config = list(
        fn = ml_fn_roc,
        trained_models = ml_resources$trained_models,
        ml_env = ml_resources$ml_env,
        model_name = models_for_roc,
        use_majority_vote = use_majority_vote
      )
    ),
    decision = list(
      method = "custom",
      config = list(
        fn = ml_fn_decision,
        trained_models = ml_resources$trained_models,
        ml_env = ml_resources$ml_env,
        model_name = chosen_model
      )
    )
  )
}
