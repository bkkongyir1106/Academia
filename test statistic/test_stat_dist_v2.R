## ============================================================
## Selective Inference in Normality Pre-testing
## Test Statistic and Adaptive Testing Simulation
##
## Purpose:
##   Compare the distribution of test statistics under:
##     1. Unconditional t-test
##     2. Conditional t-test after passing Shapiro-Wilk
##     3. Adaptive test
##
## Adaptive procedure:
##   One-sample:
##     - If Shapiro-Wilk passes: use one-sample t-test statistic
##     - If Shapiro-Wilk rejects: use standardized sign statistic
##
##   Two-sample:
##     - If both samples pass Shapiro-Wilk: use two-sample t-test statistic
##     - If either sample rejects: use standardized Wilcoxon rank-sum statistic
##
## Type I error:
##   Computed using the statistic tail-area rule:
##
##       P(T < lower critical value) + P(T > upper critical value)
##
## Output:
##   - Tail-area Type I error summary table
##   - One-sample density plots for all distributions
##   - Two-sample density plots for all distributions
##   - RData results files
## ============================================================


## ============================================================
## 0. Setup
## ============================================================

rm(list = ls())

setwd("~/Desktop/OSU/Research/Pretest-Simulation/Selective Inference")

set.seed(12345)


## ============================================================
## 1. Simulation settings
## ============================================================

distributions <- c(
  "Normal",
  "Exponential",
  "t",
  "Chi-Square",
  "Uniform",
  "LogNormal"
)

N     <- 1e5
alpha <- 0.05
n     <- 10

df_one_sample <- n - 1
df_two_sample <- 2 * n - 2


## ============================================================
## 2. Data generation
## ============================================================
## Data are centered by the theoretical median by default and
## scaled by the theoretical standard deviation.
##
## This matches the median-null formulation used in the adaptive
## one-sample simulations.
## ============================================================

generate_data <- function(n,
                          dist,
                          center_by   = "mean",
                          standardize = TRUE) {
  
  ## ----------------------------------------------------------
  ## Input validation
  ## ----------------------------------------------------------
  if (!is.numeric(n) || length(n) != 1 || n <= 0 || n != round(n)) {
    stop("n must be a single positive integer.")
  }
  
  dist      <- tolower(trimws(dist))
  center_by <- tolower(trimws(center_by))
  
  if (!center_by %in% c("median", "mean")) {
    stop("center_by must be either 'median' or 'mean'.")
  }
  
  ## Allow common spelling variations
  dist <- gsub("-", "_", dist)
  
  if (dist == "log_normal") {
    dist <- "lognormal"
  }
  
  if (dist %in% c("chi_square", "chisquare", "chi_sq")) {
    dist <- "chi_square"
  }
  
  ## ----------------------------------------------------------
  ## Distribution-specific quantities
  ## ----------------------------------------------------------
  if (dist == "normal") {
    
    samples <- rnorm(n, mean = 0, sd = 1)
    
    theoretical_mean   <- 0
    theoretical_median <- 0
    theoretical_sd     <- 1
    
  } else if (dist == "t") {
    
    df_t <- 3
    
    samples <- rt(n, df = df_t)
    
    theoretical_mean   <- 0
    theoretical_median <- 0
    theoretical_sd     <- sqrt(df_t / (df_t - 2))
    
  } else if (dist == "uniform") {
    
    a <- 0
    b <- 1
    
    samples <- runif(n, min = a, max = b)
    
    theoretical_mean   <- (a + b) / 2
    theoretical_median <- (a + b) / 2
    theoretical_sd     <- (b - a) / sqrt(12)
    
  } else if (dist == "exponential") {
    
    samples <- rexp(n, rate = 1)
    
    theoretical_mean   <- 1
    theoretical_median <- log(2)
    theoretical_sd     <- 1
    
  } else if (dist == "chi_square") {
    
    df_chi <- 7
    
    samples <- rchisq(n, df = df_chi)
    
    theoretical_mean   <- df_chi
    theoretical_median <- qchisq(0.5, df = df_chi)
    theoretical_sd     <- sqrt(2 * df_chi)
    
  } else if (dist == "lognormal") {
    
    samples <- rlnorm(n, meanlog = 0, sdlog = 1)
    
    theoretical_mean   <- exp(0.5)
    theoretical_median <- 1
    theoretical_sd     <- sqrt((exp(1) - 1) * exp(1))
    
  } else {
    
    stop(
      "Unsupported distribution: ", dist,
      ". Use 'Normal', 't', 'Uniform', 'Exponential', ",
      "'Chi-Square', or 'LogNormal'."
    )
  }
  
  ## ----------------------------------------------------------
  ## Return raw data if standardization is not requested
  ## ----------------------------------------------------------
  if (!standardize) {
    
    attr(samples, "distribution") <- dist
    attr(samples, "standardized") <- FALSE
    
    return(samples)
  }
  
  ## ----------------------------------------------------------
  ## Center and scale
  ## ----------------------------------------------------------
  center_value <- if (center_by == "median") {
    theoretical_median
  } else {
    theoretical_mean
  }
  
  samples <- (samples - center_value) / theoretical_sd
  
  ## ----------------------------------------------------------
  ## Attach metadata
  ## ----------------------------------------------------------
  attr(samples, "distribution") <- dist
  attr(samples, "center_by")    <- center_by
  attr(samples, "center_value") <- center_value
  attr(samples, "scale_factor") <- theoretical_sd
  attr(samples, "standardized") <- TRUE
  
  samples
}


## ============================================================
## 3. Core test helpers
## ============================================================

passes_shapiro <- function(x, alpha = 0.05) {
  shapiro.test(x)$p.value > alpha
}


get_one_sample_tstat <- function(x) {
  unname(t.test(x)$statistic)
}


get_one_sample_t_pvalue <- function(x) {
  t.test(x)$p.value
}


get_two_sample_tstat <- function(x, y) {
  
  ## Pooled two-sample t-test.
  ## This aligns the statistic with the reference df = 2n - 2.
  unname(t.test(x, y, var.equal = TRUE)$statistic)
}


get_two_sample_t_pvalue <- function(x, y) {
  t.test(x, y, var.equal = TRUE)$p.value
}


## ============================================================
## 4. Randomized sign test p-value
## ============================================================
## Used for the adaptive one-sample p-value branch.
## The plotted adaptive statistic still uses the standardized
## sign statistic, because the density plot requires a statistic.
## ============================================================

randomized_sign_test <- function(x, mu0 = 0) {
  
  centered_x <- x - mu0
  nonzero_n  <- sum(centered_x != 0)
  
  if (nonzero_n == 0L) {
    return(runif(1))
  }
  
  positive_count   <- sum(centered_x > 0)
  random_u         <- runif(1)
  probability_mass <- dbinom(positive_count, nonzero_n, 0.5)
  
  left_tail <- pbinom(
    positive_count - 1L,
    nonzero_n,
    0.5,
    lower.tail = TRUE
  ) + random_u * probability_mass
  
  right_tail <- pbinom(
    positive_count,
    nonzero_n,
    0.5,
    lower.tail = FALSE
  ) + random_u * probability_mass
  
  min(1, 2 * min(left_tail, right_tail))
}


get_two_sample_wilcox_pvalue <- function(x, y) {
  wilcox.test(x, y, exact = FALSE)$p.value
}


## ============================================================
## 5. Adaptive test statistic helpers
## ============================================================
## The adaptive statistic is a mixed statistic:
##   - t-statistic in the parametric branch
##   - standardized sign statistic in the one-sample fallback
##   - standardized Wilcoxon statistic in the two-sample fallback
##
## These fallback statistics are used only so that the adaptive
## procedure has a numeric statistic for density plotting.
## ============================================================

# Sign test statistic helpers
get_signed_sign_statistic <- function(x, mu0 = 0) {
  
  centered_x <- x - mu0
  nonzero_x  <- centered_x[centered_x != 0]
  n_nonzero  <- length(nonzero_x)
  
  if (n_nonzero == 0L) {
    return(0)
  }
  
  positive_count <- sum(nonzero_x > 0)
  
  ## Standardized signed sign statistic
  (positive_count - n_nonzero / 2) / sqrt(n_nonzero / 4)
}


# standardized Wilcoxon rank-sum statistic
get_wilcoxon_rank_sum_statistic <- function(x, y) {
  
  nx <- length(x)
  ny <- length(y)
  
  ranks <- rank(c(x, y))
  w_sum <- sum(ranks[seq_len(nx)])
  
  expected_w <- nx * (nx + ny + 1) / 2
  var_w      <- nx * ny * (nx + ny + 1) / 12
  
  ## Standardized Wilcoxon rank-sum statistic
  (w_sum - expected_w) / sqrt(var_w)
}


## ============================================================
## 6. One-sample simulation functions
## ============================================================

simulate_one_sample_unconditional <- function(N, n, dist_name) {
  
  tstats  <- numeric(N)
  pvalues <- numeric(N)
  
  pb <- txtProgressBar(min = 0, max = N, style = 3)
  
  for (sim in seq_len(N)) {
    
    x <- generate_data(n, dist_name)
    
    tstats[sim]  <- get_one_sample_tstat(x)
    pvalues[sim] <- get_one_sample_t_pvalue(x)
    
    setTxtProgressBar(pb, sim)
  }
  
  close(pb)
  
  list(
    tstats  = tstats,
    pvalues = pvalues
  )
}


simulate_one_sample_conditional <- function(N, n, dist_name, alpha = 0.05) {
  
  tstats  <- numeric(N)
  pvalues <- numeric(N)
  
  count <- 0L
  pb    <- txtProgressBar(min = 0, max = N, style = 3)
  
  while (count < N) {
    
    x <- generate_data(n, dist_name)
    
    if (passes_shapiro(x, alpha)) {
      
      count <- count + 1L
      
      tstats[count]  <- get_one_sample_tstat(x)
      pvalues[count] <- get_one_sample_t_pvalue(x)
      
      setTxtProgressBar(pb, count)
    }
  }
  
  close(pb)
  
  list(
    tstats  = tstats,
    pvalues = pvalues
  )
}


simulate_one_sample_adaptive <- function(N, n, dist_name, alpha = 0.05) {
  
  tstats      <- numeric(N)
  pvalues     <- numeric(N)
  test_branch <- character(N)
  
  pb <- txtProgressBar(min = 0, max = N, style = 3)
  
  for (sim in seq_len(N)) {
    
    x <- generate_data(n, dist_name, center_by = "median")
    
    if (passes_shapiro(x, alpha)) {
      
      tstats[sim]      <- get_one_sample_tstat(x)
      pvalues[sim]     <- get_one_sample_t_pvalue(x)
      test_branch[sim] <- "t-test"
      
    } else {
      
      tstats[sim]      <- get_signed_sign_statistic(x, mu0 = 0)
      pvalues[sim]     <- randomized_sign_test(x, mu0 = 0)
      test_branch[sim] <- "sign test"
    }
    
    setTxtProgressBar(pb, sim)
  }
  
  close(pb)
  
  list(
    tstats      = tstats,
    pvalues     = pvalues,
    test_branch = test_branch
  )
}


## ============================================================
## 7. Two-sample simulation functions
## ============================================================

simulate_two_sample_unconditional <- function(N, n, dist_name) {
  
  tstats  <- numeric(N)
  pvalues <- numeric(N)
  
  pb <- txtProgressBar(min = 0, max = N, style = 3)
  
  for (sim in seq_len(N)) {
    
    x <- generate_data(n, dist_name)
    y <- generate_data(n, dist_name)
    
    tstats[sim]  <- get_two_sample_tstat(x, y)
    pvalues[sim] <- get_two_sample_t_pvalue(x, y)
    
    setTxtProgressBar(pb, sim)
  }
  
  close(pb)
  
  list(
    tstats  = tstats,
    pvalues = pvalues
  )
}


simulate_two_sample_conditional <- function(N, n, dist_name, alpha = 0.05) {
  
  tstats  <- numeric(N)
  pvalues <- numeric(N)
  
  count <- 0L
  pb    <- txtProgressBar(min = 0, max = N, style = 3)
  
  while (count < N) {
    
    x <- generate_data(n, dist_name)
    y <- generate_data(n, dist_name)
    
    if (passes_shapiro(x, alpha) && passes_shapiro(y, alpha)) {
      
      count <- count + 1L
      
      tstats[count]  <- get_two_sample_tstat(x, y)
      pvalues[count] <- get_two_sample_t_pvalue(x, y)
      
      setTxtProgressBar(pb, count)
    }
  }
  
  close(pb)
  
  list(
    tstats  = tstats,
    pvalues = pvalues
  )
}


simulate_two_sample_adaptive <- function(N, n, dist_name, alpha = 0.05) {
  
  tstats      <- numeric(N)
  pvalues     <- numeric(N)
  test_branch <- character(N)
  
  pb <- txtProgressBar(min = 0, max = N, style = 3)
  
  for (sim in seq_len(N)) {
    
    x <- generate_data(n, dist_name)
    y <- generate_data(n, dist_name)
    
    x_passes <- passes_shapiro(x, alpha)
    y_passes <- passes_shapiro(y, alpha)
    
    if (x_passes && y_passes) {
      
      tstats[sim]      <- get_two_sample_tstat(x, y)
      pvalues[sim]     <- get_two_sample_t_pvalue(x, y)
      test_branch[sim] <- "t-test"
      
    } else {
      
      tstats[sim]      <- get_wilcoxon_rank_sum_statistic(x, y)
      pvalues[sim]     <- get_two_sample_wilcox_pvalue(x, y)
      test_branch[sim] <- "Wilcoxon"
    }
    
    setTxtProgressBar(pb, sim)
  }
  
  close(pb)
  
  list(
    tstats      = tstats,
    pvalues     = pvalues,
    test_branch = test_branch
  )
}


## ============================================================
## 8. Run simulations for all distributions
## ============================================================

run_all_simulations <- function(distributions, N, n, alpha) {
  
  results <- vector("list", length(distributions))
  names(results) <- distributions
  
  for (dist_name in distributions) {
    
    cat("\n============================================================\n")
    cat("Distribution:", dist_name, "\n")
    cat("============================================================\n")
    
    cat("\nOne-sample: unconditional t-test\n")
    one_uncond <- simulate_one_sample_unconditional(
      N         = N,
      n         = n,
      dist_name = dist_name
    )
    
    cat("\nOne-sample: conditional t-test\n")
    one_cond <- simulate_one_sample_conditional(
      N         = N,
      n         = n,
      dist_name = dist_name,
      alpha     = alpha
    )
    
    cat("\nOne-sample: adaptive test\n")
    one_adapt <- simulate_one_sample_adaptive(
      N         = N,
      n         = n,
      dist_name = dist_name,
      alpha     = alpha
    )
    
    cat("\nTwo-sample: unconditional t-test\n")
    two_uncond <- simulate_two_sample_unconditional(
      N         = N,
      n         = n,
      dist_name = dist_name
    )
    
    cat("\nTwo-sample: conditional t-test\n")
    two_cond <- simulate_two_sample_conditional(
      N         = N,
      n         = n,
      dist_name = dist_name,
      alpha     = alpha
    )
    
    cat("\nTwo-sample: adaptive test\n")
    two_adapt <- simulate_two_sample_adaptive(
      N         = N,
      n         = n,
      dist_name = dist_name,
      alpha     = alpha
    )
    
    results[[dist_name]] <- list(
      
      one_sample = list(
        unconditional = one_uncond,
        conditional   = one_cond,
        adaptive      = one_adapt
      ),
      
      two_sample = list(
        unconditional = two_uncond,
        conditional   = two_cond,
        adaptive      = two_adapt
      )
    )
  }
  
  results
}


simulation_results <- run_all_simulations(
  distributions = distributions,
  N             = N,
  n             = n,
  alpha         = alpha
)


## ============================================================
## 9. Tail-area Type I error summary
## ============================================================

get_critical_values <- function(df, alpha = 0.05) {
  qt(c(alpha / 2, 1 - alpha / 2), df = df)
}


estimate_tail_error <- function(statistics, critical_values) {
  
  lower_crit <- critical_values[1]
  upper_crit <- critical_values[2]
  
  left_tail  <- mean(statistics < lower_crit)
  right_tail <- mean(statistics > upper_crit)
  
  list(
    left_tail   = left_tail,
    right_tail  = right_tail,
    total_error = left_tail + right_tail
  )
}


create_tail_error_summary <- function(simulation_results,
                                      distributions,
                                      alpha = 0.05) {
  
  crit_one <- get_critical_values(df = df_one_sample, alpha = alpha)
  crit_two <- get_critical_values(df = df_two_sample, alpha = alpha)
  
  rows   <- list()
  row_id <- 1L
  
  procedure_labels <- c(
    unconditional = "Unconditional t-test",
    conditional   = "Conditional t-test",
    adaptive      = "Adaptive test"
  )
  
  for (dist_name in distributions) {
    
    for (setting_key in c("one_sample", "two_sample")) {
      
      setting_label <- ifelse(
        setting_key == "one_sample",
        "One-sample",
        "Two-sample"
      )
      
      critical_values <- if (setting_key == "one_sample") {
        crit_one
      } else {
        crit_two
      }
      
      for (procedure_key in names(procedure_labels)) {
        
        statistics <- simulation_results[[dist_name]][[setting_key]][[procedure_key]]$tstats
        
        error_parts <- estimate_tail_error(
          statistics      = statistics,
          critical_values = critical_values
        )
        
        rows[[row_id]] <- data.frame(
          Distribution = dist_name,
          Setting      = setting_label,
          Procedure    = procedure_labels[procedure_key],
          Left_Tail    = error_parts$left_tail,
          Right_Tail   = error_parts$right_tail,
          Type_I_Error = error_parts$total_error,
          stringsAsFactors = FALSE
        )
        
        row_id <- row_id + 1L
      }
    }
  }
  
  summary_table <- do.call(rbind, rows)
  
  summary_table$Left_Tail    <- round(summary_table$Left_Tail, 4)
  summary_table$Right_Tail   <- round(summary_table$Right_Tail, 4)
  summary_table$Type_I_Error <- round(summary_table$Type_I_Error, 4)
  
  summary_table
}


tail_error_summary <- create_tail_error_summary(
  simulation_results = simulation_results,
  distributions      = distributions,
  alpha              = alpha
)


cat("\n============================================================\n")
cat("Tail-area Type I Error Summary\n")
cat("============================================================\n")
print(tail_error_summary, row.names = FALSE)


write.csv(
  tail_error_summary,
  file      = "tail_area_type1_summary.csv",
  row.names = FALSE
)


## ============================================================
## 10. Plotting helper functions
## ============================================================

get_theoretical_t_density <- function(df,
                                      xlim = c(-5, 5),
                                      length_out = 1000) {
  
  x_values <- seq(xlim[1], xlim[2], length.out = length_out)
  
  list(
    x = x_values,
    y = dt(x_values, df = df)
  )
}


plot_density_panel <- function(unconditional_tstats,
                               conditional_tstats,
                               theoretical_density,
                               critical_values,
                               panel_title,
                               xlab = "Test statistic",
                               ylab = "Density") {
  
  dens_uncond <- density(unconditional_tstats)
  dens_cond   <- density(conditional_tstats)
  
  ymax <- max(
    dens_uncond$y,
    dens_cond$y,
    theoretical_density$y,
    na.rm = TRUE
  )
  
  plot(
    theoretical_density$x,
    theoretical_density$y,
    type = "l",
    col  = "red",
    lwd  = 2.5,
    ylim = c(0, ymax),
    xlab = xlab,
    ylab = ylab,
    main = panel_title
  )
  
  lines(dens_uncond, col = "orange",    lwd = 2.2)
  lines(dens_cond,   col = "steelblue", lwd = 2.2, lty = 2)
  
  abline(v = critical_values, col = "grey30", lty = 3, lwd = 2)
}


add_density_legend <- function() {
  
  par(
    fig = c(0, 1, 0, 1),
    mar = c(0, 0, 0, 0),
    oma = c(0, 0, 0, 0),
    new = TRUE
  )
  
  plot.new()
  
  legend(
    "bottom",
    legend = c(
      "Theoretical t",
      "Unconditional t-test",
      "Conditional t-test",
      "Critical values"
    ),
    col    = c("red", "orange", "steelblue", "grey30"),
    lty    = c(1, 1, 2, 3),
    lwd    = c(3, 3, 3, 2),
    horiz  = TRUE,
    bty    = "n",
    cex    = 1.35
  )
}


## ============================================================
## 11. One-sample density plots for all distributions
## ============================================================

plot_one_sample_density_panels <- function(simulation_results,
                                           distributions,
                                           output_file = "tstat_density_one_sample_all_distributions.pdf") {
  
  critical_values <- get_critical_values(df = df_one_sample, alpha = alpha)
  t_density       <- get_theoretical_t_density(df = df_one_sample)
  
  pdf(output_file, width = 12, height = 9)
  
  par(
    mfrow    = c(3, 2),
    mar      = c(4.8, 5.0, 3.0, 1.0),
    oma      = c(2.5, 0, 3.5, 0),
    cex.lab  = 1.35,
    cex.axis = 1.15,
    cex.main = 1.20,
    font.lab = 2,
    mgp      = c(2.6, 0.8, 0)
  )
  
  for (dist_name in distributions) {
    
    one_sample <- simulation_results[[dist_name]]$one_sample
    
    plot_density_panel(
      unconditional_tstats = one_sample$unconditional$tstats,
      conditional_tstats   = one_sample$conditional$tstats,
      theoretical_density  = t_density,
      critical_values      = critical_values,
      panel_title          = dist_name
    )
  }
  
  mtext(
    "One-Sample Test Statistic Distributions",
    outer = TRUE,
    cex   = 1.8,
    line  = 1.3,
    font  = 2
  )
  
  add_density_legend()
  
  dev.off()
  
  cat("Saved one-sample density plot:", output_file, "\n")
}


## ============================================================
## 12. Two-sample density plots for all distributions
## ============================================================

plot_two_sample_density_panels <- function(simulation_results,
                                           distributions,
                                           output_file = "tstat_density_two_sample_all_distributions.pdf") {
  
  critical_values <- get_critical_values(df = df_two_sample, alpha = alpha)
  t_density       <- get_theoretical_t_density(df = df_two_sample)
  
  pdf(output_file, width = 12, height = 9)
  
  par(
    mfrow    = c(3, 2),
    mar      = c(4.8, 5.0, 3.0, 1.0),
    oma      = c(2.5, 0, 3.5, 0),
    cex.lab  = 1.35,
    cex.axis = 1.15,
    cex.main = 1.20,
    font.lab = 2,
    mgp      = c(2.6, 0.8, 0)
  )
  
  for (dist_name in distributions) {
    
    two_sample <- simulation_results[[dist_name]]$two_sample
    
    plot_density_panel(
      unconditional_tstats = two_sample$unconditional$tstats,
      conditional_tstats   = two_sample$conditional$tstats,
      theoretical_density  = t_density,
      critical_values      = critical_values,
      panel_title          = dist_name
    )
  }
  
  mtext(
    "Two-Sample Test Statistic Distributions",
    outer = TRUE,
    cex   = 1.8,
    line  = 1.3,
    font  = 2
  )
  
  add_density_legend()
  
  dev.off()
  
  cat("Saved two-sample density plot:", output_file, "\n")
}

## ============================================================
## 13. Generate plots
## ============================================================

plot_one_sample_density_panels(
  simulation_results = simulation_results,
  distributions      = distributions
)

plot_two_sample_density_panels(
  simulation_results = simulation_results,
  distributions      = distributions
)


## ============================================================
## 14. Save results
## ============================================================

save(
  simulation_results,
  tail_error_summary,
  file = "test_stat_results.RData"
)

save.image("test_stat_workspace.RData")

