## ============================================================
## Selective Inference in Normality Pre-testing
## Type I Error and Power Simulation Study
##
## Efficient  version
##
## Main features:
##   1. Keeps the regular t-test as a benchmark method.
##   2. Includes an oracle procedure.
##   3. Includes regular adaptive pre-testing.
##   4. Builds regular t, oracle, and adaptive from one
##      full-sample p-value data frame.
##   5. Splits split-sample simulations into:
##        - pretest-only simulations
##        - downstream-only simulations
##   6. Combines reusable pretest/downstream data frames to form:
##        - Split n:       0.5n pretest, 0.5n downstream
##        - Split 1.5n A:  n pretest,    0.5n downstream
##        - Split 1.5n B:  0.5n pretest, n downstream
##        - Split 2n:      n pretest,    n downstream
##   7. Saves all p-values in data frames:
##        - t-test p-value
##        - sign-test p-value
##        - Shapiro-Wilk p-value
##        - final selected p-value
##        - selected downstream test
##   8. Runs both regular and randomized sign-test versions.
##   9. Saves all plots as PDF only.
##  10. Uses gap::qqunif() for QQ plots.
##  11. Provides one main run function where the user specifies
##      simulation inputs.
## ============================================================


## ============================================================
## 0. Setup
## ============================================================

## Check that the gap package is available.
if (!requireNamespace("gap", quietly = TRUE)) {
  stop("Package 'gap' is required. Please install it using install.packages('gap').")
}

## Check that the parallel package is available.
if (!requireNamespace("parallel", quietly = TRUE)) {
  stop("Package 'parallel' is required.")
}

## Set the working directory before running this script, if needed.
#setwd("~/Desktop/OSU/Research/Pretest-Simulation/Selective Inference/Selective_Inf_exploration")


## ============================================================
## 1. Default distribution setup
## ============================================================

## Defines the distributions used in the simulation.
make_default_distribution_info <- function() {
  
  list(
    
    normal = list(
      mean   = 0,
      median = 0,
      sd     = 1,
      sample = function(n) rnorm(n, mean = 0, sd = 1)
    ),
    
    exponential = list(
      mean   = 1,
      median = log(2),
      sd     = 1,
      sample = function(n) rexp(n, rate = 1)
    ),
    
    lognormal = list(
      mean   = exp(0.5),
      median = 1,
      sd     = sqrt((exp(1) - 1) * exp(1)),
      sample = function(n) rlnorm(n, meanlog = 0, sdlog = 1)
    )
  )
}


## Defines display labels for the distributions.
make_default_distribution_labels <- function() {
  
  c(
    normal      = "Normal",
    exponential = "Exponential",
    lognormal   = "Lognormal"
  )
}


## ============================================================
## 2. Method labels and plotting symbols
## ============================================================

## Stores method names, labels, colors, and plotting symbols.
make_method_settings <- function() {
  
  method_keys <- c(
    "regular_t",
    "oracle",
    "regular_adaptive",
    "adaptive_split_n",
    "adaptive_split_15n_A",
    "adaptive_split_15n_B",
    "adaptive_split_2n"
  )
  
  method_labels <- c(
    regular_t             = "Regular t",
    oracle                = "Oracle",
    regular_adaptive      = "Regular Adaptive",
    adaptive_split_n      = "Adaptive Split n: 1:1",
    adaptive_split_15n_A  = "Adaptive Split 1.5n: 2:1",
    adaptive_split_15n_B  = "Adaptive Split 1.5n: 1:2",
    adaptive_split_2n     = "Adaptive Split 2n: 1:1"
  )
  
  method_colors <- c(
    regular_t             = "grey40",
    oracle                = "lightskyblue",
    regular_adaptive      = "orange",
    adaptive_split_n      = "lightgreen",
    adaptive_split_15n_A  = "#9370DB",
    adaptive_split_15n_B  = "#20B2AA",
    adaptive_split_2n     = "red"
  )
  
  method_pch <- c(
    regular_t             = 19,
    oracle                = 16,
    regular_adaptive      = 17,
    adaptive_split_n      = 15,
    adaptive_split_15n_A  = 23,
    adaptive_split_15n_B  = 24,
    adaptive_split_2n     = 18
  )
  
  list(
    method_keys   = method_keys,
    method_labels = method_labels,
    method_colors = method_colors,
    method_pch    = method_pch
  )
}


## ============================================================
## 3. General helper functions
## ============================================================

## Creates the output directory if it does not already exist.
ensure_output_dir <- function(output_dir) {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  if (!dir.exists(output_dir)) {
    stop("Output directory could not be created: ", output_dir)
  }
  
  invisible(normalizePath(output_dir, mustWork = TRUE))
}


## Generates data from the selected distribution and applies the effect shift.
generate_data <- function(sample_size,
                          distribution_name,
                          effect_size,
                          distribution_info) {
  
  dist_info <- distribution_info[[distribution_name]]
  
  x <- dist_info$sample(sample_size)
  
  ## Power is induced by shifting the distribution by effect_size SDs.
  x + effect_size * dist_info$sd
}


## Centers data at the population median and scales by the population SD.
center_by_median <- function(x,
                             distribution_name,
                             distribution_info) {
  
  dist_info <- distribution_info[[distribution_name]]
  
  (x - dist_info$median) / dist_info$sd
}


## Returns the pretest and downstream sample sizes for all split procedures.
get_split_sizes <- function(sample_size) {
  
  n_half <- floor(sample_size / 2)
  
  list(
    
    adaptive_split_n = c(
      n_pretest    = n_half,
      n_downstream = sample_size - n_half
    ),
    
    adaptive_split_15n_A = c(
      n_pretest    = sample_size,
      n_downstream = ceiling(sample_size / 2)
    ),
    
    adaptive_split_15n_B = c(
      n_pretest    = ceiling(sample_size / 2),
      n_downstream = sample_size
    ),
    
    adaptive_split_2n = c(
      n_pretest    = sample_size,
      n_downstream = sample_size
    )
  )
}


## ============================================================
## 4. Test functions
## ============================================================

## Computes the one-sample t-test p-value for H0: location = 0.
run_ttest <- function(x) {
  
  t.test(x, mu = 0)$p.value
}


## Computes the ordinary two-sided exact sign-test p-value.
regular_sign_test <- function(x, mu0 = 0) {
  
  centered_x <- x - mu0
  nonzero_n  <- sum(centered_x != 0)
  
  if (nonzero_n == 0L) {
    return(1)
  }
  
  positive_count <- sum(centered_x > 0)
  
  lower_tail <- pbinom(
    positive_count,
    nonzero_n,
    prob       = 0.5,
    lower.tail = TRUE
  )
  
  upper_tail <- pbinom(
    positive_count - 1L,
    nonzero_n,
    prob       = 0.5,
    lower.tail = FALSE
  )
  
  min(1, 2 * min(lower_tail, upper_tail))
}


## Computes the randomized two-sided sign-test p-value.
randomized_sign_test <- function(x, mu0 = 0) {
  
  centered_x <- x - mu0
  nonzero_n  <- sum(centered_x != 0)
  
  if (nonzero_n == 0L) {
    return(runif(1))
  }
  
  positive_count <- sum(centered_x > 0)
  random_u       <- runif(1)
  
  probability_mass <- dbinom(
    positive_count,
    nonzero_n,
    prob = 0.5
  )
  
  left_tail <- pbinom(
    positive_count - 1L,
    nonzero_n,
    prob       = 0.5,
    lower.tail = TRUE
  ) + random_u * probability_mass
  
  right_tail <- pbinom(
    positive_count,
    nonzero_n,
    prob       = 0.5,
    lower.tail = FALSE
  ) + random_u * probability_mass
  
  min(1, 2 * min(left_tail, right_tail))
}


## Chooses either the regular or randomized sign test.
run_sign_test <- function(x, sign_test_type = c("regular", "randomized")) {
  
  sign_test_type <- match.arg(sign_test_type)
  
  if (sign_test_type == "regular") {
    regular_sign_test(x)
  } else {
    randomized_sign_test(x)
  }
}


## Computes the Shapiro-Wilk normality-test p-value.
get_shapiro_pvalue <- function(x) {
  
  shapiro.test(x)$p.value
}


## ============================================================
## 5. P-value data-frame helpers
## ============================================================

## Creates an empty p-value data frame for method-level results.
make_empty_pvalue_frame <- function(num_simulations) {
  
  data.frame(
    t_pvalue       = rep(NA_real_, num_simulations),
    sign_pvalue    = rep(NA_real_, num_simulations),
    shapiro_pvalue = rep(NA_real_, num_simulations),
    final_pvalue   = rep(NA_real_, num_simulations),
    selected_test  = rep(NA_character_, num_simulations),
    stringsAsFactors = FALSE
  )
}


## ============================================================
## 6. Component p-value simulations
## ============================================================

## Simulates full-sample p-values used by regular t, oracle, and regular adaptive.
simulate_full_sample_pvalues <- function(num_simulations,
                                         sample_size,
                                         distribution_name,
                                         effect_size,
                                         sign_test_type,
                                         distribution_info) {
  
  output <- make_empty_pvalue_frame(num_simulations)
  
  for (sim in seq_len(num_simulations)) {
    
    ## Generate and center one full sample.
    x <- generate_data(
      sample_size       = sample_size,
      distribution_name = distribution_name,
      effect_size       = effect_size,
      distribution_info = distribution_info
    )
    
    xc <- center_by_median(
      x                 = x,
      distribution_name = distribution_name,
      distribution_info = distribution_info
    )
    
    ## Store component p-values.
    output$t_pvalue[sim]       <- run_ttest(xc)
    output$sign_pvalue[sim]    <- run_sign_test(xc, sign_test_type)
    output$shapiro_pvalue[sim] <- get_shapiro_pvalue(xc)
  }
  
  output
}


## Simulates only Shapiro-Wilk p-values for a pretest sample size.
simulate_pretest_pvalues <- function(num_simulations,
                                     n_pretest,
                                     distribution_name,
                                     effect_size,
                                     distribution_info) {
  
  output <- data.frame(
    shapiro_pvalue = rep(NA_real_, num_simulations)
  )
  
  for (sim in seq_len(num_simulations)) {
    
    ## Generate and center one pretest sample.
    x <- generate_data(
      sample_size       = n_pretest,
      distribution_name = distribution_name,
      effect_size       = effect_size,
      distribution_info = distribution_info
    )
    
    xc <- center_by_median(
      x                 = x,
      distribution_name = distribution_name,
      distribution_info = distribution_info
    )
    
    ## Store the pretest p-value.
    output$shapiro_pvalue[sim] <- get_shapiro_pvalue(xc)
  }
  
  output
}


## Simulates t-test and sign-test p-values for a downstream sample size.
simulate_downstream_pvalues <- function(num_simulations,
                                        n_downstream,
                                        distribution_name,
                                        effect_size,
                                        sign_test_type,
                                        distribution_info) {
  
  output <- data.frame(
    t_pvalue    = rep(NA_real_, num_simulations),
    sign_pvalue = rep(NA_real_, num_simulations),
    stringsAsFactors = FALSE
  )
  
  for (sim in seq_len(num_simulations)) {
    
    ## Generate and center one downstream sample.
    x <- generate_data(
      sample_size       = n_downstream,
      distribution_name = distribution_name,
      effect_size       = effect_size,
      distribution_info = distribution_info
    )
    
    xc <- center_by_median(
      x                 = x,
      distribution_name = distribution_name,
      distribution_info = distribution_info
    )
    
    ## Store downstream p-values.
    output$t_pvalue[sim]    <- run_ttest(xc)
    output$sign_pvalue[sim] <- run_sign_test(xc, sign_test_type)
  }
  
  output
}


## ============================================================
## 7. Build methods from component p-values
## ============================================================

## Builds regular t, oracle, and regular adaptive from one full-sample p-value data frame.
build_full_sample_methods <- function(full_pvalues,
                                      distribution_name,
                                      pretest_alpha) {
  
  ## Regular t always uses the t-test p-value.
  regular_t <- full_pvalues
  regular_t$final_pvalue  <- regular_t$t_pvalue
  regular_t$selected_test <- "t"
  
  ## Oracle uses t for normal data and sign for non-normal data.
  oracle <- full_pvalues
  
  if (distribution_name == "normal") {
    oracle$final_pvalue  <- oracle$t_pvalue
    oracle$selected_test <- "t"
  } else {
    oracle$final_pvalue  <- oracle$sign_pvalue
    oracle$selected_test <- "sign"
  }
  
  ## Regular adaptive uses the same sample for pretest and downstream inference.
  regular_adaptive <- full_pvalues
  
  use_ttest <- regular_adaptive$shapiro_pvalue > pretest_alpha
  
  regular_adaptive$final_pvalue <- ifelse(
    use_ttest,
    regular_adaptive$t_pvalue,
    regular_adaptive$sign_pvalue
  )
  
  regular_adaptive$selected_test <- ifelse(
    use_ttest,
    "t",
    "sign"
  )
  
  list(
    regular_t        = regular_t,
    oracle           = oracle,
    regular_adaptive = regular_adaptive
  )
}


## Builds one adaptive split method from pretest and downstream p-value data frames.
build_split_method <- function(pretest_pvalues,
                               downstream_pvalues,
                               pretest_alpha) {
  
  use_ttest <- pretest_pvalues$shapiro_pvalue > pretest_alpha
  
  data.frame(
    t_pvalue       = downstream_pvalues$t_pvalue,
    sign_pvalue    = downstream_pvalues$sign_pvalue,
    shapiro_pvalue = pretest_pvalues$shapiro_pvalue,
    final_pvalue   = ifelse(
      use_ttest,
      downstream_pvalues$t_pvalue,
      downstream_pvalues$sign_pvalue
    ),
    selected_test  = ifelse(use_ttest, "t", "sign"),
    stringsAsFactors = FALSE
  )
}


## ============================================================
## 8. Simulation engine
## ============================================================

## Runs one distribution/sample-size task by simulating reusable component p-values.
run_one_simulation_task <- function(task,
                                    num_simulations,
                                    effect_size_value,
                                    pretest_alpha,
                                    sign_test_type,
                                    distribution_info) {
  
  distribution_name <- task$distribution_name
  sample_size       <- task$sample_size
  
  result_key  <- paste(distribution_name, sample_size, sep = "_")
  split_sizes <- get_split_sizes(sample_size)
  
  result <- list(
    distribution = distribution_name,
    sample_size  = sample_size,
    split_sizes  = split_sizes
  )
  
  ## Simulate full-sample component p-values once.
  full_pvalues <- simulate_full_sample_pvalues(
    num_simulations   = num_simulations,
    sample_size       = sample_size,
    distribution_name = distribution_name,
    effect_size       = effect_size_value,
    sign_test_type    = sign_test_type,
    distribution_info = distribution_info
  )
  
  ## Build full-sample methods outside the simulation loop.
  full_methods <- build_full_sample_methods(
    full_pvalues      = full_pvalues,
    distribution_name = distribution_name,
    pretest_alpha     = pretest_alpha
  )
  
  result$regular_t        <- full_methods$regular_t
  result$oracle           <- full_methods$oracle
  result$regular_adaptive <- full_methods$regular_adaptive
  
  ## Identify unique pretest and downstream sizes needed for the split methods.
  unique_pretest_sizes <- unique(vapply(
    split_sizes,
    function(z) z["n_pretest"],
    numeric(1)
  ))
  
  unique_downstream_sizes <- unique(vapply(
    split_sizes,
    function(z) z["n_downstream"],
    numeric(1)
  ))
  
  ## Simulate pretest p-values once for each required pretest size.
  pretest_bank <- list()
  
  for (n_pretest in unique_pretest_sizes) {
    
    pretest_bank[[as.character(n_pretest)]] <- simulate_pretest_pvalues(
      num_simulations   = num_simulations,
      n_pretest         = n_pretest,
      distribution_name = distribution_name,
      effect_size       = effect_size_value,
      distribution_info = distribution_info
    )
  }
  
  ## Simulate downstream p-values once for each required downstream size.
  downstream_bank <- list()
  
  for (n_downstream in unique_downstream_sizes) {
    
    downstream_bank[[as.character(n_downstream)]] <- simulate_downstream_pvalues(
      num_simulations   = num_simulations,
      n_downstream      = n_downstream,
      distribution_name = distribution_name,
      effect_size       = effect_size_value,
      sign_test_type    = sign_test_type,
      distribution_info = distribution_info
    )
  }
  
  ## Build each split method by combining reusable pretest and downstream data frames.
  for (split_method in names(split_sizes)) {
    
    n_pretest    <- split_sizes[[split_method]]["n_pretest"]
    n_downstream <- split_sizes[[split_method]]["n_downstream"]
    
    result[[split_method]] <- build_split_method(
      pretest_pvalues    = pretest_bank[[as.character(n_pretest)]],
      downstream_pvalues = downstream_bank[[as.character(n_downstream)]],
      pretest_alpha      = pretest_alpha
    )
  }
  
  list(
    result_key = result_key,
    result     = result
  )
}


## Runs one full Type I error or power simulation set.
run_simulation_set <- function(num_simulations,
                               distributions,
                               sample_sizes,
                               effect_size_value,
                               pretest_alpha,
                               sign_test_type,
                               method_keys,
                               distribution_info,
                               n_cores = max(1, parallel::detectCores() - 1)) {
  
  ## Create independent distribution/sample-size tasks.
  task_grid <- expand.grid(
    distribution_name = distributions,
    sample_size       = sample_sizes,
    stringsAsFactors  = FALSE
  )
  
  tasks <- split(task_grid, seq_len(nrow(task_grid)))
  
  cat(
    "\nRunning simulation set",
    "| effect size:", effect_size_value,
    "| sign test:", sign_test_type,
    "| cores:", n_cores,
    "\n"
  )
  
  ## Run sequentially or in parallel.
  if (n_cores <= 1) {
    
    task_results <- lapply(tasks, function(task) {
      run_one_simulation_task(
        task              = task,
        num_simulations   = num_simulations,
        effect_size_value = effect_size_value,
        pretest_alpha     = pretest_alpha,
        sign_test_type    = sign_test_type,
        distribution_info = distribution_info
      )
    })
    
  } else {
    
    ## Do not request more cores than there are tasks.
    n_cores <- min(n_cores, length(tasks))
    
    task_results <- parallel::mclapply(
      tasks,
      function(task) {
        run_one_simulation_task(
          task              = task,
          num_simulations   = num_simulations,
          effect_size_value = effect_size_value,
          pretest_alpha     = pretest_alpha,
          sign_test_type    = sign_test_type,
          distribution_info = distribution_info
        )
      },
      mc.cores    = n_cores,
      mc.set.seed = TRUE
    )
  }
  
  ## Put task results back into the expected named list format.
  simulation_output <- list()
  
  for (task_result in task_results) {
    simulation_output[[task_result$result_key]] <- task_result$result
  }
  
  simulation_output
}


## ============================================================
## 9. Summary table
## ============================================================

## Creates the summary table of Type I error, power, and routing frequencies.
create_error_power_table <- function(type1_results,
                                     power_results,
                                     distributions,
                                     sample_sizes,
                                     method_keys,
                                     method_labels,
                                     distribution_labels,
                                     alpha_level,
                                     effect_size,
                                     sign_test_type) {
  
  table_rows <- list()
  row_id     <- 1L
  
  for (distribution_name in distributions) {
    
    for (sample_size in sample_sizes) {
      
      result_key <- paste(distribution_name, sample_size, sep = "_")
      
      for (method_key in method_keys) {
        
        type1_df <- type1_results[[result_key]][[method_key]]
        power_df <- power_results[[result_key]][[method_key]]
        
        table_rows[[row_id]] <- data.frame(
          Sign_Test_Type        = sign_test_type,
          Distribution          = distribution_labels[distribution_name],
          Sample_Size           = sample_size,
          Method_Key            = method_key,
          Method                = method_labels[method_key],
          Type_I_Error          = mean(type1_df$final_pvalue < alpha_level, na.rm = TRUE),
          Power                 = mean(power_df$final_pvalue < alpha_level, na.rm = TRUE),
          Prop_T_Selected_Type1 = mean(type1_df$selected_test == "t", na.rm = TRUE),
          Prop_T_Selected_Power = mean(power_df$selected_test == "t", na.rm = TRUE),
          Effect_Size           = effect_size,
          stringsAsFactors      = FALSE
        )
        
        row_id <- row_id + 1L
      }
    }
  }
  
  summary_table <- do.call(rbind, table_rows)
  
  summary_table$Type_I_Error <- round(summary_table$Type_I_Error, 4)
  summary_table$Power        <- round(summary_table$Power, 4)
  
  summary_table$Prop_T_Selected_Type1 <- round(
    summary_table$Prop_T_Selected_Type1,
    4
  )
  
  summary_table$Prop_T_Selected_Power <- round(
    summary_table$Prop_T_Selected_Power,
    4
  )
  
  summary_table
}


## ============================================================
## 10. Save p-value data frames as CSV files
## ============================================================

## Saves p-value data frames as CSV files when requested.
save_pvalue_dataframes_csv <- function(results,
                                       result_type,
                                       sign_test_type,
                                       method_keys,
                                       output_dir) {
  
  pvalue_dir <- file.path(output_dir, "pvalue_dataframes")
  ensure_output_dir(pvalue_dir)
  
  for (result_key in names(results)) {
    
    result_object <- results[[result_key]]
    
    for (method_key in method_keys) {
      
      output_file <- file.path(
        pvalue_dir,
        paste(
          result_type,
          sign_test_type,
          result_key,
          method_key,
          "pvalues.csv",
          sep = "_"
        )
      )
      
      write.csv(
        result_object[[method_key]],
        file      = output_file,
        row.names = FALSE
      )
    }
  }
}


## ============================================================
## 11. Type I error and power plots
## ============================================================

## Creates Type I error and power bar plots for one distribution.
plot_error_power_by_distribution <- function(distribution_name,
                                             error_power_table,
                                             sample_sizes,
                                             method_keys,
                                             method_labels,
                                             method_colors,
                                             distribution_labels,
                                             alpha_level,
                                             effect_size,
                                             sign_test_type,
                                             output_dir) {
  
  output_file <- file.path(
    output_dir,
    paste0("type1_power_", distribution_name, "_", sign_test_type, ".pdf")
  )
  
  pdf(
    file   = output_file,
    width  = 5.8 * length(sample_sizes),
    height = 10.0
  )
  
  on.exit(dev.off(), add = TRUE)
  
  layout(
    matrix(
      seq_len(2 * length(sample_sizes)),
      nrow  = 2,
      ncol  = length(sample_sizes),
      byrow = TRUE
    )
  )
  
  par(
    oma = c(2.0, 2.0, 6.0, 0.5),
    mar = c(8.8, 5.3, 4.0, 1.0)
  )
  
  metric_names  <- c("Type_I_Error", "Power")
  metric_titles <- c("Empirical Type I Error", "Empirical Power")
  
  type1_values <- error_power_table[
    error_power_table$Distribution == distribution_labels[distribution_name],
    "Type_I_Error"
  ]
  
  type1_ymax <- max(type1_values, na.rm = TRUE) * 1.20
  type1_ymax <- max(type1_ymax, alpha_level * 1.5)
  
  for (metric_id in seq_along(metric_names)) {
    
    metric_name  <- metric_names[metric_id]
    metric_title <- metric_titles[metric_id]
    
    for (sample_size in sample_sizes) {
      
      plot_data <- error_power_table[
        error_power_table$Distribution == distribution_labels[distribution_name] &
          error_power_table$Sample_Size == sample_size,
      ]
      
      plot_data <- plot_data[
        match(method_keys, plot_data$Method_Key),
      ]
      
      values <- plot_data[[metric_name]]
      
      y_limit <- if (metric_name == "Type_I_Error") {
        c(0, type1_ymax)
      } else {
        c(0, 1)
      }
      
      bp <- barplot(
        height    = values,
        col       = method_colors[method_keys],
        border    = "grey30",
        ylim      = y_limit,
        xaxt      = "n",
        ylab      = metric_title,
        main      = paste0("n = ", sample_size),
        cex.axis  = 1.50,
        cex.lab   = 1.50,
        cex.main  = 1.80,
        font.lab  = 2,
        font.main = 2
      )
      
      if (metric_name == "Type_I_Error") {
        
        abline(h = alpha_level, lty = 2, lwd = 2, col = "grey30")
        
        text(
          x      = 0.8,
          y      = alpha_level + type1_ymax * 0.05,
          labels = paste0("alpha = ", alpha_level),
          adj    = 0,
          cex    = 1.25
        )
      }
      
      text(
        x      = bp,
        y      = par("usr")[3] - 0.045 * diff(par("usr")[3:4]),
        labels = method_labels[method_keys],
        srt    = 45,
        adj    = 1,
        xpd    = TRUE,
        cex    = 1.05
      )
    }
  }
  
  mtext(
    paste0(
      "Type I Error and Power: ",
      distribution_labels[distribution_name],
      " Distribution"
    ),
    side  = 3,
    outer = TRUE,
    line  = 3.7,
    cex   = 2.10,
    font  = 2
  )
  
  mtext(
    paste0(
      "Sign-test version: ",
      sign_test_type,
      " | Power alternative: location shift = ",
      effect_size,
      " SD"
    ),
    side  = 3,
    outer = TRUE,
    line  = 2.0,
    cex   = 1.50
  )
  
  cat("Saved plot:", output_file, "\n")
}


## ============================================================
## 12. Uniform QQ plots
## ============================================================

## Creates log-scale uniform QQ plots using gap::qqunif().
plot_qq_pvalues <- function(type1_results,
                            distributions,
                            sample_sizes,
                            method_keys,
                            method_labels,
                            method_colors,
                            method_pch,
                            distribution_labels,
                            sign_test_type,
                            output_dir,
                            max_points = 4000) {
  
  output_file <- file.path(
    output_dir,
    paste0("uniform_qq_plots_", sign_test_type, ".pdf")
  )
  
  n_dist   <- length(distributions)
  n_size   <- length(sample_sizes)
  n_panels <- n_dist * n_size
  
  layout_matrix <- rbind(
    matrix(seq_len(n_panels), nrow = n_dist, ncol = n_size, byrow = TRUE),
    rep(n_panels + 1, n_size)
  )
  
  prepare_pvalues_for_plot <- function(pvalues, max_points = 4000) {
    
    pvalues <- pvalues[is.finite(pvalues)]
    pvalues <- pvalues[!is.na(pvalues)]
    
    pvalues <- pmax(pvalues, .Machine$double.xmin)
    pvalues <- pmin(pvalues, 1)
    
    pvalues <- sort(pvalues)
    n <- length(pvalues)
    
    if (n > max_points) {
      keep <- unique(round(seq(1, n, length.out = max_points)))
      pvalues <- pvalues[keep]
    }
    
    pvalues
  }
  
  get_qqunif_coordinates <- function(pvalues) {
    
    pvalues <- sort(pvalues)
    n       <- length(pvalues)
    
    expected <- (seq_len(n) - 0.5) / n
    
    data.frame(
      expected_log = -log10(expected),
      observed_log = -log10(pvalues)
    )
  }
  
  pdf(
    file   = output_file,
    width  = 6.8 * n_size,
    height = 6.3 * n_dist + 2.1
  )
  
  on.exit(dev.off(), add = TRUE)
  
  layout(
    layout_matrix,
    heights = c(rep(1, n_dist), 0.25)
  )
  
  par(
    oma = c(1.0, 2.0, 6.4, 1.0),
    mar = c(6.3, 6.7, 4.4, 1.8),
    xpd = FALSE
  )
  
  for (distribution_name in distributions) {
    
    for (sample_size in sample_sizes) {
      
      result_key <- paste(distribution_name, sample_size, sep = "_")
      result     <- type1_results[[result_key]]
      
      first_method <- method_keys[1]
      
      first_pvalues <- prepare_pvalues_for_plot(
        result[[first_method]]$final_pvalue,
        max_points = max_points
      )
      
      ## Base QQ plot from gap::qqunif().
      gap::qqunif(
        first_pvalues,
        main = sprintf(
          "%s | n = %d",
          distribution_labels[distribution_name],
          sample_size
        ),
        col       = method_colors[first_method],
        pch       = method_pch[first_method],
        cex       = 1.05,
        lwd       = 1.4,
        ci        = TRUE,
        cex.main  = 2.50,
        cex.lab   = 2.00,
        cex.axis  = 1.80,
        font.main = 2,
        font.lab  = 2,
        las       = 1
      )
      
      grid(col = "grey88", lty = 1, lwd = 0.9)
      abline(0, 1, lty = 2, col = "grey30", lwd = 2.6)
      
      ## Re-plot first method so it remains visible above the grid.
      first_coords <- get_qqunif_coordinates(first_pvalues)
      
      points(
        first_coords$expected_log,
        first_coords$observed_log,
        col = method_colors[first_method],
        pch = method_pch[first_method],
        cex = 1.05,
        lwd = 1.4
      )
      
      ## Overlay remaining methods.
      for (method_key in method_keys[-1]) {
        
        method_pvalues <- prepare_pvalues_for_plot(
          result[[method_key]]$final_pvalue,
          max_points = max_points
        )
        
        method_coords <- get_qqunif_coordinates(method_pvalues)
        
        points(
          method_coords$expected_log,
          method_coords$observed_log,
          col = method_colors[method_key],
          pch = method_pch[method_key],
          cex = 1.05,
          lwd = 1.4
        )
      }
    }
  }
  
  ## Legend panel.
  par(mar = c(0, 0, 0, 0), xpd = TRUE)
  plot.new()
  
  legend(
    "center",
    legend    = method_labels[method_keys],
    col       = method_colors[method_keys],
    pch       = method_pch[method_keys],
    pt.cex    = 2.10,
    ncol      = min(7, length(method_keys)),
    bty       = "n",
    cex       = 1.80,
    text.font = 2,
    x.intersp = 0.90,
    y.intersp = 1.25
  )
  
  mtext(
    "Log-scale Uniform Q-Q Plots of Type I Error p-values",
    outer = TRUE,
    side  = 3,
    line  = 3.9,
    cex   = 2.55,
    font  = 2
  )
  
  mtext(
    paste0("Sign-test version: ", sign_test_type),
    outer = TRUE,
    side  = 3,
    line  = 2.0,
    cex   = 1.75
  )
  
  cat("Saved QQ plot:", output_file, "\n")
}


## ============================================================
## 13. Workflow for one sign-test version
## ============================================================

## Runs the complete workflow for one sign-test version.
run_one_sign_test_workflow <- function(sign_test_type,
                                       num_simulations,
                                       distributions,
                                       sample_sizes,
                                       effect_size,
                                       alpha_level,
                                       method_settings,
                                       distribution_info,
                                       distribution_labels,
                                       output_dir,
                                       save_pvalue_csv,
                                       max_qq_points,
                                       n_cores) {
  
  sign_test_type <- match.arg(sign_test_type, c("regular", "randomized"))
  
  method_keys   <- method_settings$method_keys
  method_labels <- method_settings$method_labels
  method_colors <- method_settings$method_colors
  method_pch    <- method_settings$method_pch
  
  sign_output_dir <- file.path(
    output_dir,
    paste0("sign_test_", sign_test_type)
  )
  
  ensure_output_dir(sign_output_dir)
  
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat("Running workflow with ", sign_test_type, " sign test\n", sep = "")
  cat(strrep("=", 70), "\n", sep = "")
  
  ## Run Type I error simulations.
  type1_results <- run_simulation_set(
    num_simulations   = num_simulations,
    distributions     = distributions,
    sample_sizes      = sample_sizes,
    effect_size_value = 0,
    pretest_alpha     = alpha_level,
    sign_test_type    = sign_test_type,
    method_keys       = method_keys,
    distribution_info = distribution_info,
    n_cores           = n_cores
  )
  
  ## Run power simulations.
  power_results <- run_simulation_set(
    num_simulations   = num_simulations,
    distributions     = distributions,
    sample_sizes      = sample_sizes,
    effect_size_value = effect_size,
    pretest_alpha     = alpha_level,
    sign_test_type    = sign_test_type,
    method_keys       = method_keys,
    distribution_info = distribution_info,
    n_cores           = n_cores
  )
  
  ## Create summary table.
  error_power_table <- create_error_power_table(
    type1_results       = type1_results,
    power_results       = power_results,
    distributions       = distributions,
    sample_sizes        = sample_sizes,
    method_keys         = method_keys,
    method_labels       = method_labels,
    distribution_labels = distribution_labels,
    alpha_level         = alpha_level,
    effect_size         = effect_size,
    sign_test_type      = sign_test_type
  )
  
  summary_file <- file.path(
    sign_output_dir,
    paste0("type1_error_and_power_summary_", sign_test_type, ".csv")
  )
  
  write.csv(
    error_power_table,
    file      = summary_file,
    row.names = FALSE
  )
  
  ## Optionally export p-value data frames as CSV files.
  if (isTRUE(save_pvalue_csv)) {
    
    save_pvalue_dataframes_csv(
      results        = type1_results,
      result_type    = "type1",
      sign_test_type = sign_test_type,
      method_keys    = method_keys,
      output_dir     = sign_output_dir
    )
    
    save_pvalue_dataframes_csv(
      results        = power_results,
      result_type    = "power",
      sign_test_type = sign_test_type,
      method_keys    = method_keys,
      output_dir     = sign_output_dir
    )
  }
  
  ## Create Type I error and power plots.
  for (distribution_name in distributions) {
    
    plot_error_power_by_distribution(
      distribution_name    = distribution_name,
      error_power_table    = error_power_table,
      sample_sizes         = sample_sizes,
      method_keys          = method_keys,
      method_labels        = method_labels,
      method_colors        = method_colors,
      distribution_labels  = distribution_labels,
      alpha_level          = alpha_level,
      effect_size          = effect_size,
      sign_test_type       = sign_test_type,
      output_dir           = sign_output_dir
    )
  }
  
  ## Create QQ plots for Type I error p-values.
  plot_qq_pvalues(
    type1_results       = type1_results,
    distributions       = distributions,
    sample_sizes        = sample_sizes,
    method_keys         = method_keys,
    method_labels       = method_labels,
    method_colors       = method_colors,
    method_pch          = method_pch,
    distribution_labels = distribution_labels,
    sign_test_type      = sign_test_type,
    output_dir          = sign_output_dir,
    max_points          = max_qq_points
  )
  
  ## Save sign-test-specific workspace.
  workspace_file <- file.path(
    sign_output_dir,
    paste0("type1_power_results_", sign_test_type, ".RData")
  )
  
  simulation_settings <- list(
    sign_test_type  = sign_test_type,
    num_simulations = num_simulations,
    distributions   = distributions,
    sample_sizes    = sample_sizes,
    effect_size     = effect_size,
    alpha_level     = alpha_level
  )
  
  save(
    type1_results,
    power_results,
    error_power_table,
    simulation_settings,
    file     = workspace_file,
    compress = "gzip"
  )
  
  cat("\nSaved summary table:", summary_file, "\n")
  cat("Saved workspace:", workspace_file, "\n")
  
  invisible(list(
    type1_results     = type1_results,
    power_results     = power_results,
    error_power_table = error_power_table,
    output_dir        = sign_output_dir,
    workspace_file    = workspace_file
  ))
}


## ============================================================
## 14. Main run function
## ============================================================

## Main user-facing function for running the full simulation study.
run_selective_inference_simulation <- function(
    num_simulations     = 1e6,
    distributions       = c("normal", "exponential", "lognormal"),
    sample_sizes        = c(10, 20, 30),
    effect_size         = 0.5,
    alpha_level         = 0.05,
    sign_test_types     = c("regular", "randomized"),
    output_dir          = "results",
    seed                = 12345,
    distribution_info   = make_default_distribution_info(),
    distribution_labels = make_default_distribution_labels(),
    save_pvalue_csv     = FALSE,
    max_qq_points       = 4000,
    n_cores             = max(1, parallel::detectCores() - 1)
) {
  
  sign_test_types <- match.arg(
    sign_test_types,
    choices = c("regular", "randomized"),
    several.ok = TRUE
  )
  
  ensure_output_dir(output_dir)
  
  method_settings <- make_method_settings()
  
  RNGversion("4.2.0")
  
  set.seed(
    seed = seed,
    kind        = "Mersenne-Twister",
    normal.kind = "Inversion",
    sample.kind = "Rejection"
  )
  
  summary_tables <- list()
  
  for (sign_test_type in sign_test_types) {
    
    one_result <- run_one_sign_test_workflow(
      sign_test_type       = sign_test_type,
      num_simulations      = num_simulations,
      distributions        = distributions,
      sample_sizes         = sample_sizes,
      effect_size          = effect_size,
      alpha_level          = alpha_level,
      method_settings      = method_settings,
      distribution_info    = distribution_info,
      distribution_labels  = distribution_labels,
      output_dir           = output_dir,
      save_pvalue_csv      = save_pvalue_csv,
      max_qq_points        = max_qq_points,
      n_cores              = n_cores
    )
    
    summary_tables[[sign_test_type]] <- one_result$error_power_table
    
    rm(one_result)
    gc()
  }
  
  combined_error_power_table <- do.call(rbind, summary_tables)
  
  combined_summary_file <- file.path(
    output_dir,
    "combined_type1_error_and_power_summary.csv"
  )
  
  write.csv(
    combined_error_power_table,
    file      = combined_summary_file,
    row.names = FALSE
  )
  
  
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat("Pipeline complete.\n")
  cat("Results saved in: ", output_dir, "\n", sep = "")
  cat("Combined summary file: ", combined_summary_file, "\n", sep = "")
  cat("Sign-test-specific workspaces saved in their respective folders.\n")
  cat(strrep("=", 70), "\n", sep = "")
  
  invisible(list(
    combined_error_power_table = combined_error_power_table,
    output_dir                 = output_dir,
    combined_summary_file      = combined_summary_file
  ))
}


## ============================================================
## 15. Run simulation
## ============================================================

final_results <- run_selective_inference_simulation(
  num_simulations = 1e6,
  distributions   = c("normal", "exponential", "lognormal"),
  sample_sizes    = c(10, 20, 30),
  effect_size     = 0.5,
  alpha_level     = 0.05,
  sign_test_types = c("regular", "randomized"),
  output_dir      = "selective_inf_results",
  seed            = 12345,
  save_pvalue_csv = FALSE,
  max_qq_points   = 4000,
  n_cores         = max(1, parallel::detectCores() - 1)
)