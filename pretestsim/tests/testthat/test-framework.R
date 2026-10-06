test_that("bundled framework loads without remote assets or user objects", {
  fw <- load_framework()
  expect_true(is.function(fw$run_simulation))
  expect_true(is.function(fw$fisher_combined))
  expect_false(exists("load_ml_resources", fw, inherits = FALSE))
  expect_false(exists("make_ml_normality_configs", fw, inherits = FALSE))
  expect_false("asset_url" %in% getNamespaceExports("pretestsim"))
  set.seed(30)
  x <- rnorm(20)
  expect_equal(unname(fw$generate_tests(x, "SW")), stats::shapiro.test(x)$p.value)
  expect_equal(fw$one_sample_t_test(x)$p.value, stats::t.test(x)$p.value)
})

test_that("Fisher custom test uses the retained framework interface", {
  skip_if_not_installed("nortest")
  fw <- load_framework()
  set.seed(30); x <- rexp(30)
  value <- fw$fisher_combined(x)
  expected <- pchisq(-2 * sum(log(c(shapiro.test(x)$p.value, nortest::ad.test(x)$p.value))), 4, lower.tail = FALSE)
  expect_equal(unname(value), expected)
  expect_named(value, "p.value")
})

test_that("contaminated normal IQR inverts the entire mixture CDF", {
  fw <- load_framework()
  value <- fw$contaminated_normal_iqr(.75, 0, 1, 5)
  expect_true(is.finite(value) && value > 0)
  q <- value / 2
  expect_equal(.75 * pnorm(q) + .25 * pnorm(q, sd = 5), .75, tolerance = 1e-9)
  set.seed(12)
  expect_length(fw$generate_data(20, "contaminated", c(.75, 0, 1, 5)), 20)
})
