library(atmle)
library(sl3)
source("tests/testthat/utils.R")

test_learn_g <- function(test_args) {
  set.seed(123)
  data <- sim_data(test_args$n,
                   test_args$controls_only,
                   test_args$family,
                   test_args$prop_miss)
  S <- data$S; W <- data[, c("W1", "W2")]; A <- data$A; Y <- data$Y
  delta <- as.integer(!is.na(Y))
  cv_strata <- paste0(S, "-", A)
  suppressWarnings({
    folds <- make_folds(
      n = test_args$n, V = 5,
      strata_ids = as.integer(factor(cv_strata))
    )
  })

  res <- learn_g(W = W,
                 A = A,
                 method = test_args$method,
                 folds = folds,
                 g_bounds = test_args$g_bounds)

  return(res)
}

test_that("learn_g basic functionalities", {
  # tabulate test scenarios
  lrnrs <- list(Lrnr_mean$new(), Lrnr_glm$new())
  test_args <- list.expand(n = list(500),
                           controls_only = list(FALSE),
                           family = list("gaussian"),
                           prop_miss = list(0),
                           method = list("glm", "glmnet", lrnrs),
                           g_bounds = list(c(0.1, 0.5), c(0.01, 0.99)),
                           cross_fit_nuisance = list(TRUE, FALSE))

  for (i in 1:length(test_args)) {
    cur_args <- test_args[[i]]
    res <- test_learn_g(cur_args)

    # checks
    expect_true(all(res >= cur_args$g_bounds[1] & res <= cur_args$g_bounds[2]))
    expect_equal(length(res), cur_args$n)
  }
})

test_that("learn_g error behavior", {
  test_args <- list(n = 500,
                    controls_only = FALSE,
                    family = "gaussian",
                    prop_miss = 0,
                    method = "abc",
                    g_bounds = c(0.01, 0.99),
                    cross_fit_nuisance = TRUE)
  expect_error(test_learn_g(test_args))
})
