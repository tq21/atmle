library(atmle)
library(sl3)
source("tests/testthat/utils.R")

test_learn_Pi <- function(test_args) {
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

  res <- learn_Pi(S = S,
                  W = W,
                  A = A,
                  controls_only = test_args$controls_only,
                  method = test_args$method,
                  folds = folds,
                  Pi_bounds = test_args$Pi_bounds,
                  cross_fit_nuisance = test_args$cross_fit_nuisance)

  return(res)
}

test_that("learn_Pi basic functionalities", {
  # tabulate test scenarios
  lrnrs <- list(Lrnr_mean$new(), Lrnr_glm$new())
  test_args <- list.expand(n = list(500),
                           controls_only = list(TRUE, FALSE),
                           family = "gaussian",
                           prop_miss = list(0),
                           method = list("glm", "glmnet", lrnrs),
                           Pi_bounds = list(c(0.1, 0.2), c(0.01, 0.99)),
                           cross_fit_nuisance = list(TRUE, FALSE))

  for (i in 1:length(test_args)) {
    cur_args <- test_args[[i]]
    res <- test_learn_Pi(cur_args)

    # checks
    expect_length(res$pred, cur_args$n)
    expect_length(res$A0, cur_args$n)
    expect_length(res$A1, cur_args$n)
    expect_true(all(res$A0 >= cur_args$Pi_bounds[1] & res$A0 <= cur_args$Pi_bounds[2]))
    if (cur_args$controls_only) {
      expect_true(all(res$A1 == 1))
    } else {
      expect_true(all(res$A1 >= cur_args$Pi_bounds[1] & res$A1 <= cur_args$Pi_bounds[2]))
    }
  }
})

test_that("learn_Pi error behavior", {
  test_args <- list(n = 500,
                    controls_only = FALSE,
                    family = "gaussian",
                    prop_miss = 0,
                    method = "abc",
                    Pi_bounds = c(0.01, 0.99),
                    cross_fit_nuisance = TRUE)
  expect_error(test_learn_Pi(test_args))
})
