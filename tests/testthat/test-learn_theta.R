library(atmle)
library(sl3)
source("tests/testthat/utils.R")

test_learn_theta <- function(test_args) {
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

  res <- learn_theta(W = W,
                     A = A,
                     Y = Y,
                     delta = delta,
                     controls_only = test_args$controls_only,
                     method = test_args$method,
                     folds = folds,
                     family = test_args$family,
                     theta_bounds = test_args$theta_bounds,
                     cross_fit_nuisance = test_args$cross_fit_nuisance)

  return(res)
}

test_that("learn_theta basic functionalities", {
  # tabulate test scenarios
  lrnrs <- list(Lrnr_mean$new(), Lrnr_glm$new())
  test_args <- list.expand(n = list(500),
                           controls_only = list(TRUE, FALSE),
                           family = list("gaussian", "binomial"),
                           prop_miss = list(0, 0.1),
                           method = list("glm", "glmnet", lrnrs),
                           theta_bounds = list(c(0, 0.5), c(-4, 4)),
                           cross_fit_nuisance = list(TRUE, FALSE))

  for (i in 1:length(test_args)) {
    cur_args <- test_args[[i]]
    res <- test_learn_theta(cur_args)
    res[is.na(res)] <- 0 # TODO: change later

    # checks
    expect_true(all(res >= cur_args$theta_bounds[1] & res <= cur_args$theta_bounds[2]))
    expect_equal(length(res), cur_args$n)
    if (cur_args$family == "binomial") {
      expect_true(all(res >= 0 & res <= 1))
    }
  }
})

test_that("learn_theta error behavior", {
  test_args <- list(n = 500,
                    controls_only = FALSE,
                    family = "gaussian",
                    prop_miss = 0,
                    method = "abc",
                    theta_bounds = c(0, 0.5),
                    cross_fit_nuisance = TRUE)
  expect_error(test_learn_theta(test_args))
})
