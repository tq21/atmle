make_pr_data <- function(family) {
  set.seed(215)
  n <- 400
  data <- data.frame(W1 = rnorm(n), W2 = rbinom(n, 1, .5),
                     S = rbinom(n, 1, .4), A = rbinom(n, 1, .5))
  mu <- with(data, -.6 + .3*W1 + .2*W2 + .5*A + .3*S)
  data$Y <- if (family == "binomial") rbinom(n, 1, plogis(mu)) else mu + rnorm(n)
  data
}

test_that("tmle_pr is exported", {
  expect_true("tmle_pr" %in% getNamespaceExports("atmle"))
})

test_that("targeting and inference reproduce frozen original implementation outputs", {
  # Generated directly from the original atmleFusion source, not the port.
  fixtures <- readRDS(test_path("fixtures", "tmle-pr-original.rds"))
  for (x in fixtures) {
    fit <- tmle_pr$new(x$data, "S", "W", "A", "Y", x$family, n_folds = 3)
    fit$Q <- x$initial$Q
    fit$g11W <- x$initial$g
    fit$g01W <- 1 - fit$g11W
    fit$Pi_bar <- x$initial$Pi
    fit$Q_star <- fit$target_Q(FALSE, g_bound = x$g_bound)
    fit$Q_star_S1 <- fit$target_Q(TRUE, g_bound = x$g_bound)
    fit$inference(g_bound = x$g_bound)
    expect_equal(fit$Q_star, x$Q_star, tolerance = 1e-12)
    expect_equal(fit$Q_star_S1, x$Q_star_S1, tolerance = 1e-12)
    expect_equal(fit$results, x$results, tolerance = 1e-12)
    expect_equal(fit$eic, x$eic, tolerance = 1e-12)
  }
})

test_that("both outcome families run with all dependencies resolved by the package", {
  for (family in c("binomial", "gaussian")) {
    data <- make_pr_data(family)
    method <- list(learners = list(sl3::Lrnr_glm$new()))
    fit <- tmle_pr$new(data, "S", c("W1", "W2"), "A", "Y", family,
                       n_folds = 3, seed = 42)
    fit$run(method, method, method)
    expect_equal(fit$results$param,
                 c("ATE (avg. over pooled)", "ATE (avg. over RCT)"))
    expect_true(all(is.finite(as.matrix(fit$results[c("psi", "se", "lower", "upper")]))))
    expect_true(all(abs(fit$results$PnEIC) < 1e-6))
    expect_equal(fit$Q$Q1WA, ifelse(data$A == 1, fit$Q$Q1W1, fit$Q$Q1W0))
    expect_equal(fit$results$psi[2],
                 mean((fit$Q_star_S1$Q1W1 - fit$Q_star_S1$Q1W0)[data$S == 1]))
    expect_equal(fit$eic$rct_W[data$S == 0], rep(0, sum(data$S == 0)))
    expect_equal(fit$results$eic[[1]], fit$eic$pooled_W)
    expect_equal(fit$results$eic[[2]], fit$eic$rct_W)
    estimate <- fit$results$psi
    standard_error <- fit$results$se
    fit$inference(alpha = .1)
    expect_equal(fit$results$psi, estimate)
    expect_equal(fit$results$se, standard_error)
    expect_equal(fit$results$lower, estimate - qnorm(.95)*standard_error)
    expect_equal(fit$results$upper, estimate + qnorm(.95)*standard_error)
  }
})

test_that("discrete Super Learner works through the existing regression helper", {
  data <- make_pr_data("binomial")
  method <- function() list(
    learners = list(sl3::Lrnr_glm$new(), sl3::Lrnr_mean$new()),
    metalearner = sl3::Lrnr_cv_selector$new(sl3::loss_loglik_binomial)
  )
  fit <- tmle_pr$new(data, "S", c("W1", "W2"), "A", "Y", "binomial", 3)
  fit$run(method(), method(), method())
  expect_true(all(is.finite(fit$results$psi)))
  expect_true(all(is.finite(fit$results$se)))
  expect_true(all(abs(fit$results$PnEIC) < 1e-6))
})

test_that("known treatment probabilities skip estimation and family objects work", {
  data <- make_pr_data("binomial")
  fit <- tmle_pr$new(data, "S", c("W1", "W2"), "A", "Y", binomial(), 3)
  fit$g11W <- rep(.5, nrow(data))
  method <- list(learners = list(sl3::Lrnr_glm$new()))
  fit$run(method, g_method = NULL, Pi_bar_method = method)
  expect_null(fit$g_fit)
  expect_equal(fit$g11W, rep(.5, nrow(data)))
  expect_equal(fit$g01W, rep(.5, nrow(data)))
  expect_true(all(is.finite(fit$results$psi)))
  expect_error(tmle_pr$new(data, "S", "W1", "A", "Y", "poisson"), "family")
})
