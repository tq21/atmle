test_that("floor weight and combined variance use paired influence curves", {
  df <- c(-3, -1, 1, 3)
  du <- c(-1, 1, -1, 1)
  u <- .12; f <- .1; n <- length(df)
  vf <- var(df)/n; vu <- var(du)/n; cf <- cov(du, df)/n
  vd <- var(du-df)/n
  w <- min(.5, max(0, (vf-cf)/(vd+max((u-f)^2, vd))))
  x <- .atmle_floor_pair(u, f, du, df)
  expect_equal(x$weight, w)
  expect_equal(x$psi, f+w*(u-f))
  expect_equal(x$se^2, (1-w)^2*vf+w^2*vu+2*w*(1-w)*cf)
  expect_equal(x$influence, (1-w)*df+w*du)
  expect_equal(x$difference_variance, vd)
  expect_equal(.atmle_floor_pair(10, 0, du, df)$bias_squared, 100)
})

test_that("degenerate pairs and convex boundaries are safe", {
  d <- c(-3, -1, 1, 3)
  expect_equal(.atmle_floor_pair(.2, .1, d, d)$weight, 0)
  expect_equal(.atmle_floor_pair(0, 0, 2*d, d)$weight, 0)
  expect_equal(.atmle_floor_pair(0, 0, .9*d, d)$weight, .5)
  expect_equal(.atmle_floor_pair(0, 0, .9*d, d, 1)$weight, 1)
  expect_equal(.atmle_floor_pair(0, 0, d, d, 0)$psi, 0)
  expect_equal(.atmle_floor_pair(0, 0, rep(0, 4), rep(0, 4))$se, 0)
  expect_error(.atmle_floor_pair(0, 0, d, d[-1]), "aligned")
  expect_error(.atmle_floor_pair(NA, 0, d, d), "finite")
  expect_equal(mat_inverse(matrix(0, 1, 1)), matrix(0, 1, 1))
})

test_that("reports add exactly one combination per population", {
  fit <- fixture_fit()
  expect_equal(nrow(fit$results), 4)
  expect_equal(as.integer(table(fit$results$param)), c(2L, 2L))
  expect_true(all(fit$results$converged))
  for (i in 3:4) {
    population <- fit$results$param[i]
    j <- which(fit$components$forced$results$param == population)
    f <- fit$components$forced$results$psi[j]
    u <- fit$components$unforced$results$psi[j]
    w <- fit$results$unforced_weight[i]
    expect_equal(fit$results$psi[i], f+w*(u-f))
    expect_equal(fit$results$se[i], sqrt(var(fit$influence[[i]])/nrow(fit$data)))
    expect_equal(fit$results$lower[i], fit$results$psi[i]-qnorm(.975)*fit$results$se[i])
  }
})

test_that("primary choice does not change the paired combination", {
  f <- fixture_fit()
  u <- make_fixture(forced = FALSE)
  expect_match(u$results$estimator[1], "unforced")
  expect_equal(f$results$psi[3:4], u$results$psi[3:4], tolerance = 1e-10)
  expect_equal(f$results$unforced_weight[3:4], u$results$unforced_weight[3:4], tolerance = 1e-10)
})

test_that("extra lambda rows are preserved without extra combinations", {
  fit <- make_fixture(n_lambda = 2L)
  base <- nrow(fit$components$forced$results)
  expect_equal(nrow(fit$results), base + 2L)
  expect_equal(sum(fit$results$estimator == "Variance-floor convex combination"), 2L)
  expect_equal(length(fit$influence), nrow(fit$results))
})

test_that("forced basis is retained even with a zero initial coefficient", {
  fit <- fixture_fit()
  tau <- fit$tau_S
  forced <- tau$forced_columns
  expect_length(forced, 1L)
  tau$fit$glmnet.fit$beta[forced-1L, ] <- 0
  models <- .atmle_working_models(tau, 1L, TRUE)
  expect_true(any(vapply(seq_len(ncol(models[[1]]$phi_WA)), function(j) {
    isTRUE(all.equal(as.numeric(models[[1]]$phi_WA[, j]), fit$A))
  }, logical(1))))
})

test_that("inference updates Wald limits while preserving estimates and variance", {
  fit <- fixture_fit()
  estimates <- fit$results$psi
  standard_errors <- fit$results$se
  influence <- fit$influence
  fit$inference(alpha = .1)
  expect_identical(fit$results$psi, estimates)
  expect_identical(fit$results$se, standard_errors)
  expect_identical(fit$influence, influence)
  expect_equal(fit$results$lower, estimates-qnorm(.95)*standard_errors)
  expect_equal(fit$results$upper, estimates+qnorm(.95)*standard_errors)
})

test_that("binary outcomes, external controls, and relaxed targeting work", {
  for (config in list(list(family = "binomial", target_gwt = FALSE),
                       list(controls_only = TRUE, target_gwt = FALSE),
                       list(controls_only = TRUE),
                       list(family = "binomial", controls_only = TRUE),
                       list(target_method = "relaxed"))) {
    fit <- do.call(make_fixture, config)
    expect_true(all(is.finite(fit$results$psi)))
    expect_true(all(fit$results$converged))
    expect_true(all(is.finite(fit$results$lower)))
  }
})

test_that("invalid inference controls fail before estimation", {
  expect_error(.atmle_validate_inference(1, .5), "alpha")
  expect_error(.atmle_validate_inference(.05, 2), "weight_cap")
})
