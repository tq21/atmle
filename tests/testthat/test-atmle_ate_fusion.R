library(atmle)

test_that("atmle_ate_fusion runs with HAL working models", {
  set.seed(147)
  n <- 160
  W1 <- rnorm(n)
  W2 <- rnorm(n)
  S <- rbinom(n, 1, plogis(-0.3 + 0.25 * W1))
  A <- rbinom(n, 1, plogis(0.1 * S + 0.2 * W1 - 0.1 * W2))
  Y <- 0.5 + 0.4 * W1 + 0.2 * W2 + A + 0.2 * (1 - S) * W2 + rnorm(n)
  data <- data.frame(S, W1, W2, A, Y)

  method <- list(learners = list(sl3::Lrnr_glm$new()))

  fit <- atmle_ate_fusion$new(
    data = data,
    S_node = "S",
    W_nodes = c("W1", "W2"),
    A_node = "A",
    Y_node = "Y",
    family = "gaussian",
    n_folds = 3
  )

  fit$run(
    g_bar_method = method,
    theta_method = method,
    Pi_method = method,
    Q_bar_method = method,
    A_cate_args = list(max_degree = 1L, smoothness_orders = 1L, num_knots = 5L),
    S_cate_args = list(max_degree = 1L, smoothness_orders = 1L, num_knots = 5L),
    target_method = "tmle",
    target_gwt = FALSE,
    max_iter = 1,
    n_lambda = 1,
    parallel = FALSE,
    verbose = FALSE
  )

  expect_s3_class(fit$results, "data.frame")
  expect_equal(nrow(fit$results), 2)
  expect_true(all(c("param", "psi", "se", "lower", "upper") %in% names(fit$results)))
  expect_true(all(is.finite(fit$results$psi)))
  expect_false("diagnose" %in% names(fit$.__enclos_env__$public_methods))
  expect_false("cate_fit_method" %in% names(formals(fit$run)))
})
