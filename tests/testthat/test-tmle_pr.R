library(atmle)

test_that("tmle_pr runs with pooled regression nuisances", {
  set.seed(20240519)
  n <- 180
  W1 <- rnorm(n)
  W2 <- rbinom(n, 1, 0.45)
  S <- rbinom(n, 1, plogis(-0.6 + 0.3 * W1 - 0.2 * W2))
  A <- rbinom(
    n,
    1,
    ifelse(S == 1, 0.5, plogis(-0.2 + 0.4 * W1 - 0.2 * W2))
  )
  Y <- 0.2 + 0.6 * A + 0.4 * W1 - 0.25 * W2 + rnorm(n, sd = 0.8)
  data <- data.frame(S, W1, W2, A, Y)

  method <- list(learners = list(sl3::Lrnr_glm$new()))

  fit <- tmle_pr$new(
    data = data,
    S_node = "S",
    W_nodes = c("W1", "W2"),
    A_node = "A",
    Y_node = "Y",
    family = "gaussian",
    n_folds = 3,
    seed = 20240519
  )

  fit$run(
    Q_method = method,
    g_method = method,
    Pi_bar_method = method
  )

  expect_s3_class(fit$results, "data.frame")
  expect_equal(nrow(fit$results), 2)
  expect_true(all(c("param", "psi", "se", "lower", "upper", "PnEIC") %in% names(fit$results)))
  expect_true(all(is.finite(fit$results$psi)))
  expect_true(all(is.finite(fit$results$se)))
  expect_true(length(fit$eic$pooled_W) == n)
  expect_true(length(fit$eic$rct_W) == n)
})
