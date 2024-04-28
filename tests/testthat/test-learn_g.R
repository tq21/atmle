test_that("learn_g works when external has both treated and controls", {
  # simulate data
  set.seed(123)
  n <- 500
  S <- rbinom(n, 1, 0.5)
  W1 <- rnorm(n)
  W2 <- rnorm(n)
  W <- cbind(W1, W2)
  A <- numeric(n)
  A[S == 1] <- rbinom(sum(S), 1, 0.67)
  A[S == 0] <- rbinom(n - sum(S), 1, plogis(1.2 * W1[S == 0] - 0.9 * W2[S == 0]))

  # make cv folds
  cv_strata <- paste0(S, "-", A)
  suppressWarnings({
    folds <- make_folds(
      n = n, V = 5,
      strata_ids = as.integer(factor(cv_strata))
    )
  })

  # glm
  g <- learn_g(
    S = S,
    W = W,
    A = A,
    g_rct = 0.67,
    controls_only = FALSE,
    method = "glm",
    folds = folds,
    g_bounds = c(0.01, 0.5)
  )
  expect_true(min(g) >= 0.01 & max(g) <= 0.5)
  expect_equal(length(g), n)

  # glmnet
  g <- learn_g(
    S = S,
    W = W,
    A = A,
    g_rct = 0.67,
    controls_only = FALSE,
    method = "glmnet",
    folds = folds,
    g_bounds = c(0.01, 0.99)
  ) # glmnet
  expect_true(min(g) >= 0.01 & max(g) <= 0.99)
  expect_equal(length(g), n)

  # # sl3 with default learners
  # g <- learn_g(S = S,
  #              W = W,
  #              A = A,
  #              g_rct = 0.67,
  #              controls_only = FALSE,
  #              method = "sl3",
  #              folds = folds,
  #              g_bounds = c(0.01, 0.99))
  # expect_true(min(g) >= 0.01 & max(g) <= 0.99)
  # expect_equal(length(g), n)
  #
  # # sl3 with custom learners
  # lrnrs <- list(Lrnr_mean$new(), Lrnr_glm$new())
  # g <- learn_g(S = S,
  #              W = W,
  #              A = A,
  #              g_rct = 0.67,
  #              controls_only = FALSE,
  #              method = lrnrs,
  #              folds = folds,
  #              g_bounds = c(0.01, 0.99))
  # expect_true(min(g) >= 0.01 & max(g) <= 0.99)
  # expect_equal(length(g), n)

  expect_error(learn_g(
    S = S,
    W = W,
    A = A,
    g_rct = 0.67,
    controls_only = FALSE,
    method = "abc",
    folds = folds,
    g_bounds = c(0.01, 0.99)
  ))
})

test_that("learn_g works when external has controls only", {
  # simulate data
  set.seed(123)
  n <- 500
  S <- rbinom(n, 1, 0.5)
  W1 <- rnorm(n)
  W2 <- rnorm(n)
  W <- cbind(W1, W2)
  A <- numeric(n)
  A[S == 1] <- rbinom(sum(S), 1, 0.67)
  A[S == 0] <- 0

  # make cv folds
  cv_strata <- paste0(S, "-", A)
  suppressWarnings({
    folds <- make_folds(
      n = n, V = 5,
      strata_ids = as.integer(factor(cv_strata))
    )
  })

  # glm
  g <- learn_g(
    S = S,
    W = W,
    A = A,
    g_rct = 0.67,
    controls_only = TRUE,
    method = "glm",
    folds = folds,
    g_bounds = c(0.01, 0.5)
  )
  expect_true(min(g) >= 0.01 & max(g) <= 0.5)
  expect_equal(length(g), n)

  # glmnet
  g <- learn_g(
    S = S,
    W = W,
    A = A,
    g_rct = 0.67,
    controls_only = TRUE,
    method = "glmnet",
    folds = folds,
    g_bounds = c(0.01, 0.99)
  ) # glmnet
  expect_true(min(g) >= 0.01 & max(g) <= 0.99)
  expect_equal(length(g), n)

  # # sl3 with default learners
  # g <- learn_g(S = S,
  #              W = W,
  #              A = A,
  #              g_rct = 0.67,
  #              controls_only = TRUE,
  #              method = "sl3",
  #              folds = folds,
  #              g_bounds = c(0.01, 0.99))
  # expect_true(min(g) >= 0.01 & max(g) <= 0.99)
  # expect_equal(length(g), n)
  #
  # # sl3 with custom learners
  # lrnrs <- list(Lrnr_mean$new(), Lrnr_glm$new())
  # g <- learn_g(S = S,
  #              W = W,
  #              A = A,
  #              g_rct = 0.67,
  #              controls_only = TRUE,
  #              method = lrnrs,
  #              folds = folds,
  #              g_bounds = c(0.01, 0.99))
  # expect_true(min(g) >= 0.01 & max(g) <= 0.99)
  # expect_equal(length(g), n)

  expect_error(learn_g(
    S = S,
    W = W,
    A = A,
    g_rct = 0.67,
    controls_only = FALSE,
    method = "abc",
    folds = folds,
    g_bounds = c(0.01, 0.99)
  ))
})

# TODO: TEST CROSS FITTING FALSE
