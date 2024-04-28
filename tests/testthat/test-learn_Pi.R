test_that("learn_Pi works when external has both treated and controls", {
  # simulate data
  set.seed(123)
  n <- 500
  S <- rbinom(n, 1, 0.5)
  W1 <- rnorm(n); W2 <- rnorm(n); W <- cbind(W1, W2)
  A <- numeric(n)
  A[S == 1] <- rbinom(sum(S), 1, 0.67)
  A[S == 0] <- rbinom(n-sum(S), 1, plogis(1.2*W1[S == 0]-0.9*W2[S == 0]))

  # make cv folds
  cv_strata <- paste0(S, "-", A)
  suppressWarnings({
    folds <- make_folds(n = n, V = 5,
                        strata_ids = as.integer(factor(cv_strata)))
  })

  # glm
  Pi <- learn_Pi(S = S,
                 W = W,
                 A = A,
                 controls_only = FALSE,
                 method = "glm",
                 folds = folds,
                 Pi_bounds = c(0.01, 0.2))
  expect_true(all(Pi$pred >= 0.01 & Pi$pred <= 0.2))
  expect_true(all(Pi$A0 >= 0.01 & Pi$A0 <= 0.2))
  expect_true(all(Pi$A1 >= 0.01 & Pi$A1 <= 0.2))
  expect_length(Pi$pred, n)
  expect_length(Pi$A0, n)
  expect_length(Pi$A1, n)

  # glmnet
  Pi <- learn_Pi(S = S,
                 W = W,
                 A = A,
                 controls_only = FALSE,
                 method = "glmnet",
                 folds = folds,
                 Pi_bounds = c(0.01, 0.99))
  expect_true(all(Pi$pred >= 0.01 & Pi$pred <= 0.99))
  expect_true(all(Pi$A0 >= 0.01 & Pi$A0 <= 0.99))
  expect_true(all(Pi$A1 >= 0.01 & Pi$A1 <= 0.99))
  expect_length(Pi$pred, n)
  expect_length(Pi$A0, n)
  expect_length(Pi$A1, n)

  # # sl3 with default learners
  # Pi <- learn_Pi(S = S,
  #                W = W,
  #                A = A,
  #                controls_only = FALSE,
  #                method = "sl3",
  #                folds = folds,
  #                Pi_bounds = c(0.01, 0.99))
  # expect_true(all(Pi$pred >= 0.01 & Pi$pred <= 0.99))
  # expect_true(all(Pi$A0 >= 0.01 & Pi$A0 <= 0.99))
  # expect_true(all(Pi$A1 >= 0.01 & Pi$A1 <= 0.99))
  # expect_length(Pi$pred, n)
  # expect_length(Pi$A0, n)
  # expect_length(Pi$A1, n)
#
  # # sl3 with custom learners
  # lrnrs <- list(Lrnr_mean$new(), Lrnr_glm$new())
  # Pi <- learn_Pi(S = S,
  #                W = W,
  #                A = A,
  #                controls_only = FALSE,
  #                method = lrnrs,
  #                folds = folds,
  #                Pi_bounds = c(0.01, 0.99))
  # expect_true(all(Pi$pred >= 0.01 & Pi$pred <= 0.99))
  # expect_true(all(Pi$A0 >= 0.01 & Pi$A0 <= 0.99))
  # expect_true(all(Pi$A1 >= 0.01 & Pi$A1 <= 0.99))
  # expect_length(Pi$pred, n)
  # expect_length(Pi$A0, n)
  # expect_length(Pi$A1, n)

  expect_error(learn_Pi(S = S,
                        W = W,
                        A = A,
                        controls_only = FALSE,
                        method = "abc",
                        folds = folds,
                        Pi_bounds = c(0.01, 0.99)))
})

test_that("learn_Pi works when external has controls only", {
  # simulate data
  set.seed(123)
  n <- 500
  S <- rbinom(n, 1, 0.5)
  W1 <- rnorm(n); W2 <- rnorm(n); W <- cbind(W1, W2)
  A <- numeric(n)
  A[S == 1] <- rbinom(sum(S), 1, 0.67)
  A[S == 0] <- 0

  # make cv folds
  cv_strata <- paste0(S, "-", A)
  suppressWarnings({
    folds <- make_folds(n = n, V = 5,
                        strata_ids = as.integer(factor(cv_strata)))
  })

  # glm
  Pi <- learn_Pi(S = S,
                 W = W,
                 A = A,
                 controls_only = TRUE,
                 method = "glm",
                 folds = folds,
                 Pi_bounds = c(0.01, 0.2))
  expect_true(all(Pi$A0[A == 0] >= 0.01 & Pi$A0[A == 0] <= 0.2))
  expect_true(all(Pi$pred[A == 1] == 1))
  expect_length(Pi$pred, n)
  expect_length(Pi$A0, n)
  expect_length(Pi$A1, n)

  # glmnet
  Pi <- learn_Pi(S = S,
                 W = W,
                 A = A,
                 controls_only = TRUE,
                 method = "glmnet",
                 folds = folds,
                 Pi_bounds = c(0.01, 0.99))
  expect_true(all(Pi$A0[A == 0] >= 0.01 & Pi$A0[A == 0] <= 0.99))
  expect_true(all(Pi$pred[A == 1] == 1))
  expect_length(Pi$pred, n)
  expect_length(Pi$A0, n)
  expect_length(Pi$A1, n)

  # # sl3 with default learners
  # Pi <- learn_Pi(S = S,
  #                W = W,
  #                A = A,
  #                controls_only = TRUE,
  #                method = "sl3",
  #                folds = folds,
  #                Pi_bounds = c(0.01, 0.99))
  # expect_true(all(Pi$A0[A == 0] >= 0.01 & Pi$A0[A == 0] <= 0.99))
  # expect_true(all(Pi$pred[A == 1] == 1))
  # expect_length(Pi$pred, n)
  # expect_length(Pi$A0, n)
  # expect_length(Pi$A1, n)
#
  # # sl3 with custom learners
  # lrnrs <- list(Lrnr_mean$new(), Lrnr_glm$new())
  # Pi <- learn_Pi(S = S,
  #                W = W,
  #                A = A,
  #                controls_only = TRUE,
  #                method = lrnrs,
  #                folds = folds,
  #                Pi_bounds = c(0.01, 0.99))
  # expect_true(all(Pi$A0[A == 0] >= 0.01 & Pi$A0[A == 0] <= 0.99))
  # expect_true(all(Pi$pred[A == 1] == 1))
  # expect_length(Pi$pred, n)
  # expect_length(Pi$A0, n)
  # expect_length(Pi$A1, n)

  expect_error(learn_Pi(S = S,
                        W = W,
                        A = A,
                        controls_only = TRUE,
                        method = "abc",
                        folds = folds,
                        Pi_bounds = c(0.01, 0.99)))
})
