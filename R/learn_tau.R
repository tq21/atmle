# function to learn tau, R-loss
learn_tau <- function(S, W, A, Y, Pi, theta,
                      controls_only,
                      method = "glmnet") {

  pred <- NULL
  A1 <- numeric(length = length(A))
  A0 <- numeric(length = length(A))
  X <- NULL
  X_A1_counter <- NULL
  X_A0_counter <- NULL
  x_basis <- NULL
  x_basis_A1 <- NULL
  x_basis_A0 <- NULL
  coefs <- NULL

  weights <- 1
  pseudo_outcome <- NULL
  pseudo_weights <- NULL

  if (controls_only) {
    pred <- numeric(length = length(A))

    # R-transformations, only controls
    pseudo_outcome <- (Y[A == 0]-theta$pred[A == 0])/(S[A == 0]-Pi$pred[A == 0])
    pseudo_weights <- (S[A == 0]-Pi$pred[A == 0])^2*weights

    # design matrix, only controls (X: W)
    X <- W[A == 0, ]

    # counterfactual design matrices
    X_A1_counter <- X_A0_counter <- cbind(1, W)

  } else {
    pred <- numeric(length = length(A))

    # R-transformations
    pseudo_outcome <- (Y-theta$pred)/(S-Pi$pred)
    pseudo_weights <- (S-Pi$pred)^2*weights

    # design matrix, both treated and controls (X: W*I(A == 1), I(A == 1), W*I(A == 0))
    X <- cbind(W * A, A, W * (1-A))

    # counterfactual design matrices
    zero_W <- matrix(0, nrow = nrow(W), ncol = ncol(W))
    X_A1_counter <- cbind(1, W, 1, zero_W)
    X_A0_counter <- cbind(1, zero_W, 0, W)
  }

  if (method == "glmnet") {
    fit <- cv.glmnet(x = as.matrix(X), y = pseudo_outcome, intercept = TRUE,
                     family = "gaussian", weights = pseudo_weights,
                     keep = TRUE, nfolds = 5, alpha = 1, relax = TRUE)
    non_zero <- which(as.numeric(coef(fit, s = "lambda.min")) != 0)
    coefs <- coef(fit, s = "lambda.min")[non_zero]

    if (controls_only) {
      # selected counterfactual bases
      x_basis <- x_basis_A0 <- as.matrix(cbind(1, W)[, non_zero, drop = FALSE])

      # predictions
      pred <- A0 <- as.numeric(x_basis_A0 %*% matrix(coefs))

    } else {
      # selected counterfactual bases
      x_basis <- as.matrix(cbind(1, X)[, non_zero, drop = FALSE])
      x_basis_A1 <- as.matrix(X_A1_counter[, non_zero, drop = FALSE])
      x_basis_A0 <- as.matrix(X_A0_counter[, non_zero, drop = FALSE])

      # predictions
      A1 <- as.numeric(x_basis_A1 %*% matrix(coefs))
      A0 <- as.numeric(x_basis_A0 %*% matrix(coefs))
      pred[A == 1] <- A1[A == 1]
      pred[A == 0] <- A0[A == 0]
    }
    print("bases selected: " %+% length(coefs))


  } else if (method == "HAL") {
    # TODO: not fully working right now
    X <- NULL
    if (controls_only) {
      X <- as.matrix(W[A == 0, ])
    } else {
      X <- as.matrix(cbind(W, A))
    }

    fit <- fit_relaxed_hal(X = X, Y = pseudo_outcome,
                           family = "gaussian",
                           weights = pseudo_weights)

    if (controls_only) {
      # design matrices
      x_basis <- x_basis_A0 <- make_counter_design_matrix(fit$basis_list, as.matrix(W))

      # predictions
      pred <- A0 <- as.numeric(x_basis_A0 %*% matrix(fit$beta))
    } else {
      # design matrices
      x_basis <- make_counter_design_matrix(fit$basis_list, X)
      x_basis_A1 <- make_counter_design_matrix(fit$basis_list, as.matrix(cbind(W, A = 1)))
      x_basis_A0 <- make_counter_design_matrix(fit$basis_list, as.matrix(cbind(W, A = 0)))

      # predictions
      A1 <- as.numeric(x_basis_A1 %*% matrix(fit$beta))
      A0 <- as.numeric(x_basis_A0 %*% matrix(fit$beta))
      pred[A == 1] <- A1[A == 1]
      pred[A == 0] <- A0[A == 0]
    }

    coefs <- fit$beta
    print("bases selected: " %+% length(coefs))
  }

  return(list(A1 = A1,
              A0 = A0,
              x_basis = x_basis,
              x_basis_A1 = x_basis_A1,
              x_basis_A0 = x_basis_A0,
              pred = pred,
              coefs = coefs))
}
