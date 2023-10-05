# function to learn tau, R-loss
learn_tau <- function(S, W, A, Y, Pi, theta, controls_only, method = "glmnet") {

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
    fit <- cv.glmnet(x = as.matrix(X) , y = pseudo_outcome, intercept = TRUE,
                     family = "gaussian", weights = pseudo_weights,
                     keep = TRUE, nfolds = 5, alpha = 1, relax = TRUE)
    coefs <- coef(fit, s = "lambda.min")

    if (controls_only) {
      A0 <- as.numeric(as.matrix(cbind(1, W)) %*% matrix(coefs))
      pred <- A0

      # design matrices, only controls
      non_zero <- which(as.numeric(coef(fit, s = "lambda.min")) != 0)
      x_basis <- as.matrix(cbind(1, W)[, non_zero, drop = FALSE])
      x_basis_A0 <- x_basis

    } else {
      A1 <- as.numeric(as.matrix(X_A1_counter) %*% matrix(coefs))
      A0 <- as.numeric(as.matrix(X_A0_counter) %*% matrix(coefs))
      pred[A == 1] <- A1[A == 1]
      pred[A == 0] <- A0[A == 0]

      # design matrices
      non_zero <- which(as.numeric(coef(fit, s = "lambda.min")) != 0)
      x_basis <- as.matrix(cbind(1, X)[, non_zero, drop = FALSE])
      x_basis_A1 <- as.matrix(X_A1_counter[, non_zero, drop = FALSE])
      x_basis_A0 <- as.matrix(X_A0_counter[, non_zero, drop = FALSE])
    }

  } else if (method == "HAL") {
    # TODO: not working right now
    fit <- fit_relaxed_hal(X = X, Y = pseudo_outcome,
                           family = "gaussian",
                           weights = pseudo_weights)
    A1 <- as.numeric(make_counter_design_matrix(fit_A1$basis_list, as.matrix(W), add_main_terms = add_main_terms) %*% matrix(fit_A1$beta))
  }

  return(list(A1 = A1,
              A0 = A0,
              x_basis = x_basis,
              x_basis_A1 = x_basis_A1,
              x_basis_A0 = x_basis_A0,
              pred = pred,
              coefs = coefs))
}
