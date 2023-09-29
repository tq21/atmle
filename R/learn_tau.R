# function to learn tau, R-loss
learn_tau <- function(S, W, A, Y, Pi, theta, method = "lasso") {
  weights <- 1
  pseudo_outcome <- (Y-theta$pred)/(S-Pi$pred)
  pseudo_weights <- (S-Pi$pred)^2*weights

  pred <- numeric(length = length(A))
  A1 <- numeric(length = length(A))
  A0 <- numeric(length = length(A))
  x_basis <- NULL
  x_basis_A1 <- NULL
  x_basis_A0 <- NULL

  # design matrix: (intercept, W, A, W * A)
  X_A1A0 <- cbind(W * A, A, W * (1-A))

  # counterfactual design matrices
  zero_W <- matrix(0, nrow = nrow(W), ncol = ncol(W))
  X_A1_counter <- cbind(1, W, 1, zero_W)
  X_A0_counter <- cbind(1, zero_W, 0, W)

  if (method == "lasso") {
    fit <- cv.glmnet(x = as.matrix(X_A1A0) , y = pseudo_outcome, intercept = TRUE,
                     family = "gaussian", weights = pseudo_weights,
                     keep = TRUE, nfolds = 5, alpha = 1, relax = TRUE)
    A1 <- as.numeric(as.matrix(X_A1_counter) %*% matrix(coef(fit, s = "lambda.min")))
    A0 <- as.numeric(as.matrix(X_A0_counter) %*% matrix(coef(fit, s = "lambda.min")))
    pred[A == 1] <- A1[A == 1]
    pred[A == 0] <- A0[A == 0]

    # design matrices
    non_zero <- which(as.numeric(coef(fit, s = "lambda.min")) != 0)
    x_basis <- as.matrix(cbind(1, X_A1A0)[, non_zero, drop = FALSE])
    x_basis_A1 <- as.matrix(X_A1_counter[, non_zero, drop = FALSE])
    x_basis_A0 <- as.matrix(X_A0_counter[, non_zero, drop = FALSE])
  } else if (method == "HAL") {
    fit <- fit_relaxed_hal(X = X_A1A0, Y = pseudo_outcome,
                           family = "gaussian",
                           weights = pseudo_weights)
    A1 <- as.numeric(make_counter_design_matrix(fit_A1$basis_list, as.matrix(W), add_main_terms = add_main_terms) %*% matrix(fit_A1$beta))
  } else if (method == "glm") {
    fit <- glm.fit(x = X_A1A0, y = pseudo_outcome,
                   weights = pseudo_weights, intercept = FALSE)
    A1 <- as.numeric(as.matrix(X_A1_counter) %*% matrix(coef(fit)))
    A0 <- as.numeric(as.matrix(X_A0_counter) %*% matrix(coef(fit)))

    x_basis <- as.matrix(X_A1A0)
    x_basis_A1 <- as.matrix(X_A1_counter)
    x_basis_A0 <- as.matrix(X_A0_counter)
  }

  return(list(A1 = A1,
              A0 = A0,
              x_basis = x_basis,
              x_basis_A1 = x_basis_A1,
              x_basis_A0 = x_basis_A0,
              pred = pred))
}
