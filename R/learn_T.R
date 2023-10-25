learn_T <- function(W, A, Y,
                    g, theta_tilde,
                    method,
                    v_folds) {
  # R-transformations
  weights <- 1 # TODO: currently not used
  pseudo_outcome <- (Y-theta_tilde)/(A-g)
  pseudo_weights <- (A-g)^2*weights

  # initialize
  pred <- NULL
  x_basis <- NULL
  coefs <- NULL

  if (method == "glmnet") {
    fit <- cv.glmnet(x = as.matrix(W), y = pseudo_outcome,
                     family = "gaussian", weights = pseudo_weights,
                     keep = TRUE, nfolds = v_folds, alpha = 1, relax = TRUE)
    coefs <- coef(fit, s = "lambda.min", gamma = 0)
    pred <- as.numeric(as.matrix(cbind(1, W)) %*% matrix(coefs))

    # design matrix
    x_basis <- as.matrix(cbind(1, W))
  } else if (method == "HAL") {
    # TODO: currently not working
    fit <- fit_relaxed_hal(as.matrix(data.table(W)), pseudo_outcome,
                           "gaussian", weights = pseudo_weights)
    pred <- as.numeric(fit$pred)
    x_basis <- make_counter_design_matrix(fit$basis_list, as.matrix(data.table(W)))
  }

  return(list(pred = pred,
              x_basis = x_basis,
              coefs = coefs))
}
