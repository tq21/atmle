# function to learn psi_tilde, R-learner
learn_psi_tilde <- function(W, A, Y, g, theta_tilde,
                            method = "glmnet") {
  weights <- 1
  pseudo_outcome <- (Y-theta_tilde)/(A-g)
  pseudo_weights <- (A-g)^2*weights

  pred <- NULL
  x_basis <- NULL
  coefs <- NULL
  if (method == "glmnet") {
    fit <- cv.glmnet(x = as.matrix(W),
                     y = pseudo_outcome, family = "gaussian", weights = pseudo_weights,
                     keep = TRUE, nfolds = 5, alpha = 1, relax = TRUE)
    coefs <- coef(fit, s = "lambda.min", gamma = 0)
    pred <- as.numeric(as.matrix(cbind(1, W)) %*% matrix(coefs))

    # design matrix
    x_basis <- as.matrix(cbind(1, W))
  } else if (method == "HAL") {
    fit <- fit_relaxed_hal(as.matrix(data.table(W)), pseudo_outcome,
                           "gaussian", weights = pseudo_weights)
    pred <- as.numeric(fit$pred)
    x_basis <- make_counter_design_matrix(fit$basis_list, as.matrix(data.table(W)))
  } else if (method == "glm") {
    fit <- lm(pseudo_outcome ~ ., data = data.frame(W))
    pred <- as.numeric(predict(fit))
    x_basis <- cbind(1, as.matrix(data.table(W)))
  }

  return(list(pred = pred,
              x_basis = x_basis,
              coefs = coefs))
}
