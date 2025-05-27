# Code is adapted from rlearner:
# R-learner for Quasi-Oracle Estimation of Heterogeneous Treatment Effects
rlasso <- function(W,
                   A,
                   Y,
                   g1W,
                   theta,
                   foldid,
                   weights = NULL,
                   max_degree = 1,
                   use_weight = TRUE, # much faster
                   browse = FALSE) {
  if (browse) browser()

  if (is.null(weights)) weights <- rep(1, length(Y))

  # make design matrix
  if (max_degree > 1) {
    phi_W <- model.matrix(as.formula(paste0("~ -1+(.)^", max_degree)), data = W)
  } else {
    phi_W <- as.matrix(W)
  }

  # option 1: pseudo outcome (Y-theta)/(A-g1W), pseudo weight (A-g)^2
  # option 2: pseudo outcome (Y-theta), pseudo design (A-g1W)phi_W
  if (use_weight) {
    pseudo_outcome <- ifelse(abs(A-g1W) < 1e-10, 0, (Y-theta)/(A-g1W))
    pseudo_weights <- (A-g1W)^2*weights
    pseudo_designs <- phi_W
  } else {
    pseudo_outcome <- Y-theta
    pseudo_weights <- weights
    pseudo_designs <- as.matrix((A-g1W)*cbind(1, phi_W))
  }

  # fit CATE and predict
  tau_fit <- cv.glmnet(x = pseudo_designs,
                       y = pseudo_outcome,
                       weights = pseudo_weights,
                       foldid = foldid,
                       alpha = 1,
                       standardize = FALSE)
  if (use_weight) {
    non_zero <- which(coef(tau_fit, s = "lambda.min")[-1] != 0)
    beta <- coef(tau_fit, s = "lambda.min")[c(1, non_zero+1)]
  } else {
    non_zero <- which(coef(tau_fit, s = "lambda.min")[-c(1, 2)] != 0)
    beta <- coef(tau_fit, s = "lambda.min")[c(2, non_zero+2)]
  }
  phi_W <- as.matrix(cbind(1, phi_W[, non_zero, drop=FALSE]))
  pred <- as.vector(phi_W %*% beta)

  ret_obj = list(tau_fit = tau_fit,
                 beta = beta,
                 phi_train = phi_W,
                 max_degree = max_degree,
                 non_zero = non_zero,
                 pred = pred,
                 pseudo_outcome = pseudo_outcome,
                 pseudo_weights = pseudo_weights)
  class(ret_obj) <- "rlasso"
  return(ret_obj)
}

predict.rlasso <- function(object,
                           newx = NULL) {
  if (!is.null(newx)) {
    if (object$max_degree > 1) {
      newx_aug <- model.matrix(as.formula(paste0("~ -1+(.)^", max_degree)), data = newx)
    } else {
      newx_aug <- as.matrix(newx)
    }
    return(as.vector(cbind(1, newx) %*% object$beta))
  } else {
    pred <- return(object$pred)
  }
}
