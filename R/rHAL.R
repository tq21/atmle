# Code is adapted from rlearner:
# R-learner for Quasi-Oracle Estimation of Heterogeneous Treatment Effects
rHAL <- function(W,
                 A,
                 Y,
                 g1W,
                 theta,
                 foldid,
                 weights = NULL,
                 enumerate_basis_args = list(max_degree = 2,
                                             smoothness_orders = 1,
                                             num_knots = c(20, 5)),
                 use_weight = TRUE, # much faster
                 browse = FALSE) {
  if (browse) browser()

  if (is.null(weights)) weights <- rep(1, length(Y))

  # make HAL design matrix
  basis_list <- enumerate_basis(
    x = as.matrix(W),
    max_degree = enumerate_basis_args$max_degree,
    smoothness_orders = enumerate_basis_args$smoothness_orders,
    num_knots = enumerate_basis_args$num_knots)
  phi_W <- make_design_matrix(X = as.matrix(W), blist = basis_list)

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
    basis_list <- basis_list[non_zero]
    beta <- coef(tau_fit, s = "lambda.min")[c(1, non_zero+1)]
  } else {
    non_zero <- which(coef(tau_fit, s = "lambda.min")[-c(1, 2)] != 0)
    basis_list <- basis_list[non_zero]
    beta <- coef(tau_fit, s = "lambda.min")[c(2, non_zero+2)]
  }
  phi_W <- as.matrix(cbind(1, phi_W[, non_zero, drop=FALSE]))
  pred <- as.vector(phi_W %*% beta)

  ret_obj = list(tau_fit = tau_fit,
                 beta = beta,
                 phi_train = phi_W,
                 basis_list = basis_list,
                 pred = pred,
                 pseudo_outcome = pseudo_outcome,
                 pseudo_weights = pseudo_weights)
  class(ret_obj) <- "rHAL"
  return(ret_obj)
}

predict.rHAL <- function(object,
                         newx = NULL) {
  if (!is.null(newx)) {
    newx_hal <- make_design_matrix(X = as.matrix(newx), blist = object$basis_list)
    return(as.vector(cbind(1, newx_hal) %*% object$beta))
  } else {
    pred <- return(object$pred)
  }
}
