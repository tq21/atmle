# Code is adapted from rlearner:
# R-learner for Quasi-Oracle Estimation of Heterogeneous Treatment Effects

sim_data <- function(n,
                     counter_A = NULL) {
  # error
  UY <- rnorm(n, 0, 0.2)

  # baseline covariates
  W1 <- round(runif(n, -1, 1), 3)
  W2 <- round(runif(n, -1, 1), 3)
  W3 <- round(runif(n, -1, 1), 3)

  # treatment
  if (is.null(counter_A)) {
    A <- rbinom(n, 1, plogis(-0.25*W1+0.5*W2))
  } else {
    A <- rep(counter_A, n)
  }

  # outcome
  cate <- 0.5*as.numeric(W1 >= 0.2)*(W1-0.2)+
    1.1*as.numeric(W2 >= 0.5)*(W2-0.5)+
    2.1*as.numeric(W1 >= -0.5)*(W2+0.5)*as.numeric(W3 >= 0)*W3+
    1.4*as.numeric(W2 >= 0.9)*(W2-0.9)*as.numeric(W3 >= 0.1)*(W3-0.1)
  Y <- 1.9-0.4*A+0.9*W1+1.4*W2+2.1*W3+A*cate+UY

  data <- data.frame(W1 = W1,
                     W2 = W2,
                     W3 = W3,
                     A = A,
                     Y = Y)

  return(data)
}

get_truth <- function() {
  data_A1 <- sim_data(1e7, counter_A = 1)
  data_A0 <- sim_data(1e7, counter_A = 0)
  return(mean(data_A1$Y- data_A0$Y))
}

theta_lasso <- function(W,
                        Y,
                        foldid,
                        family) {

  # fit theta=E(Y|W) and obtain cross-fitted predictions
  theta_fit <- cv.glmnet(x = as.matrix(W),
                         y = Y,
                         foldid = foldid,
                         keep = TRUE,
                         alpha = 1,
                         family = family)
  non_na_idx <- !is.na(colSums(theta_fit$fit.preval))
  lambda_min <- theta_fit$lambda[which.min(theta_fit$cvm[non_na_idx])]
  pred <- theta_fit$fit.preval[, non_na_idx][, theta_fit$lambda[non_na_idx] == lambda_min]
  if (family == "binomial") {
    pred <- plogis(pred)
  }

  return(pred)
}

g_lasso <- function(W,
                    A,
                    foldid) {

  # fit g1W=P(A=1|W) and obtain cross-fitted predictions
  g_fit <- cv.glmnet(x = as.matrix(W),
                     y = A,
                     foldid = foldid,
                     family = "binomial",
                     keep = TRUE,
                     alpha = 1)
  non_na_idx <- !is.na(colSums(g_fit$fit.preval))
  lambda_min <- g_fit$lambda[which.min(g_fit$cvm[non_na_idx])]
  pred <- g_fit$fit.preval[, non_na_idx][, g_fit$lambda[non_na_idx] == lambda_min]
  pred <- plogis(pred)

  return(pred)
}

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
  pred <- as.vector(cbind(1, phi_W[, non_zero, drop=FALSE]) %*% beta)

  ret_obj = list(tau_fit = tau_fit,
                 beta = beta,
                 basis_list = basis_list,
                 pred = pred)
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
