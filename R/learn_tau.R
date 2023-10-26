#' Learn working model for the conditional effect of the study indicator on the
#' outcome
#'
#' @description Function to learn the conditional effect of the study indicator
#' on the outcome,
#' \eqn{\tau(W,A)=\mathbb{E}(Y\mid S=1,W,A)-\mathbb{E}(Y\mid S=0,W,A)}.
#'
#' @export
#'
#' @importFrom glmnet cv.glmnet
#'
#' @param S A vector of study indicators, \eqn{S=1} for RCT, \eqn{S=0} for RWD.
#' @param W A matrix of baseline covariates.
#' @param A A vector of treatment indicators, \eqn{A=1} for treatment-arm,
#' \eqn{A=0} for control-arm.
#' @param Y A vector of outcomes.
#' @param Pi A vector of estimated trial enrollment probabilities,
#' \eqn{\Pi(W,A)=\mathbb{P}(S=1\mid W,A)}.
#' @param theta A vector of estimated conditional mean of outcome given
#' baseline covariates and treatment, \eqn{\theta(W,A)=\mathbb{E}(Y\mid W,A)}.
#' @param controls_only A logical indicating whether to learn only among
#' controls. This applies when the external data only has control-arm
#' observations.
#' @param method Working model type. Either \code{"glmnet"} for lasso-based
#' working model or \code{"HAL"} for highly adaptive lasso-based working model.
#' Default is \code{"glmnet"}.
#' @param v_folds A numeric of number of folds for cross-validation
#' (when necessary).
#'
#' @returns A \code{list} containing the following elements:
#' \item{A1}{A numeric vector of the estimated counterfactual conditional
#' effects under treatment;}
#' \item{A0}{A numeric vector of the estimated counterfactual conditional
#' effects under control;}
#' \item{x_basis}{A numeric matrix of the working model bases;}
#' \item{x_basis_A1}{A numeric matrix of the counterfactual working model bases
#' under treatment;}
#' \item{x_basis_A0}{A numeric matrix of the counterfactual working model bases
#' under control;}
#' \item{pred}{A numeric vector of estimated conditional effects;}
#' \item{coefs}{A numeric vector of the working model coefficients.}
learn_tau <- function(S, W, A, Y,
                      Pi, theta,
                      controls_only,
                      method,
                      v_folds) {

  # initialize
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
    pseudo_outcome <- (Y[A == 0]-theta[A == 0])/(S[A == 0]-Pi$pred[A == 0])
    pseudo_weights <- (S[A == 0]-Pi$pred[A == 0])^2*weights

    # design matrix, only controls (X: W)
    X <- W[A == 0, ]

    # counterfactual design matrices
    X_A1_counter <- X_A0_counter <- cbind(1, W)

  } else {
    pred <- numeric(length = length(A))

    # R-transformations
    pseudo_outcome <- (Y-theta)/(S-Pi$pred)
    pseudo_weights <- (S-Pi$pred)^2*weights

    # design matrix
    X <- cbind(W, A, W*A)

    # counterfactual design matrices
    X_A1_counter <- cbind(1, W, 1, W)
    X_A0_counter <- cbind(1, W, 0, W*0)
  }

  if (method == "glmnet") {
    fit <- cv.glmnet(x = as.matrix(X), y = pseudo_outcome, intercept = TRUE,
                     family = "gaussian", weights = pseudo_weights,
                     keep = TRUE, nfolds = v_folds, alpha = 1, relax = TRUE)
    non_zero <- which(as.numeric(coef(fit, s = "lambda.min", gamma = 0)) != 0)
    coefs <- coef(fit, s = "lambda.min", gamma = 0)[non_zero]

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
    #print("bases selected: " %+% length(coefs))

  } else if (method == "HAL") {
    # TODO: not fully working right now
    X <- NULL
    if (controls_only) {
      X <- as.matrix(W[A == 0, ])
    } else {
      X <- as.matrix(cbind(W, A))
    }

    # add in main terms + interactions, without penalization
    # X_unpenalized <- model.matrix(Y ~ -1 + (.)^2, data = data.frame(X))
    # X_unpenalized_A1 <- model.matrix(Y ~ -1 + (.)^2, data = cbind(W, A = 1))
    # X_unpenalized_A0 <- model.matrix(Y ~ -1 + (.)^2, data = cbind(W, A = 0))
    #X_unpenalized <- X
    # X_unpenalized <- X * (1-A)
    # X_unpenalized_A1 <- as.matrix(cbind(W, A = 1))
    # X_unpenalized_A0 <- as.matrix(cbind(W, A = 0))

    # X <- as.matrix(cbind(W[, c(1, 2)], A))
    # X_unpenalized <- model.matrix(Y ~ -1 + (.)^2, data = data.frame(W[, c(1, 2)] * (1-A)))

    # fit HAL
    fit <- fit_relaxed_hal(X = X, Y = pseudo_outcome,
                           family = "gaussian",
                           weights = pseudo_weights)

    if (controls_only) {
      # design matrices
      x_basis <- x_basis_A0 <- make_counter_design_matrix(fit$hal_basis_list, as.matrix(W))

      # predictions
      pred <- A0 <- as.numeric(x_basis_A0 %*% matrix(fit$beta))
    } else {
      # design matrices
      x_basis <- fit$x_basis
      x_basis_A1 <- make_counter_design_matrix(fit$hal_basis_list, as.matrix(cbind(W, A = 1)))
      x_basis_A0 <- make_counter_design_matrix(fit$hal_basis_list, as.matrix(cbind(W, A = 0)))

      # predictions
      A1 <- as.numeric(x_basis_A1 %*% matrix(fit$beta))
      A0 <- as.numeric(x_basis_A0 %*% matrix(fit$beta))
      pred[A == 1] <- A1[A == 1]
      pred[A == 0] <- A0[A == 0]
    }

    coefs <- fit$beta
    #print("bases selected: " %+% length(coefs))
  }

  return(list(A1 = A1,
              A0 = A0,
              x_basis = x_basis,
              x_basis_A1 = x_basis_A1,
              x_basis_A0 = x_basis_A0,
              pred = pred,
              coefs = coefs))
}
