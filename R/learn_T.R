#' Learn working model for the conditional effect of treatment on the outcome
#'
#' @description Function to learn the conditional effect of treatment on the
#' outcome, \eqn{T(W)=\mathbb{E}(Y\mid W,A=1)-\mathbb{E}(Y\mid W,A=0)}.
#'
#' @export
#'
#' @importFrom glmnet cv.glmnet
#' @importFrom data.table data.table
#'
#' @param W A matrix of baseline covariates.
#' @param A A vector of treatment indicators, \eqn{A=1} for treatment-arm,
#' \eqn{A=0} for control-arm.
#' @param Y A vector of outcomes.
#' @param g A vector of estimated treatment probabilities,
#' \eqn{g(A\mid W)=\mathbb{P}(A=1\mid W)}.
#' @param theta_tilde A vector of estimated conditional mean of outcome given
#' baseline covariates, \eqn{\theta(W)=\mathbb{E}(Y\mid W)}.
#' @param method Working model type. Either \code{"glmnet"} for lasso-based
#' working model or \code{"HAL"} for highly adaptive lasso-based working model.
#' Default is \code{"glmnet"}.
#' @param v_folds A numeric of number of folds for cross-validation
#' (when necessary).
#'
#' @returns A \code{list} containing the following elements:
#' \item{pred}{A numeric vector of the estimated conditional effects;}
#' \item{x_basis}{A numeric matrix of the working model bases;}
#' \item{coefs}{A numeric vector of the working model coefficients.}
learn_T <- function(W, A, Y,
                    g, theta_tilde,
                    method,
                    v_folds) {
  # R-transformations
  weights <- 1 # TODO: currently not used
  pseudo_outcome <- ifelse(abs(A-g) < 1e-5, 0, (Y-theta_tilde)/(A-g))
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
    # TODO: not fully working right now
    X <- as.matrix(W)

    # fit HAL
    fit <- fit_relaxed_hal(X = X, Y = pseudo_outcome,
                           family = "gaussian",
                           weights = pseudo_weights)

    # design matrices
    x_basis <- fit$x_basis

    # predictions
    pred <- as.numeric(x_basis %*% matrix(fit$beta))

    coefs <- fit$beta
    #print("bases selected: " %+% length(coefs))
  }

  return(list(pred = pred,
              x_basis = x_basis,
              coefs = coefs))
}
