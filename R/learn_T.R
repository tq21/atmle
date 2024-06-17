#' @title Learn working model for the conditional effect of treatment on the
#' outcome
#'
#' @description Function to learn the conditional effect of treatment on the
#' outcome, \eqn{T(W)=\mathbb{E}(Y\mid W,A=1)-\mathbb{E}(Y\mid W,A=0)}.
#'
#' @keywords working model
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
learn_T <- function(data,
                    W,
                    A,
                    Y,
                    delta,
                    g,
                    theta_tilde,
                    weights,
                    method,
                    folds,
                    enumerate_basis_args,
                    fit_hal_args) {

  # R-transformations
  pseudo_outcome <- ifelse(abs(data[[A]] - g) < 1e-10, 0, (data[[Y]] - theta_tilde) / (data[[A]] - g))
  pseudo_weights <- (data[[A]] - g)^2 * weights

  if (method == "glmnet") {
    # learn working model using main-term glmnet
    fit <- cv.glmnet(x = as.matrix(data[delta == 1, ..W]),
                     y = pseudo_outcome[delta == 1],
                     family = "gaussian",
                     weights = pseudo_weights[delta == 1], keep = TRUE,
                     nfolds = length(folds), alpha = 1, relax = TRUE)
    non_zero <- which(as.numeric(coef(fit, s = "lambda.min", gamma = 0)) != 0)
    coefs <- coef(fit, s = "lambda.min", gamma = 0)[non_zero]
    x_basis <- as.matrix(as.data.frame(cbind(1, data[, ..W]))[, non_zero, drop = FALSE])
    pred <- as.numeric(x_basis %*% matrix(coefs))

  } else if (method == "HAL") {
    # learn working model using HAL
    fit <- fit_relaxed_hal(X = data[delta == 1, ..W],
                           Y = pseudo_outcome[delta == 1],
                           family = "gaussian",
                           weights = pseudo_weights[delta == 1],
                           enumerate_basis_args = enumerate_basis_args,
                           fit_hal_args = fit_hal_args)
    x_basis <- fit$x_basis
    coefs <- fit$beta
    pred <- as.numeric(x_basis %*% matrix(coefs))
  }

  return(list(pred = pred,
              x_basis = x_basis,
              coefs = coefs))
}
