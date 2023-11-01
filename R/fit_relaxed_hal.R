#' @title Fit relaxed HAL
#'
#' @param X A numeric matrix of covariates.
#' @param Y A numeric vector of outcomes.
#' @param family A character string specifying the family of the outcome.
#' @param weights A numeric vector of weights.
#' @param X_no_basis A numeric matrix of covariates that should not be used
#' to construct HAL basis.
#' @param X_unpenalized A numeric matrix of covariates that should not be
#' penalized by HAL.
#'
#' @return A list with the following components:
#' \item{beta}{A numeric vector of coefficients.}
#' \item{pred}{A numeric vector of predictions.}
#' \item{hal_basis_list}{A list of HAL basis.}
#' \item{x_basis}{A numeric matrix of basis.}
fit_relaxed_hal <- function(X, Y,
                            family,
                            weights = NULL,
                            X_no_basis = NULL,
                            X_unpenalized = NULL) {
  x_basis <- NULL
  hal_basis_list <- NULL

  hal_fit <- fit_hal(X = X,
                     Y = Y,
                     family = family,
                     weights = weights,
                     max_degree = 2,
                     X_no_basis = X_no_basis,
                     X_unpenalized = X_unpenalized,
                     smoothness_orders = 1,
                     return_x_basis = TRUE,
                     fit_control = list(use_min = FALSE))
  hal_basis_list <- list()

  if (!is.null(X_unpenalized)) {
    hal_bases_idx <- 2:(length(hal_fit$coefs)-ncol(X_unpenalized))
    hal_basis_list <- hal_fit$basis_list[hal_fit$coefs[hal_bases_idx] != 0]
  } else {
    hal_basis_list <- hal_fit$basis_list[hal_fit$coefs[-1] != 0]
  }

  # if (!is.null(X_no_basis)) {
  #   hal_bases_idx <- 1:(length(hal_fit$coefs)-ncol(X_no_basis)-1)
  #   hal_basis_list <- hal_fit$basis_list[hal_fit$coefs[hal_bases_coefs] != 0]
  # } else {
  #   hal_basis_list <- hal_fit$basis_list[hal_fit$coefs[-1] != 0]
  # }

  selected_idx <- which(hal_fit$coefs != 0)
  x_basis <- cbind(1, hal_fit$x_basis)[, selected_idx, drop = FALSE]

  # basis list empty, intercept only
  # if (length(hal_basis_list) == 0) {
  #   x_basis <- matrix(rep(1, length(Y)))
  # }

  # relaxed HAL fit
  beta <- NULL
  pred <- NULL
  if (family == "binomial") {
    hal_relaxed_fit <- glm.fit(x = x_basis, y = Y, family = binomial(), weights = weights, intercept = FALSE)
    beta <- coef(hal_relaxed_fit)
    beta[is.na(beta)] <- 0
    pred <- x_basis %*% beta
    pred <- as.vector(1 / (1 + exp(-pred)))
  } else if (family == "gaussian") {
    hal_relaxed_fit <- glm.fit(x = x_basis, y = Y, family = gaussian(), weights = weights, intercept = FALSE)
    beta <- coef(hal_relaxed_fit)
    beta[is.na(beta)] <- 0
    pred <- as.vector(x_basis %*% beta)
  }

  return(list(beta = beta,
              hal_basis_list = hal_basis_list,
              pred = pred,
              x_basis = as.matrix(x_basis)))
}
