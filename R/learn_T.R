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
learn_T <- function(W,
                    A,
                    Y,
                    g,
                    delta,
                    theta_tilde,
                    method,
                    min_working_model,
                    v_folds,
                    weights,
                    enumerate_basis_args,
                    fit_hal_args,
                    pooled_working_model_formula) {

  # R-transformations
  pseudo_outcome <- ifelse(abs(A - g) < 1e-10, 0, (Y - theta_tilde) / (A - g))
  pseudo_weights <- (A - g)^2 * weights

  pred <- NULL
  x_basis <- NULL
  coefs <- NULL
  non_zero <- NULL

  if (!is.null(pooled_working_model_formula)) {
    cov_only_formula <- as.formula(paste0("~ ", pooled_working_model_formula))
    data_formula <- as.formula(paste0("~ ", pooled_working_model_formula, " + Y"))
    train_formula <- as.formula(paste0("Y ~ ", pooled_working_model_formula))
    df_train <- model.matrix(data_formula,
                             data = cbind(W[delta == 1, , drop = FALSE],
                                          Y = pseudo_outcome[delta == 1]))
    fit <- glm(formula = train_formula,
               family = "gaussian",
               data = as.data.frame(df_train),
               weights = pseudo_weights[delta == 1])
    coefs <- as.numeric(coef(fit))
    x_basis <- as.matrix(model.matrix(cov_only_formula, data = W))
    pred <- as.numeric(x_basis %*% matrix(coefs))
  } else if (method == "glmnet") {
    fit <- cv.glmnet(
      x = as.matrix(W[delta == 1, ]), y = pseudo_outcome[delta == 1],
      family = "gaussian", weights = pseudo_weights[delta == 1],
      keep = TRUE, nfolds = v_folds, alpha = 1, relax = TRUE
    )
    non_zero <- which(as.numeric(coef(fit, s = "lambda.min", gamma = 0)) != 0)
    coefs <- coef(fit, s = "lambda.min", gamma = 0)[non_zero]
    x_basis <- as.matrix(cbind(1, W)[, non_zero, drop = FALSE])
    pred <- as.numeric(x_basis %*% matrix(coefs))
  } else if (method == "HAL") {
    X <- as.matrix(W)

    if (min_working_model) {
      X_unpenalized <- as.matrix(cbind(1, W))
    } else {
      X_unpenalized <- NULL
    }

    # fit HAL
    fit <- fit_relaxed_hal(
      X = X[delta == 1, ], Y = pseudo_outcome[delta == 1],
      X_unpenalized = X_unpenalized,
      family = "gaussian",
      weights = pseudo_weights[delta == 1],
      relaxed = TRUE,
      enumerate_basis_args = enumerate_basis_args,
      fit_hal_args = fit_hal_args
    )

    # design matrices
    x_basis <- fit$x_basis

    # predictions
    pred <- as.numeric(x_basis %*% matrix(fit$beta))
    coefs <- fit$beta
  }

  return(list(
    pred = pred,
    x_basis = x_basis,
    coefs = coefs,
    non_zero = non_zero,
    pseudo_outcome = pseudo_outcome,
    pseudo_weights = pseudo_weights
  ))
}
