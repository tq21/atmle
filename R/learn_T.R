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
                    max_degree,
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

    # check arguments
    enumerate_basis_default_args <- list(
      max_degree = 2,
      smoothness_orders = 1,
      num_knots = 20
    )
    enumerate_basis_args <- modifyList(
      enumerate_basis_default_args,
      enumerate_basis_args
    )

    # make data formula
    if (max_degree > 1) {
      aug_formula <- as.formula("~-1+(.)^" %+% max_degree)
    } else {
      aug_formula <- as.formula("~-1+(.)")
    }

    X <- model.matrix(aug_formula, data = data.frame(W))

    # make design matrix
    basis_list <- enumerate_basis(x = as.matrix(X[delta == 1, ]),
                                  max_degree = enumerate_basis_args$max_degree,
                                  smoothness_orders = enumerate_basis_args$smoothness_orders,
                                  num_knots = enumerate_basis_args$num_knots)
    X_hal <- make_design_matrix(X = as.matrix(X[delta == 1, ]),
                                blist = basis_list)

    # fit penalized HAL
    fit <- cv.glmnet(x = as.matrix(X_hal),
                     y = pseudo_outcome[delta == 1],
                     weights = pseudo_weights[delta == 1],
                     family = "gaussian",
                     alpha = 1,
                     nfolds = v_folds,
                     keep = TRUE)

    # non-zero bases
    non_zero <- which(as.numeric(coef(fit, s = "lambda.min"))[-1] != 0)
    basis_list <- basis_list[non_zero]
    X_hal_selected <- X_hal[, non_zero, drop = FALSE]

    if (length(non_zero) > 0) {
      # fit relaxed HAL
      fit <- glm(pseudo_outcome[delta == 1] ~ .,
                 family = "gaussian",
                 data = data.frame(as.matrix(X_hal_selected)),
                 weights = pseudo_weights[delta == 1])
      coefs <- as.numeric(coef(fit))
      na_idx <- which(is.na(coefs[-1]))
      if (length(na_idx) > 0) {
        coefs <- coefs[!is.na(coefs)]
        basis_list <- basis_list[-na_idx]
        X_hal_selected <- X_hal_selected[, -na_idx, drop = FALSE]
      }
    } else {
      coefs <- mean(pseudo_outcome[delta == 1])
    }

    x_basis <- as.matrix(cbind(1, X_hal_selected))
    pred <- as.numeric(x_basis %*% matrix(coefs))
  }

  return(list(
    pred = pred,
    x_basis = x_basis,
    coefs = coefs,
    non_zero = non_zero,
    pseudo_outcome = pseudo_outcome,
    pseudo_weights = pseudo_weights,
    x_basis_all = as.matrix(cbind(1, W))
  ))
}
