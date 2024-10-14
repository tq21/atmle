#' @title Learn working model for the conditional effect of the study indicator \
#' on the outcome
#'
#' @description Function to learn the conditional effect of the study indicator
#' on the outcome,
#' \eqn{\tau(W,A)=\mathbb{E}(Y\mid S=1,W,A)-\mathbb{E}(Y\mid S=0,W,A)}.
#'
#' @keywords working model
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
#' @param controls_only A logical indicating whether the external data has only
#' control-arm observations.
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
learn_tau <- function(S,
                      W,
                      A,
                      Y,
                      Pi,
                      theta,
                      g,
                      delta,
                      controls_only,
                      method,
                      v_folds,
                      max_degree,
                      min_working_model,
                      target_gwt,
                      Pi_bounds,
                      enumerate_basis_args,
                      fit_hal_args,
                      weights,
                      bias_working_model_formula) {

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
  non_zero <- NULL

  pseudo_outcome <- NULL
  pseudo_weights <- NULL

  if (max_degree > ncol(W)) {
    max_degree <- ncol(W)
  }

  if (controls_only) {
    pred <- numeric(length = length(A))

    # R-transformations, only controls
    pseudo_outcome <- (Y[A == 0] - theta[A == 0]) / (S[A == 0] - Pi$pred[A == 0])
    pseudo_weights <- (S[A == 0] - Pi$pred[A == 0])^2 * weights[A == 0]

    # augment design matrix if needed
    if (max_degree > 1) {
      W_aug <- model.matrix(as.formula(paste0("~ -1+(.)^", max_degree)), data = W)
    } else {
      W_aug <- W
    }

    # design matrix
    X <- W_aug[A == 0, ]

    # counterfactual design matrices
    X_A1_counter <- X_A0_counter <- cbind(1, W_aug)
  } else {
    pred <- numeric(length = length(A))

    # R-transformations
    pseudo_outcome <- (Y - theta) / (S - Pi$pred)
    pseudo_weights <- (S - Pi$pred)^2 * weights

    # augment design matrix if needed
    if (max_degree > 1) {
      W_aug <- model.matrix(as.formula(paste0("~ -1+(.)^", max_degree)), data = W)
    } else {
      W_aug <- W
    }

    # design matrix
    X <- data.frame(W_aug, A, W_aug * A)

    # counterfactual design matrices
    X_A1_counter <- data.frame(1, W_aug, 1, W_aug)
    X_A0_counter <- data.frame(1, W_aug, 0, W_aug * 0)
  }

  if (!is.null(bias_working_model_formula)) {
    cov_only_formula <- as.formula(paste0("~ ", bias_working_model_formula))
    cov_only_no_intercept_formula <- as.formula(paste0("~ ", bias_working_model_formula, " - 1"))
    data_formula <- as.formula(paste0("~ -1 + ", bias_working_model_formula, " + Y"))
    train_formula <- as.formula(paste0("Y ~ ", bias_working_model_formula))
    if (controls_only) {
      df_train <- model.matrix(data_formula,
                               data = cbind(X[delta[A == 0] == 1, ],
                                            Y = pseudo_outcome[delta[A == 0] == 1]))
      fit <- glm(formula = train_formula,
                 family = "gaussian",
                 data = as.data.frame(df_train),
                 weights = pseudo_weights[delta[A == 0] == 1])
    } else {
      df_train <- model.matrix(data_formula,
                               data = cbind(X[delta == 1, ],
                                            Y = pseudo_outcome[delta == 1]))
      fit <- glm(formula = train_formula,
                 family = "gaussian",
                 data = as.data.frame(df_train),
                 weights = pseudo_weights[delta == 1])
    }

    coefs <- as.numeric(coef(fit))
    coefs[is.na(coefs)] <- 0

    if (controls_only) {
      x_basis <- x_basis_A0 <- as.matrix(model.matrix(cov_only_formula, data = X))
      pred <- A0 <- as.numeric(x_basis_A0 %*% matrix(coefs))
    } else {
      x_basis <- as.matrix(model.matrix(cov_only_formula, data = X))
      x_basis_A1 <- as.matrix(model.matrix(cov_only_no_intercept_formula, data = X_A1_counter))
      x_basis_A0 <- as.matrix(model.matrix(cov_only_no_intercept_formula, data = X_A0_counter))

      # predictions
      A1 <- as.numeric(x_basis_A1 %*% matrix(coefs))
      A0 <- as.numeric(x_basis_A0 %*% matrix(coefs))
      pred[A == 1] <- A1[A == 1]
      pred[A == 0] <- A0[A == 0]
    }

  } else if (method == "glmnet") {
    if (controls_only) {
      fit <- cv.glmnet(
        x = as.matrix(X[delta[A == 0] == 1, , drop = FALSE]),
        y = pseudo_outcome[delta[A == 0] == 1], intercept = TRUE,
        family = "gaussian", weights = pseudo_weights[delta[A == 0] == 1],
        keep = TRUE, nfolds = v_folds, alpha = 1, relax = TRUE
      )
    } else {
      fit <- cv.glmnet(
        x = as.matrix(X[delta == 1, , drop = FALSE]),
        y = pseudo_outcome[delta == 1], intercept = TRUE,
        family = "gaussian", weights = pseudo_weights[delta == 1],
        keep = TRUE, nfolds = v_folds, alpha = 1, relax = TRUE
      )
    }

    non_zero <- which(as.numeric(coef(fit, s = "lambda.min", gamma = 0)) != 0)
    coefs <- coef(fit, s = "lambda.min", gamma = 0)[non_zero]

    if (controls_only) {
      # selected counterfactual bases
      x_basis <- x_basis_A0 <- as.matrix(cbind(1, W_aug)[, non_zero, drop = FALSE])

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
  } else if (method == "HAL") {

    browser()

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

    if (controls_only) {
      X <- data.frame(W[A == 0, , drop = FALSE])
      X_A0 <- data.frame(W)
      X <- model.matrix(aug_formula, data = X)
      X_A0_counter <- model.matrix(aug_formula, data = X_A0)
    } else {
      X <- data.frame(W, A)
      X_A0 <- data.frame(W, A = 0)
      X_A1 <- data.frame(W, A = 1)
      X <- model.matrix(aug_formula, data = X)
      X_A0_counter <- model.matrix(aug_formula, data = X_A0)
      X_A1_counter <- model.matrix(aug_formula, data = X_A1)
    }

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

    # selected bases and predictions
    if (controls_only) {
      x_basis <- x_basis_A0 <- make_counter_design_matrix(
        basis_list = basis_list,
        X_counterfactual = X_A0_counter,
        X_unpenalized = NULL
      )
      pred <- A0 <- as.numeric(x_basis_A0 %*% matrix(coefs))
    } else {
      x_basis <- as.matrix(cbind(1, X_hal_selected))
      x_basis_A1 <- make_counter_design_matrix(
        basis_list = basis_list,
        X_counterfactual = as.matrix(X_A1_counter),
        X_unpenalized = NULL
      )
      x_basis_A0 <- make_counter_design_matrix(
        basis_list = basis_list,
        X_counterfactual = as.matrix(X_A0_counter),
        X_unpenalized = NULL
      )
      A1 <- as.numeric(x_basis_A1 %*% matrix(coefs))
      A0 <- as.numeric(x_basis_A0 %*% matrix(coefs))
      pred[A == 1] <- A1[A == 1]
      pred[A == 0] <- A0[A == 0]
    }

  } else {
    stop("bias_working_model must be either 'glmnet' or 'HAL'")
  }

  return(list(
    A1 = A1,
    A0 = A0,
    x_basis = x_basis,
    x_basis_A1 = x_basis_A1,
    x_basis_A0 = x_basis_A0,
    pred = pred,
    coefs = coefs,
    non_zero = non_zero,
    pseudo_outcome = pseudo_outcome,
    pseudo_weights = pseudo_weights
  ))
}
