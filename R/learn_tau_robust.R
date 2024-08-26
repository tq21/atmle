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
learn_tau_robust <- function(S,
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

    clever_cov <- (S*(2*A-1))/(S-Pi$pred)
    clever_cov_pred_cate <- (2*A-1)/(1-Pi$pred)
    clever_cov_A1 <- S/(S-Pi$A1)
    clever_cov_pred_cate_A1 <- 1/(1-Pi$A1)
    clever_cov_A0 <- -S/(S-Pi$A0)
    clever_cov_pred_cate_A0 <- -1/(1-Pi$A0)

    # design matrix
    X <- data.frame(1, clever_cov, W_aug, A, W_aug * A)
    X_pred_cate <- data.frame(1, clever_cov_pred_cate, W_aug, 1, W_aug)

    # counterfactual design matrices
    X_A1_counter <- data.frame(1, clever_cov_A1, W_aug, 1, W_aug)
    X_A1_counter_pred_cate <- data.frame(1, clever_cov_pred_cate_A1, W_aug, 1, W_aug)
    X_A0_counter <- data.frame(1, clever_cov_A0, W_aug, 0, W_aug * 0)
    X_A0_counter_pred_cate <- data.frame(1, clever_cov_pred_cate_A0, W_aug, 0, W_aug * 0)
  }

  if (!is.null(bias_working_model_formula)) {
    cov_only_formula <- as.formula(paste0("~ ", bias_working_model_formula))
    data_formula <- as.formula(paste0("~ ", bias_working_model_formula, " + Y"))
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

    if (controls_only) {
      x_basis <- x_basis_A0 <- as.matrix(model.matrix(cov_only_formula, data = X))
      pred <- A0 <- as.numeric(x_basis_A0 %*% matrix(coefs))
    } else {
      x_basis <- as.matrix(model.matrix(cov_only_formula, data = X))
      x_basis_A1 <- as.matrix(model.matrix(cov_only_formula, data = X_A1_counter))
      x_basis_A0 <- as.matrix(model.matrix(cov_only_formula, data = X_A0_counter))

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
      # penalized fit
      fit <- cv.glmnet(
        x = as.matrix(X[delta == 1, -1, drop = FALSE]),
        y = pseudo_outcome[delta == 1], intercept = TRUE,
        family = "gaussian", weights = pseudo_weights[delta == 1],
        keep = TRUE, nfolds = v_folds, alpha = 1
      )
    }
    non_zero <- which(as.numeric(coef(fit, s = "lambda.min")) != 0)
    non_zero <- sort(union(2, non_zero))
    x_basis <- X[, non_zero, drop = FALSE]
    x_basis_pred_cate <- X_pred_cate[, non_zero, drop = FALSE]
    x_basis_A1 <- X_A1_counter[, non_zero, drop = FALSE]
    x_basis_pred_cate_A1 <- X_A1_counter_pred_cate[, non_zero, drop = FALSE]
    x_basis_A0 <- X_A0_counter[, non_zero, drop = FALSE]
    x_basis_pred_cate_A0 <- X_A0_counter_pred_cate[, non_zero, drop = FALSE]

    # relaxed fit
    relax_fit <- glm(pseudo_outcome[delta == 1] ~ .,
                     weights = pseudo_weights[delta == 1],
                     data = x_basis[, -1, drop = FALSE],
                     family = "gaussian")
    coefs <- as.numeric(coef(relax_fit))
    na_idx <- which(is.na(coefs))
    if (length(na_idx) > 0) {
      coefs <- coefs[-na_idx]
      x_basis <- x_basis[, -na_idx, drop = FALSE]
      x_basis_pred_cate <- x_basis_pred_cate[, -na_idx, drop = FALSE]
      x_basis_A1 <- x_basis_A1[, -na_idx, drop = FALSE]
      x_basis_pred_cate_A1 <- x_basis_pred_cate_A1[, -na_idx, drop = FALSE]
      x_basis_A0 <- x_basis_A0[, -na_idx, drop = FALSE]
      x_basis_pred_cate_A0 <- x_basis_pred_cate_A0[, -na_idx, drop = FALSE]
    }

    if (controls_only) {
      # selected counterfactual bases
      x_basis <- x_basis_A0 <- as.matrix(cbind(1, W_aug)[, non_zero, drop = FALSE])

      # predictions
      pred <- A0 <- as.numeric(x_basis_A0 %*% matrix(coefs))
    } else {
      # predictions
      A1 <- as.numeric(as.matrix(x_basis_pred_cate_A1) %*% matrix(coefs))
      A0 <- as.numeric(as.matrix(x_basis_pred_cate_A0) %*% matrix(coefs))
      pred[A == 1] <- A1[A == 1]
      pred[A == 0] <- A0[A == 0]
    }
  } else if (method == "HAL") {
    non_zero <- NULL
    if (controls_only) {
      # external data has only controls
      X <- data.frame(W[A == 0, , drop = FALSE])

      if (min_working_model) {
        # enforce a minimal (main-term only) working model
        if (max_degree > 1) {
          aug_formula <- as.formula("~-1+(.)^" %+% max_degree)
        } else {
          aug_formula <- as.formula("~-1+(.)")
        }
        X_aug <- X_unpenalized <- model.matrix(aug_formula, data = data.frame(W))
      } else {
        # no minimal working model
        X_unpenalized <- X_unpenalized_A0 <- X_unpenalized_A1 <- NULL
      }

      # fit relaxed HAL
      fit <- fit_relaxed_hal(
        X = X[delta == 1, ], Y = pseudo_outcome[delta == 1],
        X_unpenalized = X_unpenalized,
        family = "gaussian", # ALWAYS GAUSSIAN, EVEN FOR BINARIES!
        weights = pseudo_weights[delta == 1],
        relaxed = TRUE
      )

      # selected bases
      x_basis <- x_basis_A0 <- make_counter_design_matrix(
        basis_list = fit$hal_basis_list,
        X_counterfactual = as.matrix(W),
        X_unpenalized = X_unpenalized
      )

      # predictions
      pred <- A0 <- as.numeric(x_basis_A0 %*% matrix(fit$beta))
    } else {
      # external data has both controls and treated
      X <- data.frame(W, A)
      X_A0 <- data.frame(W, A = 0)
      X_A1 <- data.frame(W, A = 1)

      # enforcing a minimal working model --------------------------------------
      if (min_working_model) {
        if (max_degree > 1) {
          aug_formula <- as.formula("~-1+(.)^" %+% max_degree %+% "+(A:.)^" %+% max_degree)
        } else {
          aug_formula <- as.formula("~-1+(.)+(A:.)")
        }
        X_aug <- X_unpenalized <- model.matrix(aug_formula, data = data.frame(W, A = A))
        X_A0_aug <- X_unpenalized_A0 <- model.matrix(aug_formula, data = data.frame(W, A = 0))
        X_A1_aug <- X_unpenalized_A1 <- model.matrix(aug_formula, data = data.frame(W, A = 1))
      } else {
        # no minimal working model
        X_unpenalized <- X_unpenalized_A0 <- X_unpenalized_A1 <- NULL
      }

      # fit relaxed HAL
      fit <- fit_relaxed_hal(
        X = X[delta == 1, ], Y = pseudo_outcome[delta == 1],
        X_unpenalized = X_unpenalized,
        family = "gaussian", # ALWAYS GAUSSIAN, EVEN FOR BINARIES!
        weights = pseudo_weights[delta == 1],
        relaxed = TRUE,
        enumerate_basis_args = enumerate_basis_args,
        fit_hal_args = fit_hal_args
      )

      # selected bases
      x_basis <- fit$x_basis
      x_basis_A1 <- make_counter_design_matrix(
        basis_list = fit$hal_basis_list,
        X_counterfactual = as.matrix(X_A1),
        X_unpenalized = X_unpenalized_A1[, fit$selected_unpenalized_idx, drop = FALSE]
      )
      x_basis_A0 <- make_counter_design_matrix(
        basis_list = fit$hal_basis_list,
        X_counterfactual = as.matrix(X_A0),
        X_unpenalized = X_unpenalized_A0[, fit$selected_unpenalized_idx, drop = FALSE]
      )

      # predictions
      A1 <- as.numeric(x_basis_A1 %*% matrix(fit$beta))
      A0 <- as.numeric(x_basis_A0 %*% matrix(fit$beta))
      pred[A == 1] <- A1[A == 1]
      pred[A == 0] <- A0[A == 0]
    }
  } else {
    stop("bias_working_model must be either 'glmnet' or 'HAL'")
  }

  return(list(
    A1 = A1,
    A0 = A0,
    x_basis = as.matrix(x_basis),
    x_basis_pred_cate = as.matrix(x_basis_pred_cate),
    x_basis_A1 = as.matrix(x_basis_A1),
    x_basis_pred_cate_A1 = as.matrix(x_basis_pred_cate_A1),
    x_basis_A0 = as.matrix(x_basis_A0),
    x_basis_pred_cate_A0 = as.matrix(x_basis_pred_cate_A0),
    pred = pred,
    coefs = coefs,
    non_zero = non_zero,
    pseudo_outcome = pseudo_outcome,
    pseudo_weights = pseudo_weights
  ))
}
