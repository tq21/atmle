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
                      controls_only,
                      method,
                      v_folds,
                      max_degree,
                      min_working_model,
                      min_working_model_screen,
                      undersmooth,
                      target_gwt,
                      Pi_bounds,
                      enumerate_basis_args = list(),
                      fit_hal_args = list(),
                      weak_penalize = FALSE) {

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

  weights <- 1 # TODO: not used yet
  pseudo_outcome <- NULL
  pseudo_weights <- NULL

  if (max_degree > ncol(W)) {
    max_degree <- ncol(W)
  }

  if (controls_only) {
    pred <- numeric(length = length(A))

    # R-transformations, only controls
    pseudo_outcome <- (Y[A == 0]-theta[A == 0])/(S[A == 0]-Pi$pred[A == 0])
    pseudo_weights <- (S[A == 0]-Pi$pred[A == 0])^2*weights

    # design matrix, only controls (X: W)
    if (max_degree > 1) {
      W_aug <- model.matrix(as.formula(paste0("~ -1+(.)^", max_degree)), data = W)
      X <- W_aug[A == 0, ]

      # counterfactual design matrices
      X_A1_counter <- X_A0_counter <- cbind(1, W_aug)
    } else {
      X <- W[A == 0, ]

      # counterfactual design matrices
      X_A1_counter <- X_A0_counter <- cbind(1, W)
    }

  } else {
    pred <- numeric(length = length(A))

    # R-transformations
    pseudo_outcome <- (Y-theta)/(S-Pi$pred)
    pseudo_weights <- (S-Pi$pred)^2*weights

    # design matrix
    # if (max_degree > 1) {
    #   W_aug <- model.matrix(as.formula(paste0("~ -1+(.)^", max_degree)), data = W)
    #   X <- cbind(W_aug, A, W_aug*A)
    # } else {
    #   X <- cbind(W, A, W*A) # good
    #   #X <- cbind(W*A, A, W*(1-A))
    # }

    X <- cbind(W, A, W*A)

    # counterfactual design matrices
    X_A1_counter <- cbind(1, W, 1, W)
    X_A0_counter <- cbind(1, W, 0, W*0)
    #X_A1_counter <- cbind(1, W*1, 1, W*0)
    #X_A0_counter <- cbind(1, W*0, 0, W*1)
  }

  if (method == "binomial loss") {
    # TODO: testing stage right now. Do not use this.
    X <- S * cbind(W, A, W*A)
    fit <- cv.glmnet(x = as.matrix(X), y = Y, offset = theta,
                     intercept = TRUE, family = "binomial", keep = TRUE,
                     nfolds = v_folds, alpha = 1, relax = TRUE)
    non_zero <- which(as.numeric(coef(fit, s = "lambda.min", gamma = 0)) != 0)
    coefs <- coef(fit, s = "lambda.min", gamma = 0)[non_zero]

    X_A1_counter <- cbind(1, X * as.numeric(A == 1))
    X_A0_counter <- cbind(1, X * as.numeric(A == 0))

    # TODO: controls only
    # selected counterfactual bases
    x_basis <- as.matrix(cbind(1, X)[, non_zero, drop = FALSE])
    x_basis_A1 <- as.matrix(X_A1_counter[, non_zero, drop = FALSE])
    x_basis_A0 <- as.matrix(X_A0_counter[, non_zero, drop = FALSE])

    # predictions
    A1 <- as.numeric(x_basis_A1 %*% matrix(coefs))
    A0 <- as.numeric(x_basis_A0 %*% matrix(coefs))
    pred[A == 1] <- A1[A == 1]
    pred[A == 0] <- A0[A == 0]

  } else if (method == "glmnet") {
    fit <- cv.glmnet(x = as.matrix(X), y = pseudo_outcome, intercept = TRUE,
                     family = "gaussian", weights = pseudo_weights,
                     keep = TRUE, nfolds = v_folds, alpha = 1, relax = TRUE)
    non_zero <- which(as.numeric(coef(fit, s = "lambda.min", gamma = 0)) != 0)
    coefs <- coef(fit, s = "lambda.min", gamma = 0)[non_zero]
    #print(non_zero)

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

  } else if (method == "HAL") {
    non_zero <- NULL
    if (controls_only) {
      # external data has only controls
      X <- data.frame(W[A == 0, , drop = FALSE])

      ##########################################################################
      # EXPERIMENTAL FEATURE, NOT USED RIGHT NOW!
      if (min_working_model) {
        # enforce a minimal (main-term only) working model
        if (max_degree > 1) {
          aug_formula <- as.formula("~-1+(.)^" %+% max_degree)
        } else {
          aug_formula <- as.formula("~-1+(.)")
        }
        X_aug <- X_unpenalized <- model.matrix(aug_formula, data = data.frame(W))

        if (min_working_model_screen) {
          # use lasso to screen bases for minimal working model
          lasso_screen <- cv.glmnet(x = as.matrix(X_aug), y = pseudo_outcome, intercept = TRUE,
                                    family = "gaussian", weights = pseudo_weights,
                                    keep = TRUE, nfolds = v_folds, alpha = 1)
          X_unpenalized_idx <- which(as.numeric(coef(lasso_screen, s = "lambda.min"))[-1] != 0)
          if (length(X_unpenalized_idx) == 0) {
            X_unpenalized <- X_unpenalized_A0 <- X_unpenalized_A1 <- NULL
          } else {
            X_unpenalized <- X_unpenalized_A0 <- X_unpenalized_A1 <- X_aug[, X_unpenalized_idx, drop = FALSE]
          }
        }

      } else {
        # no minimal working model
        X_unpenalized <- X_unpenalized_A0 <- X_unpenalized_A1 <- NULL
      }

      # fit relaxed HAL
      fit <- fit_relaxed_hal(X = X, Y = pseudo_outcome,
                             X_unpenalized = X_unpenalized,
                             X_weak_penalized = NULL,
                             X_weak_penalized_level = 0,
                             family = "gaussian", # ALWAYS GAUSSIAN, EVEN FOR BINARIES!
                             weights = pseudo_weights,
                             relaxed = TRUE,
                             v_folds = v_folds)

      # selected bases
      x_basis <- x_basis_A0 <- make_counter_design_matrix(
        basis_list = fit$hal_basis_list,
        X_counterfactual = as.matrix(W),
        X_unpenalized = X_unpenalized)

      # predictions
      pred <- A0 <- as.numeric(x_basis_A0 %*% matrix(fit$beta))

    } else {
      # external data has both controls and treated
      X <- data.frame(W, A)
      X_A0 <- data.frame(W, A = 0)
      X_A1 <- data.frame(W, A = 1)

      ##########################################################################
      # EXPERIMENTAL FEATURE, NOT USED RIGHT NOW!
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

        if (min_working_model_screen) {
          # use lasso to screen bases for minimal working model
          lasso_screen <- cv.glmnet(x = as.matrix(X_aug), y = pseudo_outcome, intercept = TRUE,
                                    family = "gaussian", weights = pseudo_weights,
                                    keep = TRUE, nfolds = v_folds, alpha = 1)
          X_unpenalized_idx <- which(as.numeric(coef(lasso_screen, s = "lambda.min"))[-1] != 0)
          if (length(X_unpenalized_idx) == 0) {
            X_unpenalized <- X_unpenalized_A0 <- X_unpenalized_A1 <- NULL
          } else {
            X_unpenalized <- X_aug[, X_unpenalized_idx, drop = FALSE]
            X_unpenalized_A0 <- X_A0_aug[, X_unpenalized_idx, drop = FALSE]
            X_unpenalized_A1 <- X_A1_aug[, X_unpenalized_idx, drop = FALSE]
          }
        }

      } else {
        # no minimal working model
        X_unpenalized <- X_unpenalized_A0 <- X_unpenalized_A1 <- NULL
      }

      # weakly penalize terms --------------------------------------------------
      weak_penalize <- FALSE
      if (weak_penalize) {
        X_weak_penalized <- model.matrix(as.formula("~-1+.+A:."), data = X)
      } else {
        X_weak_penalized <- NULL
      }
      # EXPERIMENTAL FEATURE, NOT USED RIGHT NOW!
      ##########################################################################

      # undersmoothing ---------------------------------------------------------
      # EXPERIMENTAL FEATURE, NOT USED RIGHT NOW!
      undersmooth <- 0
      if (undersmooth == 0) {
        # fit relaxed HAL, no undersmoothing at all
        fit <- fit_relaxed_hal(X = X, Y = pseudo_outcome,
                               X_unpenalized = X_unpenalized,
                               X_weak_penalized = X_weak_penalized,
                               X_weak_penalized_level = 0,
                               family = "gaussian", # ALWAYS GAUSSIAN, EVEN FOR BINARIES!
                               weights = pseudo_weights,
                               relaxed = TRUE,
                               v_folds = v_folds,
                               enumerate_basis_args = enumerate_basis_args,
                               fit_hal_args = fit_hal_args)

        # selected bases
        x_basis <- fit$x_basis
        x_basis_A1 <- make_counter_design_matrix(
          basis_list = fit$hal_basis_list,
          X_counterfactual = as.matrix(X_A1),
          X_unpenalized = X_unpenalized_A1[, fit$selected_unpenalized_idx, drop = FALSE])
        x_basis_A0 <- make_counter_design_matrix(
          basis_list = fit$hal_basis_list,
          X_counterfactual = as.matrix(X_A0),
          X_unpenalized = X_unpenalized_A0[, fit$selected_unpenalized_idx, drop = FALSE])

        # predictions
        A1 <- as.numeric(x_basis_A1 %*% matrix(fit$beta))
        A0 <- as.numeric(x_basis_A0 %*% matrix(fit$beta))
        pred[A == 1] <- A1[A == 1]
        pred[A == 0] <- A0[A == 0]

      } else if (undersmooth == 1) {
        # undersmoothing version 1:
        # relaxed HAL, undersmooth to increase the model,
        # get the largest model such that PnD* is still solved
        tau <- undersmooth_1(S = S, W = W, A = A, X = X, Y = pseudo_outcome,
                             Pi = Pi, theta = theta, g = g,
                             pseudo_weights = pseudo_weights,
                             X_unpenalized = X_unpenalized,
                             controls_only = controls_only,
                             v_folds = v_folds,
                             target_gwt = target_gwt,
                             Pi_bounds = Pi_bounds)
        return(tau)

      } else if (undersmooth == 2) {
        # undersmoothing version 2:
        # use regular HAL, undersmooth till PnD* is solved
        tau <- undersmooth_2(S = S, W = W, A = A, X = X, Y = pseudo_outcome,
                             Pi = Pi, theta = theta, g = g,
                             X_A1 = X_A1, X_A0 = X_A0,
                             pseudo_weights = pseudo_weights,
                             X_unpenalized = X_unpenalized,
                             controls_only = controls_only,
                             v_folds = v_folds,
                             target_gwt = target_gwt,
                             Pi_bounds = Pi_bounds)
        return(tau)

      } else if (undersmooth == 3 | undersmooth == 4) {
        # undersmoothing version 3:
        # use larger model for targeting Pi only,
        # use CV selected model for point estimation and inference

        # undersmoothing version 4:
        # the same as version 3, but use the larger model for inference

        # tau objects
        two_taus <- undersmooth_3(S = S, W = W, A = A, X = X, Y = pseudo_outcome,
                                  Pi = Pi, theta = theta, g = g,
                                  pseudo_weights = pseudo_weights,
                                  X_unpenalized = X_unpenalized,
                                  controls_only = controls_only,
                                  v_folds = v_folds,
                                  target_gwt = target_gwt,
                                  Pi_bounds = Pi_bounds)
        return(two_taus)
      }
    }
  }

  return(list(A1 = A1,
              A0 = A0,
              x_basis = x_basis,
              x_basis_A1 = x_basis_A1,
              x_basis_A0 = x_basis_A0,
              pred = pred,
              coefs = fit$beta,
              non_zero = non_zero,
              pseudo_outcome = pseudo_outcome,
              pseudo_weights = pseudo_weights))
}
