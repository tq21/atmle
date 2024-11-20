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
learn_tau_seq <- function(S,
                          W,
                          A,
                          Y,
                          Pi,
                          theta,
                          g,
                          delta,
                          controls_only,
                          v_folds,
                          max_degree,
                          min_working_model,
                          enumerate_basis_args,
                          fit_hal_args,
                          weights,
                          bias_working_model_formula) {

  pred <- numeric(length = length(A))

  if (max_degree > ncol(W)) {
    max_degree <- ncol(W)
  }

  if (controls_only) {
    pseudo_outcome <- (Y[A == 0] - theta[A == 0]) / (S[A == 0] - Pi$pred[A == 0])
    pseudo_weights <- (S[A == 0] - Pi$pred[A == 0])^2 * weights[A == 0]
    if (max_degree > 1) {
      W_aug <- model.matrix(as.formula(paste0("~ -1+(.)^", max_degree)), data = W)
    } else {
      W_aug <- W
    }
    X <- W_aug[A == 0, ]
    X_A1_counter <- X_A0_counter <- cbind(1, W_aug)

  } else {
    pseudo_outcome <- (Y - theta) / (S - Pi$pred)
    pseudo_weights <- (S - Pi$pred)^2 * weights
    if (max_degree > 1) {
      W_aug <- model.matrix(as.formula(paste0("~ -1+(.)^", max_degree)), data = W)
    } else {
      W_aug <- W
    }
    X <- data.frame(W_aug, A, W_aug * A)
    X_A1_counter <- data.frame(1, W_aug, 1, W_aug)
    X_A0_counter <- data.frame(1, W_aug, 0, W_aug * 0)
  }

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

  # construct sequence of estimators
  lambda_seq <- fit$lambda
  res <- map(lambda_seq, function(.s) {
    pred <- numeric(nrow(X_hal))

    # non-zero bases
    non_zero <- which(as.numeric(coef(fit, s = .s))[-1] != 0)
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

    return(list(A1 = A1,
                A0 = A0,
                x_basis = x_basis,
                x_basis_A1 = x_basis_A1,
                x_basis_A0 = x_basis_A0,
                pred = pred,
                coefs = coefs,
                non_zero = non_zero))
  })

  return(res)
}
