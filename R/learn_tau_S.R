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
learn_tau_S <- function(S,
                        W,
                        A,
                        Y,
                        Pi,
                        Pi_star,
                        theta,
                        g1W,
                        delta,
                        controls_only,
                        method,
                        v_folds,
                        min_working_model,
                        target_gwt,
                        Pi_bounds,
                        enumerate_basis_args,
                        weights,
                        bias_working_model_formula,
                        target_method = "oneshot",
                        verbose = TRUE) {

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

    W_aug <- W

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
      # TODO: check
      fit <- cv.glmnet(x = as.matrix(X[delta[A == 0] == 1, , drop = FALSE]),
                       y = pseudo_outcome[delta[A == 0] == 1],
                       weights = pseudo_weights[delta[A == 0] == 1],
                       family = "gaussian", keep = TRUE, nfolds = v_folds,
                       alpha = 1)
      non_zero <- which(as.numeric(coef(fit, s = "lambda.min")) != 0)
      coefs <- coef(fit, s = "lambda.min")[non_zero]
      x_basis <- x_basis_A0 <- as.matrix(cbind(1, W_aug)[, non_zero, drop = FALSE])
      X_select <- as.matrix(cbind(1, X)[delta[A == 0] == 1, non_zero, drop = FALSE])

      # targeting
      if (length(non_zero) > 1) {
        if (target_method == "relaxed") {
          # use relaxed fit as targeted fit
          relaxed_fit <- glm(pseudo_outcome[delta[A == 0] == 1] ~ -1+.,
                             family = "gaussian",
                             data = data.frame(X_select[delta[A == 0] == 1,,drop=FALSE]),
                             weights = pseudo_weights[delta[A == 0] == 1])
          coefs <- as.numeric(coef(relaxed_fit))
          na_idx <- which(is.na(coefs))
          if (length(na_idx) > 0) {
            coefs <- coefs[!is.na(coefs)]
            x_basis <- x_basis[, -na_idx, drop=FALSE]
          }
          A0 <- as.numeric(x_basis_A0 %*% matrix(coefs))
          psi_pound_est <- mean((1-Pi$A0)*A0)
          eic <- get_eic_psi_pound(Pi = Pi,
                                   tau = list(A0 = A0,
                                              x_basis = x_basis,
                                              x_basis_A0 = x_basis_A0),
                                   g = g,
                                   theta = theta,
                                   psi_pound_est = psi_pound_est,
                                   S = S,
                                   A = A,
                                   Y = Y,
                                   n = length(Y),
                                   controls_only = controls_only,
                                   weights = weights)
        } else if (target_method == "onestep") {
          # TODO
        }
      } else {
        coefs <- mean(pseudo_outcome[delta == 1])
      }
    } else {
      fit <- cv.glmnet(x = as.matrix(X[delta == 1, , drop = FALSE]),
                       y = pseudo_outcome[delta == 1],
                       weights = pseudo_weights[delta == 1],
                       family = "gaussian", keep = TRUE, nfolds = v_folds,
                       alpha = 1)
      non_zero <- which(as.numeric(coef(fit, s = "lambda.min")) != 0)
      coefs <- coef(fit, s = "lambda.min")[non_zero]
      x_basis <- as.matrix(cbind(1, X)[, non_zero, drop = FALSE])
      x_basis_A0 <- as.matrix(X_A0_counter[, non_zero, drop = FALSE])
      x_basis_A1 <- as.matrix(X_A1_counter[, non_zero, drop = FALSE])

      # targeting
      if (length(non_zero) > 1) {
        if (target_method == "relaxed") {
          # use relaxed fit as targeted fit
          relaxed_fit <- glm(pseudo_outcome[delta == 1] ~ -1+.,
                             family = "gaussian",
                             data = data.frame(x_basis[delta == 1,,drop=FALSE]),
                             weights = pseudo_weights[delta == 1])
          coefs <- as.numeric(coef(relaxed_fit))
          na_idx <- which(is.na(coefs))
          if (length(na_idx) > 0) {
            coefs <- coefs[!is.na(coefs)]
            x_basis <- x_basis[, -na_idx, drop=FALSE]
            x_basis_A1 <- x_basis_A1[, -na_idx, drop=FALSE]
            x_basis_A0 <- x_basis_A0[, -na_idx, drop=FALSE]
          }
          A0 <- as.numeric(x_basis_A0 %*% matrix(coefs))
          A1 <- as.numeric(x_basis_A1 %*% matrix(coefs))
          pred <- as.numeric(x_basis %*% matrix(coefs))
          psi_pound_est <- mean((1-Pi$A0)*A0-(1-Pi$A1)*A1)
          eic <- get_eic_psi_pound(Pi = Pi,
                                   tau = list(pred = pred,
                                              A0 = A0,
                                              A1 = A1,
                                              x_basis = x_basis,
                                              x_basis_A0 = x_basis_A0,
                                              x_basis_A1 = x_basis_A1),
                                   g = g,
                                   theta = theta,
                                   psi_pound_est = psi_pound_est,
                                   S = S,
                                   A = A,
                                   Y = Y,
                                   n = length(Y),
                                   controls_only = controls_only,
                                   weights = weights)
        } else if (target_method == "onestep") {
          # TODO: closed-form
          direction <- get_beta_h(x_basis = x_basis,
                                  x_basis_A0 = x_basis_A0,
                                  x_basis_A1 = x_basis_A1,
                                  Pi = Pi)
          n <- length(Y)
          sn <- 0
          cur_iter <- 1
          A0 <- as.numeric(x_basis_A0 %*% matrix(coefs))
          A1 <- as.numeric(x_basis_A1 %*% matrix(coefs))
          pred <- as.numeric(x_basis %*% matrix(coefs))
          psi_pound_est <- mean((1-Pi$A0)*A0-(1-Pi$A1)*A1)
          eic <- get_eic_psi_pound(Pi = Pi,
                                   tau = list(pred = pred,
                                              A0 = A0,
                                              A1 = A1,
                                              x_basis = x_basis,
                                              x_basis_A0 = x_basis_A0,
                                              x_basis_A1 = x_basis_A1),
                                   g = g,
                                   theta = theta,
                                   psi_pound_est = psi_pound_est,
                                   S = S,
                                   A = A,
                                   Y = Y,
                                   n = length(Y),
                                   controls_only = controls_only,
                                   weights = weights)
          PnEIC_prev <- mean(eic)
          prev_pos <- (PnEIC_prev >= 0)
          while (cur_iter <= max_iter && abs(PnEIC_prev) > sn) {
            coefs <- coefs + dx * sign(PnEIC_prev) * direction
            A0 <- as.numeric(x_basis_A0 %*% matrix(coefs))
            A1 <- as.numeric(x_basis_A1 %*% matrix(coefs))
            pred <- as.numeric(x_basis %*% matrix(coefs))
            psi_pound_est <- mean((1-Pi$A0)*A0-(1-Pi$A1)*A1)
            eic <- get_eic_psi_pound(Pi = Pi,
                                     tau = list(pred = pred,
                                                A0 = A0,
                                                A1 = A1,
                                                x_basis = x_basis,
                                                x_basis_A0 = x_basis_A0,
                                                x_basis_A1 = x_basis_A1),
                                     g = g,
                                     theta = theta,
                                     psi_pound_est = psi_pound_est,
                                     S = S,
                                     A = A,
                                     Y = Y,
                                     n = length(Y),
                                     controls_only = controls_only,
                                     weights = weights)
            PnEIC_cur <- mean(eic)
            cur_pos <- (PnEIC_cur >= 0)
            if (cur_pos != prev_pos) dx <- dx / 10
            if (abs(PnEIC_cur - PnEIC_prev) < 1e-6) dx <- dx * 10
            prev_pos  <- cur_pos
            PnEIC_prev <- PnEIC_cur
            sn <- 0.01 * sqrt(var(eic, na.rm = TRUE))/(sqrt(n) * log(n))
            if (verbose) message(sprintf("iter %d: PnEIC = %g", cur_iter, PnEIC_cur))
            cur_iter <- cur_iter + 1
          }
        }
      } else {
        coefs <- mean(pseudo_outcome[delta == 1])
      }
    }
  } else if (method == "HAL") {
    browser()
    # HAL-based R learner
    tau_S_obj <- rHAL(W = as.matrix(cbind(W, A=A)[delta == 1,,drop=FALSE]),
                      A = S[delta == 1],
                      Y = Y[delta == 1],,
                      g1W = Pi$pred[delta == 1],
                      theta = theta,
                      foldid = foldid,
                      weights = weights,
                      enumerate_basis_args = enumerate_basis_args,
                      use_weight = TRUE) # much faster, no need to compute (A-g1W)*phi_W

    # counterfactual predictions
    phi_WA_train <- tau_S_obj$phi_W
    phi_WA <- as.matrix(cbind(1, make_design_matrix(X = as.matrix(cbind(W, A=A)), blist = tau_S_obj$basis_list)))
    phi_W0 <- as.matrix(cbind(1, make_design_matrix(X = as.matrix(cbind(W, A=0)), blist = tau_S_obj$basis_list)))
    phi_W1 <- as.matrix(cbind(1, make_design_matrix(X = as.matrix(cbind(W, A=1)), blist = tau_S_obj$basis_list)))

    # targeting
    if (target_method == "relaxed") {
      # relaxed-fit targeting
      relaxed_fit <- glm(tau_S_obj$pseudo_outcome ~ -1+.,
                         family = "gaussian",
                         data = as.data.frame(phi_WA_train),
                         weights = tau_S_obj$pseudo_weights)
      beta <- as.numeric(coef(relaxed_fit))
      na_idx <- which(is.na(beta))
      if (length(na_idx) > 0) {
        beta <- beta[!is.na(beta)]
        phi_W <- phi_W[, -na_idx, drop = FALSE]
      }
      IM <- t(phi_WA) %*% diag(Pi$pred*(1-Pi$pred)) %*% phi_WA / length(Y)
      IM_inv <- tryCatch({
        solve(IM)
      }, error = function(e) {
        # TODO: add helpful message if verbose
        if (eic_method == "svd_pseudo_inv") {
          svd_pseudo_inv(IM)
        } else if (eic_method == "diag") {
          solve(IM + diag(1e-3, nrow(IM), ncol(IM)))
        } else {
          stop("Unknown eic_method specified.")
        }
      })
    } else if (target_method == "oneshot") {
      # TODO: implement target in a sequence of WMs
      IM <- t(phi_WA) %*% diag(Pi$pred*(1-Pi$pred)) %*% phi_WA / length(Y)
      IM_inv <- tryCatch({
        solve(IM)
      }, error = function(e) {
        # TODO: add helpful message if verbose
        if (eic_method == "svd_pseudo_inv") {
          svd_pseudo_inv(IM)
        } else if (eic_method == "diag") {
          solve(IM + diag(1e-3, nrow(IM), ncol(IM)))
        } else {
          stop("Unknown eic_method specified.")
        }
      })
      beta <- as.vector(tau_A_obj$beta)
      clever_cov <- as.vector(IM_inv %*% colMeans(phi_W))
      H <- (A-g1W)*as.vector(phi_W %*% clever_cov)
      tau <- as.vector(phi_W %*% beta)
      R <- Y-theta-(A-g1W)*tau
      epsilon <- sum(H*R)/sum(H*H)
      beta <- beta+epsilon*clever_cov
    }

    # compute EIC
    cate_WA <- as.vector(phi_WA %*% beta)
    cate_W0 <- as.vector(phi_W0 %*% beta)
    cate_W1 <- as.vector(phi_W1 %*% beta)
    eic <- eic_psi_pound_wm(S = S,
                            Y = Y,
                            A = A,
                            phi_WA = phi_WA,
                            phi_W1 = phi_W1,
                            phi_W0 = phi_W0,
                            g1W = g1W$pred,
                            theta = theta,
                            Pi = Pi,
                            cate_WA = cate_WA,
                            cate_W0 = cate_W0,
                            cate_W1 = cate_W1,
                            weights = weights,
                            controls_only = controls_only,
                            IM_inv = IM_inv)

    eic <- eic_psi_tilde_wm(Y = Y,
                            A = A,
                            phi_W = phi_W,
                            g1W = g1W,
                            theta = theta,
                            cate = cate,
                            weights = weights,
                            eic_method = eic_method,
                            IM_inv = IM_inv)


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
                     parallel = TRUE,
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
      x_basis <- make_counter_design_matrix(
        basis_list = basis_list,
        X_counterfactual = as.matrix(X),
        X_unpenalized = NULL
      )
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
    pseudo_weights = pseudo_weights,
    x_basis_all = as.matrix(cbind(1, X)),
    x_basis_A1_all = as.matrix(X_A1_counter),
    x_basis_A0_all = as.matrix(X_A0_counter),
    eic = eic
  ))
}
