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
#' @param theta_WA A vector of estimated conditional mean of outcome given
#' baseline covariates and treatment, \eqn{\theta_WA(W,A)=\mathbb{E}(Y\mid W,A)}.
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
                        theta_WA,
                        g1W,
                        delta,
                        controls_only,
                        method,
                        v_folds,
                        foldid,
                        min_working_model,
                        target_gwt,
                        Pi_bounds,
                        enumerate_basis_args,
                        weights,
                        bias_working_model_formula,
                        verbose = TRUE) {

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
    # main-term lasso-based R learner
    tau_S <- rlasso(W = as.matrix(cbind(W, A=A, W*A)[delta == 1,,drop=FALSE]),
                    A = S[delta == 1],
                    Y = Y[delta == 1],,
                    g1W = Pi$pred[delta == 1],
                    theta = theta_WA[delta == 1],
                    foldid = foldid[delta == 1],
                    weights = weights[delta == 1],
                    use_weight = TRUE) # much faster, no need to compute (S-Pi)*phi_WA

    # counterfactual predictions
    tau_S$phi_WA <- as.matrix(cbind(1, cbind(W, A=A, W*A)[,tau_S$non_zero,drop=FALSE]))
    tau_S$phi_W0 <- as.matrix(cbind(1, cbind(W, A=0, W*0)[,tau_S$non_zero,drop=FALSE]))
    tau_S$phi_W1 <- as.matrix(cbind(1, cbind(W, A=1, W*1)[,tau_S$non_zero,drop=FALSE]))
    tau_S$cate_WA <- as.vector(tau_S$phi_WA %*% tau_S$beta)
    tau_S$cate_W0 <- as.vector(tau_S$phi_W0 %*% tau_S$beta)
    tau_S$cate_W1 <- as.vector(tau_S$phi_W1 %*% tau_S$beta)
  } else if (method == "HAL") {
    # HAL-based R learner
    tau_S <- rHAL(W = as.matrix(cbind(W, A=A)[delta == 1,,drop=FALSE]),
                  A = S[delta == 1],
                  Y = Y[delta == 1],,
                  g1W = Pi$pred[delta == 1],
                  theta = theta_WA[delta == 1],
                  foldid = foldid[delta == 1],
                  weights = weights[delta == 1],
                  enumerate_basis_args = enumerate_basis_args,
                  use_weight = TRUE) # much faster, no need to compute (S-Pi)*phi_WA

    # counterfactual predictions
    tau_S$phi_WA <- as.matrix(cbind(1, make_design_matrix(X = as.matrix(cbind(W, A=A)), blist = tau_S$basis_list)))
    tau_S$phi_W0 <- as.matrix(cbind(1, make_design_matrix(X = as.matrix(cbind(W, A=0)), blist = tau_S$basis_list)))
    tau_S$phi_W1 <- as.matrix(cbind(1, make_design_matrix(X = as.matrix(cbind(W, A=1)), blist = tau_S$basis_list)))
    tau_S$cate_WA <- as.vector(tau_S$phi_WA %*% tau_S$beta)
    tau_S$cate_W0 <- as.vector(tau_S$phi_W0 %*% tau_S$beta)
    tau_S$cate_W1 <- as.vector(tau_S$phi_W1 %*% tau_S$beta)
  } else {
    stop("bias_working_model must be either 'glmnet' or 'HAL'")
  }

  return(tau_S)
}
