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
#' @param theta_W_tilde A vector of estimated conditional mean of outcome given
#' baseline covariates, \eqn{\theta_W(W)=\mathbb{E}(Y\mid W)}.
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
learn_tau_A <- function(W,
                        A,
                        Y,
                        g1W,
                        delta,
                        theta_W,
                        method,
                        foldid,
                        weights,
                        enumerate_basis_args,
                        pooled_working_model_formula,
                        target_method = "oneshot",
                        eic_method = "svd_pseudo_inv",
                        verbose = FALSE,
                        browse = FALSE) {

  if (browse) browser()

  if (!is.null(pooled_working_model_formula)) {
    # pre-specified working model formula
    pseudo_outcome <- ifelse(abs(A - g) < 1e-10, 0, (Y - theta_W_tilde) / (A - g))
    pseudo_weights <- (A - g)^2 * weights
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
    # main-term lasso-based R learner
    tau_A <- rlasso(W = as.matrix(W[delta == 1,,drop=FALSE]),
                    A = A[delta == 1],
                    Y = Y[delta == 1],,
                    g1W = g1W[delta == 1],
                    theta = theta_W[delta == 1],
                    foldid = foldid[delta == 1],
                    weights = weights[delta == 1],
                    use_weight = TRUE) # much faster, no need to compute (A-g1W)*phi_W
    tau_A$phi_W <- as.matrix(cbind(1, W[, tau_A$non_zero, drop=FALSE]))
    tau_A$cate_W <- as.vector(tau_A$phi_W %*% tau_A$beta)
  } else if (method == "HAL") {
    # HAL-based R learner
    tau_A <- rHAL(W = as.matrix(W[delta == 1,,drop=FALSE]),
                  A = A[delta == 1],
                  Y = Y[delta == 1],,
                  g1W = g1W[delta == 1],
                  theta = theta_W[delta == 1],
                  foldid = foldid[delta == 1],
                  weights = weights[delta == 1],
                  enumerate_basis_args = enumerate_basis_args,
                  use_weight = TRUE) # much faster, no need to compute (A-g1W)*phi_W
    tau_A$phi_W <- as.matrix(cbind(1, make_design_matrix(X = as.matrix(W), blist = tau_A$basis_list)))
    tau_A$cate_W <- as.vector(tau_A$phi_W %*% tau_A$beta)
  }

  return(tau_A)
}
