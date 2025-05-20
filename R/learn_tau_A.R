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
learn_tau_A <- function(W,
                        A,
                        Y,
                        g1W,
                        delta,
                        theta,
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
    pseudo_outcome <- ifelse(abs(A - g) < 1e-10, 0, (Y - theta_tilde) / (A - g))
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
    fit <- cv.glmnet(x = as.matrix(W[delta == 1, ]),
                     y = pseudo_outcome[delta == 1],
                     weights = pseudo_weights[delta == 1],
                     family = "gaussian", keep = TRUE, nfolds = v_folds,
                     alpha = 1)
    non_zero <- which(as.numeric(coef(fit, s = "lambda.min")) != 0)
    coefs <- coef(fit, s = "lambda.min")[non_zero]
    x_basis <- as.matrix(cbind(1, W)[, non_zero, drop=FALSE])

    pred <- drop(x_basis %*% coefs)
  } else if (method == "HAL") {
    # HAL-based R learner
    tau_A_obj <- rHAL(W = as.matrix(W[delta == 1,,drop=FALSE]),
                      A = A[delta == 1],
                      Y = Y[delta == 1],,
                      g1W = g1W,
                      theta = theta,
                      foldid = foldid,
                      weights = weights,
                      enumerate_basis_args = enumerate_basis_args,
                      use_weight = TRUE) # much faster, no need to compute (A-g1W)*phi_W

    # targeting
    # TODO: need to predict on all in case there is missing Y
    phi_W_train <- tau_A_obj$phi_W
    phi_W <- as.matrix(cbind(1, make_design_matrix(X = as.matrix(W), blist = tau_A_obj$basis_list)))
    IM <- t(phi_W) %*% diag(g1W*(1-g1W)) %*% phi_W / length(Y)
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
    if (target_method == "relaxed") {
      # relaxed-fit targeting
      relaxed_fit <- glm(tau_A_obj$pseudo_outcome ~ -1+.,
                         family = "gaussian",
                         data = as.data.frame(phi_W_train),
                         weights = tau_A_obj$pseudo_weights)
      beta <- as.numeric(coef(relaxed_fit))
      na_idx <- which(is.na(beta))
      if (length(na_idx) > 0) {
        beta <- beta[!is.na(beta)]
        phi_W <- phi_W[, -na_idx, drop = FALSE]
        IM <- t(phi_W) %*% diag(g1W*(1-g1W)) %*% phi_W / length(Y)
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
      }
    } else if (target_method == "oneshot") {
      # TODO: implement target in a sequence of WMs
      beta <- as.vector(tau_A_obj$beta)
      clever_cov <- as.vector(IM_inv %*% colMeans(phi_W))
      H <- (A-g1W)*as.vector(phi_W %*% clever_cov)
      tau <- as.vector(phi_W %*% beta)
      R <- Y-theta-(A-g1W)*tau
      epsilon <- sum(H*R)/sum(H*H)
      beta <- beta+epsilon*clever_cov
    }

    # compute EIC
    cate <- as.vector(phi_W %*% beta)
    eic <- eic_psi_tilde_wm(Y = Y,
                            A = A,
                            phi_W = phi_W,
                            g1W = g1W,
                            theta = theta,
                            cate = cate,
                            weights = weights,
                            eic_method = eic_method,
                            IM_inv = IM_inv)
  }

  return(list(cate = cate,
              eic = eic,
              phi_W = phi_W))
}
