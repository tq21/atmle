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
                        theta,
                        g1W,
                        delta,
                        v_folds,
                        weights,
                        enumerate_basis_args,
                        dx,
                        max_iter) {

  # R-transformations
  pseudo_outcome <- ifelse(abs(A - g1W) < 1e-10, 0, (Y - theta) / (A - g1W))
  pseudo_weights <- (A - g1W)^2 * weights

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

  X <- data.frame(W)

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
  phi_W <- cbind(1, X_hal_selected)
  beta <- as.numeric(coef(fit, s = "lambda.min")[coef(fit, s = "lambda.min") != 0])

  if (length(non_zero) > 0) {
    # regularized targeting ----------------------------------------------------
    cur_iter <- 1
    PnEIC <- Inf
    sn <- 0
    beta_reg <- beta
    while (cur_iter <= max_iter & abs(PnEIC) > sn) {
      # compute canonical gradient
      cate_pred <- as.numeric(phi_W%*%beta_reg)
      eic_np <- eic_ate(QW1 = theta+(1-g1W)*cate_pred,
                        QW0 = theta-g1W*cate_pred,
                        psi = mean(cate_pred),
                        A = A,
                        g1W = g1W,
                        Y = Y,
                        QWA = theta+(A-g1W)*cate_pred)
      PnEIC <- mean(eic_np)

      # obtain current score matrix
      score_mat <- (Y-theta-(A-g1W)*cate_pred)*(A-g1W)*phi_W

      # project gradient onto the current score space (regularized)
      proj_fit <- glmnet(x = score_mat,
                         y = eic_np,
                         intercept = FALSE,
                         lambda = 1e-5,
                         alpha = 1)
      direction <- as.numeric(coef(proj_fit, s = 1e-5)[-1])

      # update beta
      beta_reg <- beta_reg + dx*sign(PnEIC)*direction
      sn <- sqrt(var(eic_np, na.rm = TRUE))/(sqrt(length(A)) * log(length(A)))
      cur_iter <- cur_iter + 1
      #print(PnEIC)
    }

    # relaxed HAL ------------------------------------------------------------
    fit <- glm(pseudo_outcome[delta == 1] ~ .,
               family = "gaussian",
               data = data.frame(as.matrix(X_hal_selected)),
               weights = pseudo_weights[delta == 1])
    beta_relax <- as.numeric(coef(fit))
    beta_relax[is.na(beta_relax)] <- 0

  } else {
    beta_reg <- beta_relax <- mean(pseudo_outcome[delta == 1])
  }

  x_basis <- make_counter_design_matrix(basis_list = basis_list,
                                        X_counterfactual = as.matrix(X),
                                        X_unpenalized = NULL)
  pred_reg <- as.numeric(x_basis%*%matrix(beta_reg))
  pred_relax <- as.numeric(x_basis%*%matrix(beta_relax))

  return(list(reg = list(pred = pred_reg,
                         x_basis = x_basis,
                         coefs = beta_reg,
                         non_zero = non_zero),
              relax = list(pred = pred_relax,
                           x_basis = x_basis,
                           coefs = beta_relax,
                           non_zero = non_zero)))
}
