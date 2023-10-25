#' Adaptive-TMLE
#'
#' @description Adaptive-TMLE (A-TMLE) for estimating the average treatment
#' effect based on combined randomized controlled trial data and real-world
#' data.
#'
#' @export
#'
#' @param data A \code{data.frame} containing baseline covariates \eqn{W}, binary
#' treatment indicator \eqn{A} (\eqn{A=1} for active treatment),
#' outcome \eqn{Y}, and binary indicator of whether the observation is from the
#' randomized controlled trial \eqn{S=1} or from the external data \eqn{S=0}.
#' @param S_node The column indices of the \eqn{S} node in \code{data}.
#' @param W_node The column indices of the \eqn{W} node in \code{data}.
#' @param A_node The column indices of the \eqn{A} node in \code{data}.
#' @param Y_node The column indices of the \eqn{Y} node in \code{data}.
#' @param controls_only A logical indicating whether the external data has
#' only control-arm or both control-arm and treatment-arm.
#' @param family A character string specifying the family of the outcome
#' \eqn{Y}.
#' Currently only \code{"gaussian"} and \code{"binomial"} is supported.
#' @param atmle_pooled A logical indicating whether to use A-TMLE for the
#' pooled-ATE parameter also. If set to \code{FALSE}, use a regular TMLE.
#' Default is \code{TRUE}.
#' @param theta_method The method to estimate the nuisance function
#' \eqn{\theta(W,A)=\mathbb{E}(Y\mid W,A)}.
#' \code{"glm"} for main-term linear model, \code{"glmnet"} for lasso,
#' or a \code{list} of \code{sl3} learners for super learner-based estimation.
#' Default is \code{"glmnet"}.
#' @param Pi_method The method to estimate the nuisance function
#' \eqn{\Pi(S\mid W,A)=\mathbb{P}(S=1\mid W,A)}.
#' \code{"glm"} for main-term linear model, \code{"glmnet"} for lasso,
#' or a \code{list} of \code{sl3} learners for super learner-based estimation.
#' @param g_method The method to estimate the nuisance function
#' \eqn{g(A\mid W)=\mathbb{P}(A=1\mid W)}.
#' \code{"glm"} for main-term linear model, \code{"glmnet"} for lasso,
#' or a \code{list} of \code{sl3} learners for super learner-based estimation.
#' @param theta_tilde_method The method to estimate the nuisance function
#' \eqn{\tilde{\theta}(W,A)=\mathbb{E}(Y\mid W,A,S=1)}.
#' \code{"glm"} for main-term linear model, \code{"glmnet"} for lasso,
#' or a \code{list} of \code{sl3} learners for super learner-based estimation.
#' @param Q_method The method to estimate the nuisance function
#' \eqn{Q(A,W)=\mathbb{E}(Y\mid W,A)}.
#' \code{"glm"} for main-term linear model, \code{"glmnet"} for lasso,
#' or a list of \code{sl3} learners for super learner-based estimation.
#' @param bias_working_model The working model for the bias estimand.
#' Either \code{"glmnet"} for lasso-based working model or \code{"HAL"} for highly adaptive
#' lasso-based working model. Default is \code{"glmnet"}.
#' @param pooled_working_model The working model for the pooled-ATE estimand.
#' Either \code{"glmnet"} for lasso-based working model or \code{"HAL"} for highly adaptive
#' lasso-based working model. Default is \code{"glmnet"}.
#' @param g_rct The probability of receiving the active treatment in the
#' randomized controlled trial.
#' @param var_method The method to estimate the variance of the A-TMLE
#' estimator. Either \code{"ic"} for influence curve-based variance estimator or
#' \code{"bootstrap"} for bootstrap-based variance estimator. Default is \code{"ic"}.
#' @param max_iter The maximum number of iterations for TMLE targeting of
#' \eqn{\Pi(S\mid W,A)=\mathbb{P}(S=1\mid W,A)}. Default is \code{1}.
#' @param v_folds The number of folds for cross-validation (whenever necessary).
#' Default is \code{5}.
#' @param verbose A logical indicating whether to print out the progress.
#' Default is \code{TRUE}.
#' @returns A \code{list} containing the following elements:
#' \item{ate}{The estimated average treatment effect;}
#' \item{ate_se}{The estimated standard error of the average treatment effect;}
#' \item{ate_ci}{The estimated \eqn{95\%} confidence interval for the average treatment
#' effect.}
atmle <- function(data,
                  S_node,
                  W_node,
                  A_node,
                  Y_node,
                  controls_only,
                  family = "gaussian",
                  atmle_pooled = TRUE,
                  theta_method = "glm",
                  Pi_method = "glm",
                  g_method = "glm",
                  theta_tilde_method = "glm",
                  Q_method = "glm",
                  bias_working_model = "glmnet",
                  pooled_working_model = "glmnet",
                  g_rct,
                  var_method = "ic",
                  max_iter = 1,
                  v_folds = 5,
                  verbose = TRUE) {

  # sanity checks --------------------------------------------------------------
  check_data_and_nodes(data, S_node, W_node, A_node, Y_node)
  check_learners(theta_method, Pi_method, g_method, theta_tilde_method,
                 Q_method, bias_working_model, pooled_working_model)
  check_args(controls_only, atmle_pooled, var_method, g_rct, verbose)

  # define nodes ---------------------------------------------------------------
  S <- data[, S_node] # study indicator
  W <- data[, W_node] # covariates
  A <- data[, A_node] # treatment
  Y <- data[, Y_node] # outcome
  n <- nrow(data) # sample size

  # estimate bias psi_pound ----------------------------------------------------
  # learn nuisance parts
  if (verbose) print("learning \U03B8(W,A)=E(Y|W,A)")
  theta <- learn_theta(W, A, Y, controls_only, theta_method, v_folds)

  if (verbose) print("learning \U03A0(S=1|W,A)=P(S=1|W,A)")
  Pi <- learn_Pi(S, W, A, controls_only, Pi_method, v_folds)

  if (verbose) print("learning g(A=1|W)=P(A=1|W)")
  g <- learn_g(S, W, A, g_rct, controls_only, g_method, v_folds)

  # learn working model tau for bias
  if (verbose) print("learning \U03C4(Y|S,W,A)=E(Y|S,W,A)")
  tau <- learn_tau(S, W, A, Y, Pi, theta, controls_only, bias_working_model,
                   v_folds)

  # TMLE to target Pi
  if (verbose) print("targeting \U03A0(S=1|W,A)=P(S=1|W,A)")
  iter <- 0
  while (iter < max_iter) {
    # TODO: check empirical mean of IC, iterate until convergence
    # targeted Pi
    Pi_star <- Pi_tmle(S, W, A, g, tau, Pi, controls_only)

    # re-learn working model tau with targeted Pi
    tau_star <- learn_tau(S, W, A, Y, Pi_star, theta, controls_only,
                          bias_working_model, v_folds)

    Pi <- Pi_star
    tau <- tau_star
    iter <- iter + 1
  }

  psi_pound_est <- NULL
  if (controls_only) {
    psi_pound_est <- mean((1-Pi$A0)*tau$A0)
  } else {
    psi_pound_est <- mean((1-Pi$A0)*tau$A0-(1-Pi$A1)*tau$A1)
  }

  # estimate pooled-ATE psi_tilde ----------------------------------------------
  T_working <- NULL
  psi_tilde_est <- NULL
  psi_tilde_eic <- NULL
  if (atmle_pooled) {
    # use atmle for pooled-ATE
    if (verbose) print("learning \U03B8\U0303(W)=E(Y|W)")
    theta_tilde <- learn_theta_tilde(W, Y, theta_tilde_method, v_folds, family)

    if (verbose) print("learning \U03A4(W)=E(Y|W,A=1)-E(Y|W,A=0)")
    T_working <- learn_T(W, A, Y, g, theta_tilde, pooled_working_model, v_folds)

    # estimates
    psi_tilde_est <- mean(T_working$pred)
    psi_tilde_eic <- get_eic_psi_tilde(T_working, g, theta_tilde, Y, A, n)
  } else {
    # use regular TMLE for pooled-ATE
    Q <- learn_Q(W, A, Y, method = Q_method)
    Q_star <- tmle(Y = Y, A = A, W = W, g1W = g,
                   Q = as.matrix(data.frame(Q$A1, Q$A0)),
                   family = "gaussian")

    # estimates
    psi_tilde_est <- Q_star$estimates$ATE$psi
    psi_tilde_eic <- Q_star$estimates$IC$IC.ATE
  }

  # final estimates ------------------------------------------------------------
  est <- NULL
  lower <- NULL
  upper <- NULL
  psi_pound_lower <- NULL
  psi_pound_upper <- NULL
  psi_tilde_lower <- NULL
  psi_tilde_upper <- NULL

  if (var_method == "ic") {
    # bias parameter
    psi_pound_eic <- get_eic_psi_pound(Pi, tau, g, theta, psi_pound_est, S, A, Y, n, controls_only)
    psi_pound_se <- sqrt(var(psi_pound_eic, na.rm = TRUE)/n)
    psi_pound_lower <- psi_pound_est+qnorm(0.025)*psi_pound_se
    psi_pound_upper <- psi_pound_est+qnorm(0.975)*psi_pound_se

    # pooled-ATE parameter
    psi_tilde_se <- sqrt(var(psi_tilde_eic, na.rm = TRUE)/n)
    psi_tilde_lower <- psi_tilde_est+qnorm(0.025)*psi_tilde_se
    psi_tilde_upper <- psi_tilde_est+qnorm(0.975)*psi_tilde_se

    # RCT-ATE
    est <- psi_tilde_est - psi_pound_est
    eic <- psi_tilde_eic - psi_pound_eic
    se <- sqrt(var(eic, na.rm = TRUE)/n)
    lower <- est+qnorm(0.025)*se
    upper <- est+qnorm(0.975)*se
  } else if (var_method == "bootstrap") {
    # bias parameter
    psi_pound_se <- bootstrap_psi_pound(tau, W, Pi)
    psi_pound_lower <- psi_pound_est+qnorm(0.025)*psi_pound_se
    psi_pound_upper <- psi_pound_est+qnorm(0.975)*psi_pound_se

    # pooled-ATE parameter (use ic-based for now)
    psi_tilde_se <- NULL
    #if (atmle_pooled) {
    #  psi_tilde_se <- bootstrap_psi_tilde(W, psi_tilde)
    #  psi_tilde_lower <- psi_tilde_est+qnorm(0.025)*psi_tilde_se
    #  psi_tilde_upper <- psi_tilde_est+qnorm(0.975)*psi_tilde_se
    #} else {
    psi_tilde_se <- sqrt(var(psi_tilde_eic, na.rm = TRUE)/n)
    psi_tilde_lower <- psi_tilde_est+qnorm(0.025)*psi_tilde_se
    psi_tilde_upper <- psi_tilde_est+qnorm(0.975)*psi_tilde_se
    #}

    # RCT-ATE
    est <- psi_tilde_est - psi_pound_est
    se <- sqrt(psi_pound_se^2+psi_tilde_se^2)
    lower <- est+qnorm(0.025)*se
    upper <- est+qnorm(0.975)*se
  }

  return(list(est = est,
              lower = lower,
              upper = upper,
              psi_pound_est = psi_pound_est,
              psi_pound_lower = psi_pound_lower,
              psi_pound_upper = psi_pound_upper,
              psi_tilde_est = psi_tilde_est,
              psi_tilde_lower = psi_tilde_lower,
              psi_tilde_upper = psi_tilde_upper))
}
