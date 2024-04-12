#' @title Adaptive-TMLE for RCT-ATE
#'
#' @description Adaptive-TMLE (A-TMLE) for estimating the average treatment
#' effect based on combined randomized controlled trial data and real-world
#' data.
#'
#' @export
#'
#' @importFrom tmle tmle
#' @importFrom origami make_folds
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
#' \code{"sl3"} for default super learner, or a \code{list} of \code{sl3}
#' learners. Default is \code{"glmnet"}.
#' @param Pi_method The method to estimate the nuisance function
#' \eqn{\Pi(S\mid W,A)=\mathbb{P}(S=1\mid W,A)}.
#' \code{"glm"} for main-term linear model, \code{"glmnet"} for lasso,
#' \code{"sl3"} for default super learner, or a \code{list} of \code{sl3}
#' learners. Default is \code{"glmnet"}.
#' @param g_method The method to estimate the nuisance function
#' \eqn{g(A\mid W)=\mathbb{P}(A=1\mid W)}.
#' \code{"glm"} for main-term linear model, \code{"glmnet"} for lasso,
#' \code{"sl3"} for default super learner, or a \code{list} of \code{sl3}
#' learners. Default is \code{"glmnet"}.
#' @param theta_tilde_method The method to estimate the nuisance function
#' \eqn{\tilde{\theta}(W,A)=\mathbb{E}(Y\mid W,A,S=1)}.
#' \code{"glm"} for main-term linear model, \code{"glmnet"} for lasso,
#' \code{"sl3"} for default super learner, or a \code{list} of \code{sl3}
#' learners. Default is \code{"glmnet"}.
#' @param Q_method The method to estimate the nuisance function
#' \eqn{Q(A,W)=\mathbb{E}(Y\mid W,A)}.
#' \code{"glm"} for main-term linear model, \code{"glmnet"} for lasso,
#' \code{"sl3"} for default super learner, or a \code{list} of \code{sl3}
#' learners. Default is \code{"glmnet"}.
#' @param bias_working_model The working model for the bias estimand.
#' Either \code{"glmnet"} for lasso-based working model or \code{"HAL"} for highly adaptive
#' lasso-based working model
#' (PLEASE USE "glmnet" FOR NOW). Default is \code{"glmnet"}.
#' @param pooled_working_model The working model for the pooled-ATE estimand.
#' Either \code{"glmnet"} for lasso-based working model or \code{"HAL"} for highly adaptive
#' lasso-based working model.
#' (PLEASE USE "glmnet" FOR NOW). Default is \code{"glmnet"}.
#' @param g_rct The probability of receiving the active treatment in the
#' randomized controlled trial.
#' @param var_method The method to estimate the variance of the A-TMLE
#' estimator. Either \code{"ic"} for influence curve-based variance estimator or
#' \code{"bootstrap"} for bootstrap-based variance estimator
#' (PLEASE USE INFLUENCE CURVE-BASED FOR NOW). Default is \code{"ic"}.
#' @param max_iter The maximum number of iterations for TMLE targeting of
#' \eqn{\Pi(S\mid W,A)=\mathbb{P}(S=1\mid W,A)}
#' (PLEASE SET TO 1 FOR NOW). Default is \code{1}.
#' @param v_folds The number of folds for cross-validation (whenever necessary).
#' Default is \code{5}.
#' @param g_bounds A numeric vector of lower and upper bounds for the
#' treatment mechanism. The first element is the lower bound, and the second
#' element is the upper bound. Default is \code{c(0.01, 0.99)}.
#' @param Pi_bounds A numeric vector of lower and upper bounds for the
#' trial enrollment probabilities. The first element is the lower bound,
#' and the second element is the upper bound. Default is \code{c(0.01, 0.99)}.
#' @param theta_bounds A numeric vector of lower and upper bounds for the
#' conditional mean of outcome given baseline covariates and treatment.
#' The first element is the lower bound, and the second element is the upper
#' bound. Default is \code{c(-Inf, Inf)}.
#' @param target_gwt If \code{TRUE}, the treatment mechanism is moved from the
#' denominator of the clever covariate to the weight when fitting the TMLE
#' submodel.
#' @param verbose A logical indicating whether to print out the progress.
#' Default is \code{TRUE}.
#' @returns A \code{list} containing the following elements:
#' \item{ate}{The estimated average treatment effect;}
#' \item{lower}{The lower bound of the \eqn{95\%} confidence interval for the
#' average treatment effect;}
#' \item{upper}{The upper bound of the \eqn{95\%} confidence interval for the
#' average treatment effect;}
#'
#' @examples
#' # simulate data
#' set.seed(123)
#' n <- 2000
#' S <- rbinom(n, 1, 0.5)
#' W1 <- rnorm(n); W2 <- rnorm(n); W <- cbind(W1, W2)
#' A <- numeric(n)
#' A[S == 1] <- rbinom(sum(S), 1, 0.67)
#' A[S == 0] <- rbinom(n-sum(S), 1, plogis(1.2*W1-0.9*W2))
#' UY <- rnorm(n, 0, 1)
#' U_bias <- rnorm(n, 0, 0.5)
#' Y <- -0.5-0.8*W1-1.1*W2+1.5*A+UY+(1-S)*(0.9+2.6*W1)+(1-S)*U_bias
#' data <- data.frame(S, W1, W2, A, Y)
#'
#' atmle_res <- atmle(data,
#'                    S_node = c(1),
#'                    W_node = c(2, 3),
#'                    A_node = 4,
#'                    Y_node = 5,
#'                    controls_only = FALSE,
#'                    family = "gaussian",
#'                    atmle_pooled = TRUE,
#'                    g_rct = 0.67,
#'                    verbose = FALSE)
atmle <- function(data,
                  S_node,
                  W_node,
                  A_node,
                  Y_node,
                  controls_only,
                  family,
                  atmle_pooled = TRUE,
                  theta_method = "glmnet",
                  Pi_method = "glmnet",
                  g_method = "glmnet",
                  theta_tilde_method = "glmnet",
                  Q_method = "glmnet",
                  bias_working_model = "glmnet",
                  pooled_working_model = "glmnet",
                  min_working_model = FALSE,
                  min_working_model_screen = FALSE,
                  undersmooth = 0,
                  g_rct,
                  cross_fit_nuisance = TRUE,
                  enumerate_basis_args = list(),
                  fit_hal_args = list(),
                  var_method = "ic",
                  max_degree = 1,
                  max_iter = 1,
                  v_folds = 5,
                  g_bounds = c(0.01, 0.99),
                  Pi_bounds = c(0.01, 0.99),
                  theta_bounds = c(-Inf, Inf),
                  target_gwt = TRUE,
                  verbose = TRUE) {

  # define nodes ---------------------------------------------------------------
  S <- data[, S_node] # study indicator
  W <- data[, W_node] # covariates
  A <- data[, A_node] # treatment
  Y <- data[, Y_node] # outcome
  n <- nrow(data) # sample size

  # cross fitting schemes
  if (family == "gaussian") {
    cv_strata <- paste0(S, "-", A)
  } else if (family == "binomial") {
    cv_strata <- paste0(S, "-", A, "-", Y)
  }

  suppressWarnings({
    folds <- make_folds(n = n, V = v_folds,
                        strata_ids = as.integer(factor(cv_strata)))
  })

  # estimate bias psi_pound ----------------------------------------------------
  # learn nuisance parts
  if (verbose) print("learning \U03B8(W,A)=E(Y|W,A)")
  theta <- learn_theta(W = W,
                       A = A,
                       Y = Y,
                       controls_only = controls_only,
                       method = theta_method,
                       folds = folds,
                       family = family,
                       theta_bounds = theta_bounds,
                       cross_fit_nuisance = cross_fit_nuisance)

  if (verbose) print("learning \U03A0(S=1|W,A)=P(S=1|W,A)")
  Pi <- learn_Pi(S = S,
                 W = W,
                 A = A,
                 controls_only = controls_only,
                 method = Pi_method,
                 folds = folds,
                 Pi_bounds = Pi_bounds,
                 cross_fit_nuisance = cross_fit_nuisance)

  if (verbose) print("learning g(A=1|W)=P(A=1|W)")
  g <- learn_g(S = S,
               W = W,
               A = A,
               g_rct = g_rct,
               controls_only = controls_only,
               method = g_method,
               folds = folds,
               g_bounds = g_bounds,
               cross_fit_nuisance = cross_fit_nuisance)

  # learn working model tau for bias
  if (verbose) print("learning \U03C4(Y|S,W,A)=E(Y|S,W,A)")
  tau <- learn_tau(S = S, W = W, A = A, Y = Y, Pi = Pi, theta = theta, g = g,
                   controls_only = controls_only,
                   method = bias_working_model,
                   v_folds = v_folds,
                   max_degree = max_degree,
                   min_working_model = min_working_model,
                   min_working_model_screen = min_working_model_screen,
                   undersmooth = undersmooth,
                   target_gwt = target_gwt,
                   Pi_bounds = Pi_bounds,
                   enumerate_basis_args = enumerate_basis_args,
                   fit_hal_args = fit_hal_args)

  if (verbose) print("targeting \U03A0(S=1|W,A)=P(S=1|W,A)")
  if (undersmooth == 3) {
    # use tau_target for targeting Pi
    Pi <- Pi_tmle(S, W, A, g, tau$tau_target, Pi, controls_only, target_gwt, Pi_bounds)
    tau <- tau$tau
  } else if (undersmooth == 4) {
    # use tau_target for targeting Pi
    Pi <- Pi_tmle(S, W, A, g, tau$tau_target, Pi, controls_only, target_gwt, Pi_bounds)
  } else {
    Pi <- Pi_tmle(S, W, A, g, tau, Pi, controls_only, target_gwt, Pi_bounds)
  }

  ## LOG
  tmp_log <- list()
  tmp_log$epsilon <- Pi$epsilon

  psi_pound_est <- NULL
  if (controls_only) {
    if (undersmooth == 4) {
      psi_pound_est <- mean((1-Pi$A0)*tau$tau$A0)
    } else {
      psi_pound_est <- mean((1-Pi$A0)*tau$A0)
    }

  } else {
    if (undersmooth == 4) {
      psi_pound_est <- mean((1-Pi$A0)*tau$tau$A0-(1-Pi$A1)*tau$tau$A1)
    } else {
      psi_pound_est <- mean((1-Pi$A0)*tau$A0-(1-Pi$A1)*tau$A1)
    }
  }

  # estimate pooled-ATE psi_tilde ----------------------------------------------
  T_working <- NULL
  psi_tilde_est <- NULL
  psi_tilde_eic <- NULL
  if (atmle_pooled) {
    # use atmle for pooled-ATE
    if (verbose) print("learning \U03B8\U0303(W)=E(Y|W)")
    theta_tilde <- learn_theta_tilde(W = W,
                                     Y = Y,
                                     method = theta_tilde_method,
                                     folds = folds,
                                     family = family,
                                     theta_bounds = theta_bounds,
                                     cross_fit_nuisance = cross_fit_nuisance)

    if (verbose) print("learning \U03A4(W)=E(Y|W,A=1)-E(Y|W,A=0)")
    T_working <- learn_T(W = W,
                         A = A,
                         Y = Y,
                         g = g,
                         theta_tilde = theta_tilde,
                         method = pooled_working_model,
                         min_working_model = min_working_model,
                         v_folds = v_folds)

    # estimates
    psi_tilde_est <- mean(T_working$pred)
    psi_tilde_eic <- get_eic_psi_tilde(T_working, g, theta_tilde, Y, A, n)
  } else {
    # use regular TMLE for pooled-ATE
    Q <- learn_Q(W, A, Y, method = Q_method)
    Q_star <- tmle(Y = Y, A = A, W = W, g1W = g,
                   Q = as.matrix(data.frame(Q$A1, Q$A0)),
                   family = family)

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

  tryCatch({
    if (var_method == "ic") {
      if (undersmooth == 4) {
        psi_pound_eic <- get_eic_psi_pound(Pi, tau$tau_target, g, theta, psi_pound_est, S, A, Y, n, controls_only)
      } else {
        psi_pound_eic <- get_eic_psi_pound(Pi, tau, g, theta, psi_pound_est, S, A, Y, n, controls_only)
      }

      # psi pound
      psi_pound_se <- sqrt(var(psi_pound_eic, na.rm = TRUE)/n)
      psi_pound_lower <- psi_pound_est - 1.96*psi_pound_se
      psi_pound_upper <- psi_pound_est + 1.96*psi_pound_se

      # psi tilde
      psi_tilde_se <- sqrt(var(psi_tilde_eic, na.rm = TRUE)/n)
      psi_tilde_lower <- psi_tilde_est - 1.96*psi_tilde_se
      psi_tilde_upper <- psi_tilde_est + 1.96*psi_tilde_se

      # final estimate
      est <- psi_tilde_est - psi_pound_est
      eic <- psi_tilde_eic - psi_pound_eic
      se <- sqrt(var(eic, na.rm = TRUE)/n)
      lower <- est-1.96*se
      upper <- est+1.96*se
    } else if (var_method == "bootstrap") {
      est <- psi_tilde_est - psi_pound_est
      if (undersmooth == 4) {
        ci <- bootstrap_ci(tau = tau$tau_target,
                           T_working = T_working,
                           Pi = Pi)
      } else {
        ci <- bootstrap_ci(tau = tau,
                           T_working = T_working,
                           Pi = Pi)
      }
      lower <- ci$lower
      upper <- ci$upper
    }

    if (undersmooth == 4) {
      results <- list(est = est,
                      lower = lower,
                      upper = upper,
                      psi_pound_est = psi_pound_est,
                      psi_pound_lower = psi_pound_lower,
                      psi_pound_upper = psi_pound_upper,
                      psi_tilde_est = psi_tilde_est,
                      psi_tilde_lower = psi_tilde_lower,
                      psi_tilde_upper = psi_tilde_upper,
                      tau_A1 = tau$tau$A1,
                      tau_A0 = tau$tau$A0)
    } else {
      results <- list(est = est,
                      lower = lower,
                      upper = upper,
                      psi_pound_est = psi_pound_est,
                      psi_pound_lower = psi_pound_lower,
                      psi_pound_upper = psi_pound_upper,
                      psi_tilde_est = psi_tilde_est,
                      psi_tilde_lower = psi_tilde_lower,
                      psi_tilde_upper = psi_tilde_upper,
                      tau_A1 = tau$A1,
                      tau_A0 = tau$A0)
    }

    return(c(results, tmp_log))
  }, error = function(e) {
    print("non-invertible information matrix, using bootstrap inference...")
    est <- psi_tilde_est - psi_pound_est
    if (undersmooth == 4) {
      ci <- bootstrap_ci(tau = tau$tau_target,
                         T_working = T_working,
                         Pi = Pi)
    } else {
      ci <- bootstrap_ci(tau = tau,
                         T_working = T_working,
                         Pi = Pi)
    }
    lower <- ci$lower
    upper <- ci$upper

    if (undersmooth == 4) {
      results <- list(est = est,
                      lower = lower,
                      upper = upper,
                      psi_pound_est = psi_pound_est,
                      psi_pound_lower = psi_pound_lower,
                      psi_pound_upper = psi_pound_upper,
                      psi_tilde_est = psi_tilde_est,
                      psi_tilde_lower = psi_tilde_lower,
                      psi_tilde_upper = psi_tilde_upper,
                      tau_A1 = tau$tau$A1,
                      tau_A0 = tau$tau$A0)
    } else {
      results <- list(est = est,
                      lower = lower,
                      upper = upper,
                      psi_pound_est = psi_pound_est,
                      psi_pound_lower = psi_pound_lower,
                      psi_pound_upper = psi_pound_upper,
                      psi_tilde_est = psi_tilde_est,
                      psi_tilde_lower = psi_tilde_lower,
                      psi_tilde_upper = psi_tilde_upper,
                      tau_A1 = tau$A1,
                      tau_A0 = tau$A0)
    }

    return(c(results, tmp_log))
  })
}
