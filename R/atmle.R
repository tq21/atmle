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
#' @param data A \code{data.frame} containing baseline covariates \eqn{W},
#' binary treatment indicator \eqn{A} (\eqn{A=1} for active treatment),
#' outcome \eqn{Y}, and binary study indicator of whether the observation is
#' from the randomized controlled trial \eqn{S=1} or from the external data
#' \eqn{S=0}. If both studies are observational, then \eqn{S=1} should be the
#' reference study.
#' @param S_node The column index of the \eqn{S} node in \code{data}.
#' @param W_node The column indices of the \eqn{W} node in \code{data}.
#' @param A_node The column index of the \eqn{A} node in \code{data}.
#' @param Y_node The column index of the \eqn{Y} node in \code{data}.
#' @param controls_only A logical indicating whether the external data has
#' only control-arm or both control-arm and treatment-arm.
#' @param family A character string specifying the family of the outcome
#' \eqn{Y}. Currently only \code{"gaussian"} and \code{"binomial"} are
#' supported.
#' @param g_rct The probability of receiving the active treatment in the
#' randomized controlled trial.
#' @param atmle_pooled A logical indicating whether to also use A-TMLE for the
#' pooled-ATE estimand If set to \code{FALSE}, use a regular TMLE for the
#' pooled-ATE, but A-TMLE for the bias-estimand.
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
#' @param g_delta_method The method to estimate the nuisance function
#' \eqn{\delta(W,A)=\mathbb{P}(\Delta=1\mid W,A)}. \code{"glm"} for main-term
#' linear model, \code{"glmnet"} for lasso, \code{"sl3"} for default super
#' learner, or a \code{list} of \code{sl3} learners. Default is \code{"glmnet"}.
#' \eqn{\Delta} is the indicator of missing outcome. \code{1} - observed,
#' \code{0} - missing. Only applicable when there are missing outcomes.
#' @param theta_tilde_method The method to estimate the nuisance function
#' \eqn{\tilde{\theta}(W,A)=\mathbb{E}(Y\mid W,A,S=1)}.
#' \code{"glm"} for main-term linear model, \code{"glmnet"} for lasso,
#' \code{"sl3"} for default super learner, or a \code{list} of \code{sl3}
#' learners. Default is \code{"glmnet"}.
#' @param Q_method The method to estimate the nuisance function
#' \eqn{Q(A,W)=\mathbb{E}(Y\mid W,A)}.
#' \code{"glm"} for main-term linear model, \code{"glmnet"} for lasso,
#' \code{"sl3"} for default super learner, or a \code{list} of \code{sl3}
#' learners. Default is \code{"glmnet"}. Only applicable when \code{atmle_pooled}
#' is set to \code{FALSE}.
#' @param bias_working_model The working model for the bias estimand.
#' Either \code{"glmnet"} for lasso-based working model or \code{"HAL"} for
#' highly adaptive lasso-based working model
#' @param pooled_working_model The working model for the pooled-ATE estimand.
#' Either \code{"glmnet"} for lasso-based working model or \code{"HAL"} for
#' highly adaptive lasso-based working model.
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
#' set.seed(123)
#'
#' n <- 2000
#' S <- rbinom(n, 1, 0.2)
#' W1 <- rnorm(n)
#' W2 <- rnorm(n)
#' W3 <- rnorm(n)
#' A <- numeric(n)
#' g_rct <- 0.67
#' A[S == 1] <- rbinom(sum(S), 1, g_rct)
#' A[S == 0] <- rbinom(n - sum(S), 1, plogis(0.5 * W1[S == 0]))
#' UY <- rnorm(n, 0, 1)
#' Y <- 2.5 + 0.9 * W1 + 1.1 * W2 + 2.7 * W3 + 1.5 * A + UY + (1 - S) * (0.2 + 0.1 * W1 * (1 - A))
#' data <- data.frame(S, W1, W2, A, Y)
#' true_ate <- 1.5
#'
#' res <- atmle(data,
#'   S_node = c(1),
#'   W_node = c(2, 3),
#'   A_node = 4,
#'   Y_node = 5,
#'   controls_only = FALSE,
#'   family = "gaussian",
#'   atmle_pooled = TRUE,
#'   g_rct = g_rct,
#'   verbose = FALSE
#' )
atmle <- function(data,
                  S,
                  W,
                  A,
                  Y,
                  controls_only,
                  family,
                  atmle_pooled = TRUE,
                  theta_method = "glmnet",
                  Pi_method = "glmnet",
                  g_method = "glmnet",
                  g_delta_method = "glmnet",
                  theta_tilde_method = "glmnet",
                  Q_method = "glmnet",
                  bias_working_model = "glmnet",
                  bias_working_model_formula = NULL,
                  pooled_working_model = "glmnet",
                  pooled_working_model_formula = NULL,
                  min_working_model = FALSE,
                  max_degree = 1,
                  v_folds = 5,
                  g_bounds = c(0.01, 0.99),
                  Pi_bounds = c(0.01, 0.99),
                  theta_bounds = c(-Inf, Inf),
                  target_gwt = TRUE,
                  verbose = TRUE,
                  enumerate_basis_args = list(),
                  fit_hal_args = list()) {
  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }

  # extract variables ----------------------------------------------------------
  S <- data[[S]]
  W <- data[, W, drop = FALSE]
  A <- data[[A]]
  Y <- data[[Y]]
  delta <- as.integer(!is.na(Y))
  n <- nrow(data)

  # validate controls_only argument
  if (controls_only & 1 %in% A[S == 0]) {
    stop("The 'controls_only' argument is set to TRUE, but there are treated units in the external data.")
  } else if (!controls_only & sum(A[S == 0]) == 0) {
    stop("The 'controls_only' argument is set to FALSE, but there are only control units in the external data.")
  }

  # cv.glmnet only works when design matrix has at least 2 columns
  # append dummy column of ones if necessary
  if (ncol(W) == 1) {
    W <- cbind(W, 1)
  }

  # cross fitting schemes
  if (family == "gaussian") {
    if (sum(delta) < n) {
      cv_strata <- paste0(S, "-", A, "-", delta)
    } else {
      cv_strata <- paste0(S, "-", A)
    }
  } else if (family == "binomial") {
    if (sum(delta) < n) {
      cv_strata <- paste0(S, "-", A, "-", Y, "-", delta)
    } else {
      cv_strata <- paste0(S, "-", A, "-", Y)
    }
  }

  suppressWarnings({
    folds <- make_folds(
      n = n, V = v_folds,
      strata_ids = as.integer(factor(cv_strata))
    )
  })

  # estimate bias psi_pound ----------------------------------------------------
  # learn nuisance parts
  if (verbose) cat("learning \U03B8(W,A)=E(Y|W,A)...")
  theta <- learn_theta(
    W = W,
    A = A,
    Y = Y,
    delta = delta,
    controls_only = controls_only,
    method = theta_method,
    folds = folds,
    family = family,
    theta_bounds = theta_bounds,
    cross_fit_nuisance = cross_fit_nuisance
  )
  if (verbose) cat("Done!\n")

  if (verbose) cat("learning \U03A0(S=1|W,A)=P(S=1|W,A)...")
  Pi <- learn_Pi(
    S = S,
    W = W,
    A = A,
    controls_only = controls_only,
    method = Pi_method,
    folds = folds,
    Pi_bounds = Pi_bounds
  )
  if (verbose) cat("Done!\n")

  if (verbose) cat("learning g(A=1|W)=P(A=1|W)...")
  g <- learn_g(
    W = W,
    A = A,
    method = g_method,
    folds = folds,
    g_bounds = g_bounds
  )
  if (verbose) cat("Done!\n")

  if (sum(delta) < n) {
    # outcome has missing
    if (verbose) cat("learning g(\U0394=1|S,W,A)=P(\U0394=1|S,W,A)...")
    g_delta <- learn_g_delta(
      S = S,
      W = W,
      A = A,
      delta = delta,
      method = g_delta_method,
      folds = folds,
      g_bounds = g_bounds
    )
    if (verbose) cat("Done!\n")

    if (verbose) cat("learning g(\U0394=1|W,A)=P(\U0394=1|W,A)...")
    g_delta_tilde <- learn_g_delta_tilde(
      W = W,
      A = A,
      delta = delta,
      method = g_delta_method,
      folds = folds,
      g_bounds = g_bounds
    )
    if (verbose) cat("Done!\n")
  } else {
    # no censoring
    g_delta <- g_delta_tilde <- list(
      pred = rep(1, length(A)),
      A0 = rep(1, length(A)),
      A1 = rep(1, length(A))
    )
  }

  # censoring weights
  weights <- delta / g_delta$pred
  weights_tilde <- delta / g_delta_tilde$pred

  # learn working model tau for bias
  if (verbose) cat("learning \U03C4(W,A)=E(Y|S=1,W,A)-E(Y|S=0,W,A)...")
  tau <- learn_tau(
    S = S, W = W, A = A, Y = Y, Pi = Pi, theta = theta, g = g,
    delta = delta,
    controls_only = controls_only,
    method = bias_working_model,
    v_folds = v_folds,
    max_degree = max_degree,
    min_working_model = min_working_model,
    target_gwt = target_gwt,
    Pi_bounds = Pi_bounds,
    enumerate_basis_args = enumerate_basis_args,
    fit_hal_args = fit_hal_args,
    weights = weights,
    bias_working_model_formula = bias_working_model_formula
  )
  if (verbose) cat("Done!\n")

  # TMLE to target Pi
  if (verbose) cat("targeting \U03A0(S=1|W,A)=P(S=1|W,A)...")
  Pi <- Pi_tmle(
    S = S,
    W = W,
    A = A,
    g = g,
    tau = tau,
    Pi = Pi,
    controls_only = controls_only,
    target_gwt = target_gwt,
    Pi_bounds = Pi_bounds
  )
  if (verbose) cat("Done!\n")

  psi_pound_est <- NULL
  if (controls_only) {
    psi_pound_est <- mean((1 - Pi$A0) * tau$A0)
  } else {
    psi_pound_est <- mean((1 - Pi$A0) * tau$A0 - (1 - Pi$A1) * tau$A1)
  }

  # estimate pooled-ATE psi_tilde ----------------------------------------------
  T_working <- NULL
  psi_tilde_est <- NULL
  psi_tilde_eic <- NULL
  if (atmle_pooled) {
    # use atmle for pooled-ATE
    if (verbose) cat("learning \U03B8\U0303(W)=E(Y|W)...")
    theta_tilde <- learn_theta_tilde(
      W = W,
      Y = Y,
      delta = delta,
      method = theta_tilde_method,
      folds = folds,
      family = family,
      theta_bounds = theta_bounds
    )
    if (verbose) cat("Done!\n")

    if (verbose) cat("learning \U03A4(W)=E(Y|W,A=1)-E(Y|W,A=0)...")
    T_working <- learn_T(
      W = W,
      A = A,
      Y = Y,
      g = g,
      delta = delta,
      theta_tilde = theta_tilde,
      method = pooled_working_model,
      min_working_model = min_working_model,
      v_folds = v_folds,
      weights = weights_tilde,
      enumerate_basis_args = enumerate_basis_args,
      fit_hal_args = fit_hal_args,
      pooled_working_model_formula = pooled_working_model_formula
    )
    if (verbose) cat("Done!\n\n")

    # estimates
    psi_tilde_est <- mean(T_working$pred)
    psi_tilde_eic <- get_eic_psi_tilde(T_working, g, theta_tilde, Y, A, n, weights_tilde)
  } else {
    # use regular TMLE for pooled-ATE
    Q <- learn_Q(
      W = W,
      A = A,
      Y = Y,
      delta = delta,
      method = Q_method,
      v_folds = v_folds,
      family = family,
      theta_bounds = theta_bounds
    )
    Q_star <- tmle(
      Y = Y, A = A, W = W, g1W = g,
      Q = as.matrix(data.frame(Q$A1, Q$A0)),
      family = family
    )

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

  # bias parameter
  psi_pound_eic <- get_eic_psi_pound(
    Pi = Pi,
    tau = tau,
    g = g,
    theta = theta,
    psi_pound_est = psi_pound_est,
    S = S,
    A = A,
    Y = Y,
    n = n,
    controls_only = controls_only,
    weights = weights
  )
  psi_pound_se <- sqrt(var(psi_pound_eic, na.rm = TRUE) / n)
  psi_pound_lower <- psi_pound_est - 1.96 * psi_pound_se
  psi_pound_upper <- psi_pound_est + 1.96 * psi_pound_se

  # pooled-ATE parameter
  psi_tilde_se <- sqrt(var(psi_tilde_eic, na.rm = TRUE) / n)
  psi_tilde_lower <- psi_tilde_est - 1.96 * psi_tilde_se
  psi_tilde_upper <- psi_tilde_est + 1.96 * psi_tilde_se

  # RCT-ATE
  est <- psi_tilde_est - psi_pound_est
  eic <- psi_tilde_eic - psi_pound_eic
  se <- sqrt(var(eic, na.rm = TRUE) / n)
  lower <- est - 1.96 * se
  upper <- est + 1.96 * se

  results <- list(
    est = est,
    lower = lower,
    upper = upper,
    psi_pound_est = psi_pound_est,
    psi_pound_lower = psi_pound_lower,
    psi_pound_upper = psi_pound_upper,
    psi_tilde_est = psi_tilde_est,
    psi_tilde_lower = psi_tilde_lower,
    psi_tilde_upper = psi_tilde_upper,
    tau_A1 = tau$A1,
    tau_A0 = tau$A0
  )

  if (verbose) {
    cat("Pooled ATE: ", signif(results$psi_tilde_est, 3), " (", signif(results$psi_tilde_lower, 3), ", ", signif(results$psi_tilde_upper, 3), ")\n", sep = "")
    cat("Bias-corrected ATE: ", signif(results$est, 3), " (", signif(results$lower, 3), ", ", signif(results$upper, 3), ")\n", sep = "")
  }

  return(results)
}
