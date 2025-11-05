#' @title Adaptive-TMLE
#'
#' @description Adaptive-TMLE for estimating the average treatment effect based
#' on randomized controlled trial augmented with real-world data.
#'
#' @export
#'
#' @importFrom purrr map
#' @importFrom origami make_folds folds2foldvec fold_from_foldvec
#'
#' @param data A \code{data.frame} containing baseline covariates \eqn{W},
#' binary treatment indicator \eqn{A} (\eqn{A=1} for active treatment),
#' outcome \eqn{Y}, and binary study indicator \eqn{S} of whether the
#' observation is from the randomized controlled trial \eqn{S=1} or from the
#' external data \eqn{S=0}. If both studies are observational, then \eqn{S=1}
#' should be the reference study.
#' @param S The column name of the \eqn{S} node in \code{data}.
#' @param W The column names of the \eqn{W} node in \code{data}.
#' @param A The column name of the \eqn{A} node in \code{data}.
#' @param Y The column name of the \eqn{Y} node in \code{data}.
#' @param family A character string specifying the family of the outcome
#' \eqn{Y}. Currently only \code{"gaussian"} and \code{"binomial"} are
#' supported.
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
#' \eqn{\tilde{g}_\Delta(1\mid W,A)=\mathbb{P}(\Delta=1\mid W,A)}.
#' \code{"glm"} for main-term linear model, \code{"glmnet"} for lasso,
#' \code{"sl3"} for default super learner, or a \code{list} of \code{sl3}
#' learners. Default is \code{"glmnet"}. \eqn{\Delta} is the indicator of
#' missing outcome. \code{1} - observed, \code{0} - missing.
#' Only applicable when there are missing outcomes.
#' @param theta_tilde_method The method to estimate the nuisance function
#' \eqn{\tilde{\theta}(W)=\mathbb{E}(Y\mid W)}.
#' \code{"glm"} for main-term linear model, \code{"glmnet"} for lasso,
#' \code{"sl3"} for default super learner, or a \code{list} of \code{sl3}
#' learners. Default is \code{"glmnet"}.
#' @param bias_working_model The working model for the bias estimand.
#' Either \code{"glmnet"} for main-term lasso-based working model or
#' \code{"HAL"} for highly adaptive lasso-based working model. Recommended to
#' use \code{"HAL"}.
#' @param bias_working_model_formula An optional formula to specify the
#' working model for the bias estimand. If provided, this will override the
#' \code{bias_working_model} argument.
#' @param pooled_working_model The working model for the pooled-ATE estimand.
#' Either \code{"glmnet"} for main-term lasso-based working model or
#' \code{"HAL"} for highly adaptive lasso-based working model. Recommended to
#' use \code{"HAL"}.
#' @param pooled_working_model_formula An optional formula to specify the
#' working model for the pooled-ATE estimand. If provided, this will override
#' the \code{pooled_working_model} argument.
#' @param cross_fit_nuisance A logical indicating whether to use cross-fitting
#' for nuisance function estimation. Default is \code{TRUE}.
#' @param v_folds The number of folds for cross-validation (whenever necessary).
#' @param stratify A logical indicating whether to stratify the cross-validation
#' fold regime on binary variables. Default is \code{TRUE}.
#' @param g_bounds A numeric vector of lower and upper bounds for the
#' treatment mechanism. The first element is the lower bound, and the second
#' element is the upper bound.
#' @param Pi_bounds A numeric vector of lower and upper bounds for the
#' trial enrollment probabilities. The first element is the lower bound,
#' and the second element is the upper bound.
#' @param theta_bounds A numeric vector of lower and upper bounds for the
#' conditional mean of outcome given baseline covariates and treatment.
#' The first element is the lower bound, and the second element is the upper
#' bound.
#' @param target_gwt If \code{TRUE}, the treatment mechanism is moved from the
#' denominator of the clever covariate to the weight when fitting the TMLE
#' submodel.
#' @param verbose A logical indicating whether to print out the progress.
#' Default is \code{TRUE}.
#' @param max_iter Maximum number of iterations for the iterative targeting
#' procedure. Default is 50.
#' @param target_method The targeting method for the working model parameters.
#' Either \code{"tmle"} for TMLE-type targeting or \code{"relaxed"} for
#' relaxed targeting. Default is \code{"tmle"}.
#' @param eic_method What to do if the information matrix is near singular.
#' Either \code{"svd_pseudo_inv"} for SVD-based pseudo-inverse, or \code{"diag"}
#' for adding a small ridge penalty to the diagonal of the information matrix.
#' Default is \code{"svd_pseudo_inv"}.
#' @param alpha Significance level for the confidence interval.
#' Default is 0.05.
#' @param enumerate_basis_args A \code{list} of arguments to be passed to
#' \code{enumerate_basis} function when constructing the HAL working model
#' basis.
#' @param browse A logical indicating whether to enter the browser for
#' debugging. Default is \code{FALSE}.
#'
#' @returns A \code{list} containing the following elements:
#' \item{est}{The estimated average treatment effect;}
#' \item{lower}{The lower bound of the \eqn{(1-\alpha)\%} confidence interval
#' for the average treatment effect;}
#' \item{upper}{The upper bound of the \eqn{(1-\alpha)\%} confidence interval
#' for the average treatment effect;}
#' \item{psi_pound_est}{The estimated bias parameter;}
#' \item{psi_pound_lower}{The lower bound of the \eqn{(1-\alpha)\%} confidence
#' interval for the bias parameter;}
#' \item{psi_pound_upper}{The upper bound of the \eqn{(1-\alpha)\%} confidence
#' interval for the bias parameter;}
#' \item{psi_tilde_est}{The estimated pooled average treatment effect;}
#' \item{psi_tilde_lower}{The lower bound of the \eqn{(1-\alpha)\%} confidence
#' interval for the pooled average treatment effect;}
#' \item{psi_tilde_upper}{The upper bound of the \eqn{(1-\alpha)\%} confidence
#' interval for the pooled average treatment effect;}
#' \item{eic}{A numeric vector of the estimated efficient influence curve;}
#' \item{tau_A}{A \code{list} containing the working model fit for the
#' pooled-ATE estimand;}
#' \item{tau_S}{A \code{list} containing the working model fit for the bias
#' estimand.}
atmle <- function(data,
                  S,
                  W,
                  A,
                  Y,
                  family,
                  theta_method = "glmnet",
                  Pi_method = "glmnet",
                  g_method = "glmnet",
                  g_delta_method = "glmnet",
                  theta_tilde_method = "glmnet",
                  bias_working_model = "HAL",
                  bias_working_model_formula = NULL,
                  pooled_working_model = "HAL",
                  pooled_working_model_formula = NULL,
                  cross_fit_nuisance = TRUE,
                  v_folds = NULL,
                  stratify = TRUE,
                  g_bounds = NULL,
                  Pi_bounds = NULL,
                  theta_bounds = NULL,
                  target_gwt = TRUE,
                  verbose = TRUE,
                  max_iter = 50,
                  target_method = "relaxed",
                  eic_method = "svd_pseudo_inv",
                  alpha = 0.05,
                  enumerate_basis_args = list(max_degree = 2,
                                              smoothness_orders = 1,
                                              num_knots = c(20, 5)),
                  browse = FALSE) {

  if (browse) browser()

  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }

  # extract variables ----------------------------------------------------------
  S <- data[[S]]
  W <- data[, W, drop = FALSE]
  A <- data[[A]]
  Y <- data[[Y]]
  delta <- as.integer(!is.na(Y))
  n_eff <- sum(delta)
  n <- nrow(data)

  # data-adaptive bounds
  if (is.null(g_bounds)) g_bounds <- c(5/sqrt(n_eff)/log(n_eff), 1-5/sqrt(n_eff)/log(n_eff))
  if (is.null(Pi_bounds)) Pi_bounds <- c(5/sqrt(n_eff)/log(n_eff), 1-5/sqrt(n_eff)/log(n_eff))
  if (is.null(theta_bounds)) theta_bounds <- c(-Inf, Inf)

  # cross-validation scheme (based on tmle R package)
  if (is.null(v_folds)) {
    if (n_eff <= 30){
      v_folds <- n_eff
    } else if (n_eff <= 500) {
      v_folds <- 20
    } else if (n_eff <= 1000) {
      v_folds <- 10
    } else if (n_eff <= 10000){
      v_folds <- 5
    } else {
      v_folds <- 3
    }
  }

  # define controls_only argument
  controls_only <- all(A[S == 0] == 0)

  # cv.glmnet only works when design matrix has at least 2 columns
  # append dummy column of ones if necessary
  if (ncol(W) == 1) {
    W <- cbind(W, 1)
  }

  # cross fitting schemes
  if (family == "gaussian") {
    if (stratify) {
      if (sum(delta) < n) {
        cv_strata <- paste0(S, "-", A, "-", delta)
      } else {
        cv_strata <- paste0(S, "-", A)
      }
    } else {
      cv_strata <- rep(1, n)
    }
  } else if (family == "binomial") {
    if (stratify) {
      if (sum(delta) < n) {
        cv_strata <- paste0(S, "-", A, "-", Y, "-", delta)
      } else {
        cv_strata <- paste0(S, "-", A, "-", Y)
      }
    } else {
      cv_strata <- rep(1, n)
    }
  }
  suppressWarnings({
    folds <- make_folds(
      n = n, V = v_folds,
      strata_ids = as.integer(factor(cv_strata))
    )
  })
  foldid <- folds2foldvec(folds)
  foldid_obs <- foldid[delta == 1]
  folds_obs <- map(seq(v_folds), function(v) {
    fold_from_foldvec(v = v, folds = foldid_obs)
  })
  foldid_S1 <- foldid[S == 1]
  folds_S1 <- map(seq(v_folds), function(v) {
    fold_from_foldvec(v = v, folds = foldid_S1)
  })
  foldid_S0 <- foldid[S == 0]
  folds_S0 <- map(seq(v_folds), function(v) {
    fold_from_foldvec(v = v, folds = foldid_S0)
  })

  if (sum(delta) < n) {
    # outcome has missing
    if (verbose) cat("learning g(\U0394=1|S,W,A)=P(\U0394=1|S,W,A)...")
    g_delta <- learn_g_delta(S = S,
                             W = W,
                             A = A,
                             delta = delta,
                             method = g_delta_method,
                             folds = folds,
                             g_bounds = g_bounds,
                             cross_fit_nuisance = cross_fit_nuisance)
    if (verbose) cat("Done!\n")

    if (verbose) cat("learning g(\U0394=1|W,A)=P(\U0394=1|W,A)...")
    g_delta_tilde <- learn_g_delta_tilde(W = W,
                                         A = A,
                                         delta = delta,
                                         method = g_delta_method,
                                         folds = folds,
                                         g_bounds = g_bounds,
                                         cross_fit_nuisance = cross_fit_nuisance)
    if (verbose) cat("Done!\n")
  } else {
    # no censoring
    g_delta <- g_delta_tilde <- list(pred = rep(1, length(A)),
                                     A0 = rep(1, length(A)),
                                     A1 = rep(1, length(A)))
  }

  # censoring weights
  weights <- delta/g_delta$pred
  weights_tilde <- delta/g_delta_tilde$pred

  # estimate bias psi_pound ----------------------------------------------------
  # learn nuisance parts
  if (verbose) cat("learning \U03B8(W,A)=E(Y|W,A)...")
  theta_WA <- learn_theta_W(W = as.matrix(cbind(W, A=A)),
                            Y = Y,
                            delta = delta,
                            weights = weights,
                            method = theta_method,
                            folds = folds,
                            folds_obs = folds_obs,
                            family = family,
                            theta_bounds = theta_bounds,
                            cross_fit_nuisance = cross_fit_nuisance)
  if (verbose) cat("Done!\n")

  if (verbose) cat("learning g(1|W)=P(A=1|W)...")
  g1W <- learn_g(S = S,
                 W = W,
                 A = A,
                 method = g_method,
                 controls_only = controls_only,
                 folds = folds,
                 folds_S1 = folds_S1,
                 folds_S0 = folds_S0,
                 g_bounds = g_bounds,
                 cross_fit_nuisance = cross_fit_nuisance)
  if (verbose) cat("Done!\n")

  if (verbose) cat("learning \U03A0(S=1|W,A)=P(S=1|W,A)...")
  Pi <- learn_Pi(g = g1W,
                 A = A,
                 Pi_bounds = Pi_bounds)
  if (verbose) cat("Done!\n")

  # learn working model tau for bias
  if (verbose) cat("learning \U03C4(W,A)=E(Y|S=1,W,A)-E(Y|S=0,W,A)...")
  tau_S <- learn_tau_S(S = S,
                       W = W,
                       A = A,
                       Y = Y,
                       Pi = Pi,
                       theta = theta_WA,
                       g1W = g1W,
                       delta = delta,
                       controls_only = controls_only,
                       method = bias_working_model,
                       v_folds = v_folds,
                       foldid = foldid,
                       min_working_model = NULL,
                       target_gwt = target_gwt,
                       Pi_bounds = Pi_bounds,
                       enumerate_basis_args = enumerate_basis_args,
                       weights = weights,
                       bias_working_model_formula = NULL,
                       verbose = verbose)
  if (verbose) cat("Done!\n")

  # iterative targeting between Pi and beta_S
  cur_iter <- 1
  PnEIC <- Inf
  sn <- 0
  while (cur_iter <= max_iter & abs(PnEIC) > sn) {
    # target Pi
    Pi_and_tau_S <- target_Pi(S = S,
                              W = W,
                              A = A,
                              Y = Y,
                              delta = delta,
                              g1W = g1W$pred,
                              Pi = Pi,
                              theta_WA = theta_WA,
                              tau_S = tau_S,
                              controls_only = controls_only,
                              target_gwt = target_gwt,
                              weights = weights,
                              Pi_bounds = Pi_bounds)
    Pi <- Pi_and_tau_S$Pi; tau_S <- Pi_and_tau_S$tau_S
    psi_pound_eic <- eic_psi_pound_wm(S = S,
                                      Y = Y,
                                      A = A,
                                      g1W = g1W$pred,
                                      theta_WA = theta_WA,
                                      Pi = Pi,
                                      tau_S = tau_S,
                                      weights = weights,
                                      controls_only = controls_only)
    PnEIC <- mean(psi_pound_eic)
    sn <- 0.001*sqrt(var(psi_pound_eic, na.rm = TRUE))/(sqrt(length(Y)) * log(length(Y)))
    if (abs(PnEIC) <= sn) {
      break
    }

    # target beta_S
    tau_S <- target_beta_S(S = S,
                           W = W,
                           A = A,
                           Y = Y,
                           g1W = g1W,
                           Pi = Pi,
                           theta_WA = theta_WA,
                           tau_S = tau_S,
                           weights = weights,
                           controls_only = controls_only,
                           eic_method = eic_method,
                           target_method = target_method)
    psi_pound_eic <- eic_psi_pound_wm(S = S,
                                      Y = Y,
                                      A = A,
                                      g1W = g1W$pred,
                                      theta_WA = theta_WA,
                                      Pi = Pi,
                                      tau_S = tau_S,
                                      weights = weights,
                                      controls_only = controls_only)
    PnEIC <- mean(psi_pound_eic)
    sn <- sqrt(var(psi_pound_eic, na.rm = TRUE))/(sqrt(length(Y)) * log(length(Y)))
    cur_iter <- cur_iter + 1
    if (verbose) print(PnEIC)
  }

  if (controls_only) {
    psi_pound_est <- mean((1-Pi$A0)*tau_S$cate_W0)
  } else {
    psi_pound_est <- mean((1-Pi$A0)*tau_S$cate_W0-(1-Pi$A1)*tau_S$cate_W1)
  }

  # estimate pooled-ATE psi_tilde ----------------------------------------------
  if (verbose) cat("learning \U03B8\U0303(W)=E(Y|W)...")
  theta_W <- learn_theta_W(W = W,
                           Y = Y,
                           delta = delta,
                           weights = weights_tilde,
                           method = theta_tilde_method,
                           folds = folds,
                           folds_obs = folds_obs,
                           family = family,
                           theta_bounds = theta_bounds,
                           cross_fit_nuisance = cross_fit_nuisance)
  if (verbose) cat("Done!\n")

  if (verbose) cat("learning \U03A4(W)=E(Y|W,A=1)-E(Y|W,A=0)...")
  tau_A <- learn_tau_A(W = W,
                       A = A,
                       Y = Y,
                       g1W = g1W$pred,
                       delta = delta,
                       theta = theta_W,
                       method = pooled_working_model,
                       foldid = foldid,
                       weights = weights_tilde,
                       enumerate_basis_args = enumerate_basis_args,
                       pooled_working_model_formula = NULL,
                       target_method = target_method,
                       verbose = verbose)
  if (verbose) cat("Done!\n\n")

  if (verbose) cat("targeting beta_A...")
  tau_A <- target_beta_A(A = A,
                         Y = Y,
                         g1W = g1W$pred,
                         theta_W = theta_W,
                         tau_A = tau_A,
                         eic_method = eic_method,
                         target_method = target_method)
  if (verbose) cat("Done!\n\n")

  # estimates
  psi_tilde_est <- mean(tau_A$cate_W)
  psi_tilde_eic <- eic_psi_tilde_wm(Y = Y,
                                    A = A,
                                    g1W = g1W$pred,
                                    theta_W = theta_W,
                                    tau_A = tau_A,
                                    weights = weights_tilde,
                                    eic_method = eic_method)

  # final estimates ------------------------------------------------------------
  psi_pound_se <- sqrt(var(psi_pound_eic, na.rm = TRUE)/n)
  psi_pound_lower <- psi_pound_est+qnorm(alpha/2)*psi_pound_se
  psi_pound_upper <- psi_pound_est+qnorm(1-alpha/2)*psi_pound_se
  psi_tilde_se <- sqrt(var(psi_tilde_eic, na.rm = TRUE)/n)
  psi_tilde_lower <- psi_tilde_est+qnorm(alpha/2)*psi_tilde_se
  psi_tilde_upper <- psi_tilde_est+qnorm(1-alpha/2)*psi_tilde_se
  est <- psi_tilde_est-psi_pound_est
  eic <- psi_tilde_eic-psi_pound_eic
  se <- sqrt(var(eic, na.rm = TRUE)/n)
  lower <- est+qnorm(alpha/2)*se
  upper <- est+qnorm(1-alpha/2)*se

  results <- list(est = est,
                  lower = lower,
                  upper = upper,
                  psi_pound_est = psi_pound_est,
                  psi_pound_lower = psi_pound_lower,
                  psi_pound_upper = psi_pound_upper,
                  psi_tilde_est = psi_tilde_est,
                  psi_tilde_lower = psi_tilde_lower,
                  psi_tilde_upper = psi_tilde_upper,
                  eic = eic,
                  tau_A = tau_A,
                  tau_S = tau_S)

  if (verbose) {
    cat("Pooled ATE: ", signif(results$psi_tilde_est, 5), " (", signif(results$psi_tilde_lower, 5), ", ", signif(results$psi_tilde_upper, 5), ")\n", sep = "")
    cat("Bias: ", signif(results$psi_pound_est, 5), " (", signif(results$psi_pound_lower, 5), ", ", signif(results$psi_pound_upper, 5), ")\n", sep = "")
    cat("Bias-corrected ATE: ", signif(results$est, 5), " (", signif(results$lower, 5), ", ", signif(results$upper, 5), ")\n", sep = "")
  }

  return(results)
}
