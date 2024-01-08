library(ggplot2)
library(devtools)
library(EScvtmle)
library(sl3)
load_all()

`%+%` <- function(a, b) paste0(a, b)

#' @param B Number of runs for each sample size
#' @param n_min Minimum of sample size
#' @param n_max Maximum of sample size
#' @param n_step Step size of sample size increase
#' @param bA True ATE
#' @param nuisance_method Fitting method for nuisance parts, "lasso" or "HAL"
#' @param working_model Working model types, "lasso" or "HAL"
run_sim_n_increase <- function(B,
                               n_min,
                               n_max,
                               n_step,
                               bA,
                               bias,
                               controls_only,
                               g_rct,
                               num_covs,
                               var_method,
                               family,
                               method,
                               type,
                               verbose = TRUE) {

  n_seq <- seq(n_min, n_max, n_step)

  # results
  all_psi_est <- vector(mode = "list", length = length(n_seq))
  all_psi_coverage <- vector(mode = "list", length = length(n_seq))
  all_psi_ci_lower <- vector(mode = "list", length = length(n_seq))
  all_psi_ci_upper <- vector(mode = "list", length = length(n_seq))
  all_escvtmle_prop_selected <- vector(mode = "list", length = length(n_seq))

  for (i in 1:length(n_seq)) {
    n <- n_seq[i]
    print(n)

    psi_est <- vector(length = B)
    psi_coverage <- vector(length = B)
    psi_ci_lower <- vector(length = B)
    psi_ci_upper <- vector(length = B)
    escvtmle_prop_selected <- vector(length = B)

    for (j in 1:B) {
      # simulate data
      data <- sim_four_covs(bA, n, 0.2, g_rct, bias, controls_only)
      S_node <- 1
      W_node <- 2:5
      A_node <- 6
      Y_node <- 7

      # fit
      res <- NULL
      if (method == "oracle_atmle") {
        res <- atmle_oracle(data,
                            S_node = S_node,
                            W_node = W_node,
                            A_node = A_node,
                            Y_node = Y_node,
                            controls_only = controls_only,
                            bias = bias,
                            atmle_pooled = TRUE,
                            var_method = var_method,
                            nuisance_method = nuisance_method,
                            working_model = working_model,
                            g_rct = g_rct,
                            verbose = FALSE)
      } else if (method == "oracle_atmle_tmle") {
        res <- atmle_oracle(data,
                            S_node = S_node,
                            W_node = W_node,
                            A_node = A_node,
                            Y_node = Y_node,
                            controls_only = controls_only,
                            bias = bias,
                            atmle_pooled = FALSE,
                            var_method = var_method,
                            nuisance_method = nuisance_method,
                            working_model = working_model,
                            g_rct = g_rct,
                            verbose = FALSE)
      } else if (method == "atmle") {
        res <- atmle(data = data,
                     S_node = S_node,
                     W_node = W_node,
                     A_node = A_node,
                     Y_node = Y_node,
                     controls_only = controls_only,
                     family = family,
                     atmle_pooled = TRUE,
                     var_method = var_method,
                     theta_method = "glm",
                     Pi_method = "glmnet",
                     g_method = "glm",
                     theta_tilde_method = "glm",
                     Q_method = "glm",
                     bias_working_model = "glmnet",
                     pooled_working_model = "glmnet",
                     g_rct = g_rct,
                     verbose = FALSE)
      } else if (method == "atmle_tmle") {
        res <- atmle(data = data,
                     S_node = S_node,
                     W_node = W_node,
                     A_node = A_node,
                     Y_node = Y_node,
                     controls_only = controls_only,
                     family = family,
                     atmle_pooled = FALSE,
                     var_method = var_method,
                     theta_method = "sl3",
                     Pi_method = "glmnet",
                     g_method = "sl3",
                     theta_tilde_method = "sl3",
                     Q_method = "sl3",
                     bias_working_model = "glmnet",
                     pooled_working_model = "glmnet",
                     g_rct = g_rct,
                     verbose = FALSE)
      } else if (method == "escvtmle") {
        covariates <- NULL
        if (num_covs == 2) {
          covariates <- c("W1", "W2")
        } else if (num_covs == 4) {
          covariates <- c("W1", "W2", "W3", "W4")
        }
        tmp <- ES.cvtmle(txinrwd = !controls_only,
                         data = data,
                         study = "S",
                         covariates = covariates,
                         treatment_var = "A",
                         treatment = 1,
                         outcome = "Y",
                         pRCT = g_rct,
                         family = "gaussian",
                         Q.SL.library = c("SL.glm"),
                         g.SL.library = c("SL.glm"),
                         Q.discreteSL = TRUE,
                         g.discreteSL = TRUE,
                         V = 5)
        res <- list(est = tmp$ATE$b2v,
                    lower = as.numeric(tmp$CI$b2v[1]),
                    upper = as.numeric(tmp$CI$b2v[2]))
        escvtmle_prop_selected[j] <- tmp$proportionselected$b2v
      } else if (method == "tmle") {
        res <- nonparametric(data = data,
                             S_node = S_node,
                             W_node = W_node,
                             A_node = A_node,
                             Y_node = Y_node,
                             controls_only = controls_only,
                             family = family,
                             atmle_pooled = TRUE,
                             var_method = var_method,
                             theta_method = nuisance_method,
                             Pi_method = nuisance_method,
                             g_method = nuisance_method,
                             theta_tilde_method = nuisance_method,
                             Q_method = nuisance_method,
                             bias_working_model = working_model,
                             pooled_working_model = working_model,
                             g_rct = g_rct,
                             verbose = FALSE)
      } else if (method == "rct_only") {
        res <- rct_only(data,
                        S_node = S_node,
                        W_node = W_node,
                        A_node = A_node,
                        Y_node = Y_node,
                        nuisance_method = nuisance_method,
                        g_rct = g_rct,
                        verbose = FALSE)
      }

      if (res$lower <= bA & res$upper >= bA) {
        if (verbose) print("psi covered")
        psi_coverage[j] <- 1
      } else {
        if (verbose) print("psi not covered")
        psi_coverage[j] <- 0
      }

      psi_est[j] <- res$est
      psi_ci_lower[j] <- res$lower
      psi_ci_upper[j] <- res$upper
    }

    all_psi_est[[i]] <- psi_est
    all_psi_coverage[[i]] <- psi_coverage
    all_psi_ci_lower[[i]] <- psi_ci_lower
    all_psi_ci_upper[[i]] <- psi_ci_upper

    if (method == "escvtmle") {
      all_escvtmle_prop_selected[[i]] <- escvtmle_prop_selected
    }
  }

  return(list(all_psi_est = all_psi_est,
              all_psi_coverage = all_psi_coverage,
              all_psi_ci_lower = all_psi_ci_lower,
              all_psi_ci_upper = all_psi_ci_upper,
              escvtmle_prop_selected = escvtmle_prop_selected))
}

# A-TMLE oracle, working model always true
atmle_oracle <- function(data,
                         S_node,
                         W_node,
                         A_node,
                         Y_node,
                         controls_only,
                         bias,
                         family="gaussian",
                         atmle_pooled = TRUE,
                         var_method = "ic",
                         nuisance_method="glmnet",
                         working_model="glmnet",
                         g_rct=0.5,
                         verbose=TRUE) {

  # define nodes
  S <- data[, S_node]
  W <- data[, W_node]
  A <- data[, A_node]
  Y <- data[, Y_node]
  n <- nrow(data)

  # estimate bias psi_pound ----------------------------------------------------
  # learn nuisance parts
  if (verbose) print("learning E(Y|W,A)")
  theta <- learn_theta(W, A, Y, controls_only, family = "gaussian", nuisance_method)

  if (verbose) print("learning P(S=1|W,A)")
  Pi <- learn_Pi(S, W, A, controls_only, nuisance_method)

  if (verbose) print("learning P(A=1|W)")
  g <- learn_g(S, W, A, g_rct, nuisance_method)

  # learn working model tau for bias
  if (verbose) print("learning E(Y|S,W,A)")
  tau <- tau_oracle(S, W, A, Y, Pi, theta, controls_only, bias)

  # TMLE to target Pi
  if (verbose) print("targeting P(S=1|W,A)")
  for (i in 1:1) {
    # TMLE to target Pi
    Pi_star <- Pi_tmle(S, W, A, g, tau, Pi, controls_only)
    Pi <- Pi_star
  }

  psi_pound_est <- NULL
  if (controls_only) {
    psi_pound_est <- mean((1-Pi$A0)*tau$A0)
  } else {
    psi_pound_est <- mean((1-Pi$A0)*tau$A0-(1-Pi$A1)*tau$A1)
  }

  # estimate pooled-ATE psi_tilde ----------------------------------------------
  psi_tilde <- NULL
  psi_tilde_est <- NULL
  psi_tilde_eic <- NULL
  if (atmle_pooled) {
    # use atmle for pooled-ATE
    if (verbose) print("learning E(Y|W)")
    theta_tilde <- learn_theta_tilde(W, Y, family = "gaussian", nuisance_method)

    if (verbose) print("learning psi_tilde")
    psi_tilde <- learn_psi_tilde(W, A, Y, g, theta_tilde, "glmnet")

    # estimates
    psi_tilde_est <- mean(psi_tilde$pred)
    psi_tilde_eic <- get_eic_psi_tilde(psi_tilde, g, theta_tilde, Y, A, n)
  } else {
    # use regular TMLE for pooled-ATE
    Q <- learn_Q(W, A, Y, method = nuisance_method)
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

# oracle tau working model
tau_oracle <- function(S, W, A, Y, Pi, theta, controls_only, bias) {
  pred <- NULL
  A1 <- numeric(length = length(A))
  A0 <- numeric(length = length(A))
  X <- NULL
  X_A1_counter <- NULL
  X_A0_counter <- NULL
  x_basis <- NULL
  x_basis_A1 <- NULL
  x_basis_A0 <- NULL
  coefs <- NULL

  pred <- numeric(length = length(A))

  # design matrix, both treated and controls (X: W*I(A == 1), I(A == 1), W*I(A == 0))
  X <- cbind(W * A, A, W * (1-A))

  # counterfactual design matrices
  zero_W <- matrix(0, nrow = nrow(W), ncol = ncol(W))
  X_A1_counter <- cbind(1, W, 1, zero_W)
  X_A0_counter <- cbind(1, zero_W, 0, W)

  # oracle
  non_zero <- NULL
  if (is.numeric(bias)) {
    non_zero <- c(1)
    coefs <- -c(bias)
  } else if (bias == "param_simple") {
    non_zero <- c(1, 7)
    coefs <- -c(0, 0.8) # W1
  } else if (bias == "param_complex") {
    non_zero <- c(1, 7, 8, 9)
    coefs <- -c(0, 0.6, 1.5, 1.4) # W1, W2, W3
  }

  # selected counterfactual bases
  x_basis <- as.matrix(cbind(1, X)[, non_zero, drop = FALSE])
  x_basis_A1 <- as.matrix(X_A1_counter[, non_zero, drop = FALSE])
  x_basis_A0 <- as.matrix(X_A0_counter[, non_zero, drop = FALSE])

  # predictions
  A1 <- as.numeric(x_basis_A1 %*% matrix(coefs))
  A0 <- as.numeric(x_basis_A0 %*% matrix(coefs))
  pred[A == 1] <- A1[A == 1]
  pred[A == 0] <- A0[A == 0]

  return(list(A1 = A1,
              A0 = A0,
              x_basis = x_basis,
              x_basis_A1 = x_basis_A1,
              x_basis_A0 = x_basis_A0,
              pred = pred,
              coefs = coefs))
}
