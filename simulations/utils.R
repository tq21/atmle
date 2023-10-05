library(ggplot2)
library(devtools)
library(EScvtmle)
library(sl3)
load_all()

`%+%` <- function(a, b) paste0(a, b)

generate_realistic_data <- function(ate, n_rct, n_rwd, g_rct, bias, controls_only) {
  # total number of observations needed for RCT and RWD combined
  n <- max(n_rct, n_rwd) * 2  # ensure we get at least n_rct and n_rwd after applying the criteria

  UY <- rnorm(n, 0, 0.5)

  # baseline covariates
  W1 <- rnorm(n, 0, 1)
  W2 <- rnorm(n, 0, 1)
  W3 <- rnorm(n, 0, 1)
  W4 <- rnorm(n, 0, 1)

  # probability being in RCT based on inclusion criteria
  #prob_rct <- plogis(0.7+0.3*W1+1.1*W2-0.9*W3-0.9*W4)
  prob_rct <- 0.2
  in_rct <- rbinom(n, 1, prob_rct)

  # ensure exactly n_rct and n_rwd samples
  rct_indices <- sample(which(in_rct == 1), n_rct, replace = TRUE)
  rwd_indices <- sample(which(in_rct == 0), n_rwd, replace = TRUE)

  # assign treatment in RCT
  A_rct <- rbinom(n_rct, 1, g_rct)

  # assign treatment in RWD based on doctor's decision rule
  #decision_rule <- plogis(0.7*W1[rwd_indices]-1.1*W2[rwd_indices])
  decision_rule <- 0.4
  A_rwd <- if (controls_only) rep(0, n_rwd) else rbinom(n_rwd, 1, decision_rule)

  # outcome Y for RCT data without bias
  Y_rct <- 2.1+0.8*W1[rct_indices]+2.5*W2[rct_indices]-3.1*W3[rct_indices]+0.9*W4[rct_indices]+ate*A_rct+UY[rct_indices]

  # bias term for RWD data
  b <- NULL
  if (is.numeric(bias)) {
    if (bias == 0) {
      b <- rep(0, n_rwd)
    } else {
      b <- rnorm(n_rwd, bias, 0.1)*(1-A_rwd)
    }
  } else if (bias == "param_simple") {
    b <- (1.9+2.3*W1[rwd_indices]+rnorm(n_rwd, 0, 0.1))*(1-A_rwd)
  } else if (bias == "param_complex") {
    b <- (1.9+2.3*W1[rwd_indices]+1.2*W2[rwd_indices]+0.5*W3[rwd_indices]+rnorm(n_rwd, 0, 0.1))*(1-A_rwd)
  }

  # outcome Y for RWD data with bias
  Y_rwd <- 2.1+0.8*W1[rwd_indices]+2.5*W2[rwd_indices]-3.1*W3[rwd_indices]+0.9*W4[rwd_indices]+ate*A_rwd+b+UY[rwd_indices]

  # data frames for RCT and RWD
  rct_data <- data.frame(
    S = rep(1, n_rct),
    W1 = W1[rct_indices],
    W2 = W2[rct_indices],
    W3 = W3[rct_indices],
    W4 = W4[rct_indices],
    A = A_rct,
    Y = Y_rct
  )

  rwd_data <- data.frame(
    S = rep(0, n_rwd),
    W1 = W1[rwd_indices],
    W2 = W2[rwd_indices],
    W3 = W3[rwd_indices],
    W4 = W4[rwd_indices],
    A = A_rwd,
    Y = Y_rwd
  )

  return(rbind(rct_data, rwd_data))
}

#' @param B Number of simulations to run
#' @param n Sample size of each data
#' @param bA True ATE
#' @param nuisance_method Fitting method for nuisance parts, "lasso" or "HAL"
#' @param working_model Working model types, "lasso" or "HAL"
run_sim <- function(B,
                    n_rct,
                    n_rwd,
                    ate,
                    bias,
                    nuisance_method,
                    working_model,
                    g_rct,
                    controls_only,
                    verbose=TRUE,
                    method="atmle") {
  # results
  psi_coverage <- rep(NA, B)
  psi_est <- rep(NA, B)
  psi_ci_lower <- rep(NA, B)
  psi_ci_upper <- rep(NA, B)
  escvtmle_prop_selected <- rep(NA, B)

  for (i in 1:B) {
    if (verbose) print(i)

    # simulate data
    data <- generate_realistic_data(ate, n_rct, n_rwd, g_rct, bias, controls_only)

    # fit
    res <- NULL
    if (method == "atmle") {
      res <- atmle(data,
                   S_node = 1,
                   W_node = c(2, 3, 4, 5),
                   A_node = 6,
                   Y_node = 7,
                   controls_only = controls_only,
                   atmle_pooled = TRUE,
                   nuisance_method = nuisance_method,
                   working_model = working_model,
                   g_rct = g_rct,
                   verbose = FALSE)
    } else if (method == "atmle_tmle") {
      res <- atmle(data,
                   S_node = 1,
                   W_node = c(2, 3, 4, 5),
                   A_node = 6,
                   Y_node = 7,
                   controls_only = controls_only,
                   atmle_pooled = FALSE,
                   nuisance_method = nuisance_method,
                   working_model = working_model,
                   g_rct = g_rct,
                   verbose = FALSE)
    } else if (method == "atmle psi_tilde") {
      res <- psi_tilde_only_atmle(data,
                                  S_node = 1,
                                  W_node = c(2, 3, 4, 5),
                                  A_node = 6,
                                  Y_node = 7,
                                  nuisance_method = nuisance_method,
                                  working_model = working_model,
                                  g_rct = g_rct,
                                  verbose = FALSE)
    } else if (method == "tmle psi_tilde") {
      res <- psi_tilde_only_tmle(data,
                                 S_node = 1,
                                 W_node = c(2, 3, 4, 5),
                                 A_node = 6,
                                 Y_node = 7,
                                 nuisance_method = nuisance_method,
                                 working_model = working_model,
                                 g_rct = g_rct,
                                 verbose = FALSE)
    } else if (method == "escvtmle") {
      tmp <- ES.cvtmle(txinrwd = !controls_only,
                       data = data,
                       study = "S",
                       covariates = c("W1", "W2", "W3", "W4"),
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
      escvtmle_prop_selected[i] <- tmp$proportionselected$b2v
    } else if (method == "rct_only") {
      res <- rct_only(data,
                      S_node = 1,
                      W_node = c(2, 3, 4, 5),
                      A_node = 6,
                      Y_node = 7,
                      nuisance_method = nuisance_method,
                      g_rct = g_rct,
                      verbose = FALSE)
    } else if (method == "tmle") {
      res <- nonparametric(data,
                           S_node = 1,
                           W_node = c(2, 3, 4, 5),
                           A_node = 6,
                           Y_node = 7,
                           nuisance_method = nuisance_method,
                           working_model = working_model,
                           g_rct = g_rct,
                           verbose = FALSE)
    }

    if (res$lower <= ate & res$upper >= ate) {
      psi_coverage[i] <- 1
      if (verbose) print("psi covered")
    } else {
      psi_coverage[i] <- 0
      if (verbose) print("psi not covered")
    }

    psi_est[i] <- res$est
    psi_ci_lower[i] <- res$lower
    psi_ci_upper[i] <- res$upper
  }

  return(list(psi_est = psi_est,
              psi_coverage = psi_coverage,
              psi_ci_lower = psi_ci_lower,
              psi_ci_upper = psi_ci_upper,
              escvtmle_prop_selected = escvtmle_prop_selected))
}

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
                               nuisance_method,
                               working_model,
                               g_rct,
                               controls_only,
                               verbose=TRUE,
                               method="atmle") {

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
      n_rwd <- n
      n_rct <- round(n * 0.1)
      data <- generate_realistic_data(bA, n_rct, n_rwd, g_rct, bias, controls_only)

      # fit
      res <- NULL
      if (method == "atmle") {
        res <- atmle(data,
                     S_node = 1,
                     W_node = c(2, 3, 4, 5),
                     A_node = 6,
                     Y_node = 7,
                     atmle_pooled = TRUE,
                     nuisance_method=nuisance_method,
                     working_model=working_model,
                     g_rct = g_rct,
                     verbose = FALSE)
      } else if (method == "atmle_tmle") {
        res <- atmle(data,
                     S_node = 1,
                     W_node = c(2, 3, 4, 5),
                     A_node = 6,
                     Y_node = 7,
                     atmle_pooled = FALSE,
                     nuisance_method=nuisance_method,
                     working_model=working_model,
                     g_rct = g_rct,
                     verbose = FALSE)
      } else if (method == "atmle psi_tilde") {
        res <- psi_tilde_only_atmle(data,
                                    S_node = 1,
                                    W_node = c(2, 3, 4, 5),
                                    A_node = 6,
                                    Y_node = 7,
                                    nuisance_method=nuisance_method,
                                    working_model=working_model,
                                    g_rct = g_rct,
                                    verbose = FALSE)
      } else if (method == "tmle psi_tilde") {
        res <- psi_tilde_only_tmle(data,
                                   S_node = 1,
                                   W_node = c(2, 3, 4, 5),
                                   A_node = 6,
                                   Y_node = 7,
                                   nuisance_method=nuisance_method,
                                   working_model=working_model,
                                   g_rct = g_rct,
                                   verbose = FALSE)
      } else if (method == "escvtmle") {
        #data$S <- 1 - data$S
        tmp <- ES.cvtmle(txinrwd = !controls_only,
                         data = data,
                         study = "S",
                         covariates = c("W1", "W2", "W3", "W4"),
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
        res <- nonparametric(data,
                             S_node = 1,
                             W_node = c(2, 3, 4, 5),
                             A_node = 6,
                             Y_node = 7,
                             nuisance_method = nuisance_method,
                             working_model = working_model,
                             g_rct = g_rct,
                             verbose = FALSE)
      } else if (method == "rct_only") {
        res <- rct_only(data,
                        S_node = 1,
                        W_node = c(2, 3, 4, 5),
                        A_node = 6,
                        Y_node = 7,
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
