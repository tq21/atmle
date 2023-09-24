library(ggplot2)
library(devtools)
library(EScvtmle)
load_all()

`%+%` <- function(a, b) paste0(a, b)

generate_data <- function(N, bA, bias, pRCT){
  # U
  UY <- rnorm(N, 0, 1)

  # S
  #S <- rbinom(N, 1, 0.5)
  #S <- c(rep(1, 100), rep(0, N - 100)) # this does not work

  # W
  # W1 <- rnorm(N)
  # W2 <- rnorm(N)
  # W3 <- rnorm(N)
  # W4 <- rnorm(N)
  W1 <- runif(N)
  W2 <- runif(N)
  W3 <- runif(N)
  W4 <- runif(N)
  # W1 <- runif(N)
  # W2 <- runif(N, -2, 4)
  # W3 <- rbinom(N, 1, 0.7)
  # W4 <- rgamma(N, 3, 1)

  # S
  S <- rbinom(N, 1, 0.5)

  # A
  A <- vector(length = N)
  A <- rbinom(N, 1, pRCT)

  #A[S == 0] <- rbinom(N - sum(S), 1, 0.67)
  #A[S == 0] <- rbinom(N - sum(S), 1, plogis(-0.5*W1+0.9*W2-1.2*W3+0.3*W4))
  #A[S == 0] <- rbinom(N - sum(S), 1, plogis(-5*W1^2-5*W2^3+0.9*W2-1.2*W3+0.3*W4))
  #A[S == 0] <- rbinom(N - sum(S), 1, 0.05)
  #A[S == 1] <- rbinom(sum(S), 1, pRCT) # rct

  # Y
  Y <- vector(length = N)
  if (is.numeric(bias)) {
    b <- bias
    # if (bias != 0) {
    #   b <- rnorm(N, bias, 1)
    # }
    Y_S0 <- 0.3+0.5*W1+0.1*W2+0.3*W3-0.5*W4+bA*A-b+UY # RWD has bias
    Y_S1 <- 0.3+0.5*W1+0.1*W2+0.3*W3-0.5*W4+bA*A+UY
    Y[S == 0] <- Y_S0[S == 0]
    Y[S == 1] <- Y_S1[S == 1]
  } else if (bias == "param_simple") {
    b <- 1.3*W1+1.9+rnorm(N)
    Y_S0 <- 0.3+bA*A+0.5*W1+0.1*W2+0.3*W3-0.5*W4-b*A+UY
    Y_S1 <- 0.3+bA*A+0.5*W1+0.1*W2+0.3*W3-0.5*W4+UY
    Y[S == 0] <- Y_S0[S == 0]
    Y[S == 1] <- Y_S1[S == 1]
  } else if (bias == "param_complex") {
    b <- 1.3*W1+1.1*W2+0.9*W3-1.1*W4+1.2+rnorm(N)
    Y_S0 <- 0.3+0.5*W1+0.1*W2+0.3*W3-0.5*W4+A*bA-b*A+UY
    Y_S1 <- 0.3+0.5*W1+0.1*W2+0.3*W3-0.5*W4+A*bA+UY
    Y[S == 0] <- Y_S0[S == 0]
    Y[S == 1] <- Y_S1[S == 1]
  } else if (bias == "HAL") {
    Y_S0 <- 0.3+bA*A+0.5*W1+0.1*W2+0.3*W3-0.5*W4+UY+
      -0.5
    Y_S1 <- 0.3+bA*A+0.5*W1+0.1*W2+0.3*W3-0.5*W4+UY
    Y[S == 0] <- Y_S0[S == 0]
    Y[S == 1] <- Y_S1[S == 1]
  } else if (bias == "hard") {
    Y_S0 <- 0.3+bA*A+0.5*W1+0.1*W2+0.3*W3-0.5*W4+UY
      1.3*W1+0.9*W2^2+0.7*W3^3+1.5*W4^2+1.5*W1*W2+1.9*W1*W2*W4^2-0.6
    Y_S1 <- 0.3+bA*A+0.5*W1+0.1*W2+0.3*W3-0.5*W4+UY
    Y[S == 0] <- Y_S0[S == 0]
    Y[S == 1] <- Y_S1[S == 1]
  } else if (bias == "sinusoidal") {
    Y_S0 <- 0.3+bA*A+0.5*W1^2+0.1*W2^3+0.3*W3-0.5*W4+UY+
      3.8*W3*as.numeric(W1 < 0.5)* sin(pi/2*abs(W1))+
      4*as.numeric(W2 > 0.7)*cos(pi/2*abs(W1))+
      2.1*W4*sin(pi*W4)+3.2*W3*cos(abs(W4-W3))
    Y_S1 <- 0.3+bA*A+0.5*W1+0.1*W2+0.3*W3-0.5*W4+UY
    Y[S == 0] <- Y_S0[S == 0]
    Y[S == 1] <- Y_S1[S == 1]
  } else if (bias == "test 2") {
    Y_S0 <- 0.3+bA*A+0.5*W1+0.1*W2+0.3*W3-0.5*W4+UY+100
    Y_S1 <- 0.3+bA*A+0.5*W1+0.1*W2+0.3*W3-0.5*W4+UY
    Y[S == 0] <- Y_S0[S == 0]
    Y[S == 1] <- Y_S1[S == 1]
  }

  # data
  data <- data.frame(S, W1, W2, W3, W4, A, Y)

  return(data)
}

#' @param B Number of simulations to run
#' @param n Sample size of each data
#' @param bA True ATE
#' @param nuisance_method Fitting method for nuisance parts, "lasso" or "HAL"
#' @param working_model Working model types, "lasso" or "HAL"
run_sim <- function(B,
                    n,
                    bA,
                    bias,
                    nuisance_method,
                    working_model,
                    pRCT,
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
    data <- generate_data(n, bA, bias, pRCT)

    # fit
    res <- NULL
    if (method == "atmle") {
      res <- atmle(data,
                   S_node = 1,
                   W_node = c(2, 3, 4, 5),
                   A_node = 6,
                   Y_node = 7,
                   nuisance_method=nuisance_method,
                   working_model=working_model,
                   p_rct = pRCT,
                   verbose = FALSE)
    } else if (method == "atmle_tmle") {
      res <- atmle_tmle(data,
                        S_node = 1,
                        W_node = c(2, 3, 4, 5),
                        A_node = 6,
                        Y_node = 7,
                        nuisance_method=nuisance_method,
                        working_model=working_model,
                        p_rct = pRCT,
                        verbose = FALSE)
    } else if (method == "atmle psi_tilde") {
      res <- psi_tilde_only_atmle(data,
                                  S_node = 1,
                                  W_node = c(2, 3, 4, 5),
                                  A_node = 6,
                                  Y_node = 7,
                                  nuisance_method=nuisance_method,
                                  working_model=working_model,
                                  p_rct = pRCT,
                                  verbose = FALSE)
    } else if (method == "tmle psi_tilde") {
      res <- psi_tilde_only_tmle(data,
                                 S_node = 1,
                                 W_node = c(2, 3, 4, 5),
                                 A_node = 6,
                                 Y_node = 7,
                                 nuisance_method=nuisance_method,
                                 working_model=working_model,
                                 p_rct = pRCT,
                                 verbose = FALSE)
    } else if (method == "escvtmle") {
      data$S <- 1 - data$S
      tmp <- ES.cvtmle(txinrwd = TRUE,
                       data = data,
                       study = "S",
                       covariates = c("W1", "W2", "W3", "W4"),
                       treatment_var = "A",
                       treatment = 1,
                       outcome = "Y",
                       pRCT = pRCT,
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
                      nuisance_method=nuisance_method,
                      p_rct = pRCT,
                      verbose = FALSE)
    } else if (method == "tmle") {
      res <- nonparametric(data,
                           S_node = 1,
                           W_node = c(2, 3, 4, 5),
                           A_node = 6,
                           Y_node = 7,
                           nuisance_method = nuisance_method,
                           working_model = working_model,
                           p_rct = pRCT,
                           verbose = FALSE)
    }

    if (res$lower <= bA & res$upper >= bA) {
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
                               pRCT,
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
      data <- generate_data(n, bA, bias, pRCT)

      # fit
      res <- NULL
      if (method == "atmle") {
        res <- atmle(data,
                     S_node = 1,
                     W_node = c(2, 3, 4, 5),
                     A_node = 6,
                     Y_node = 7,
                     nuisance_method=nuisance_method,
                     working_model=working_model,
                     p_rct = pRCT,
                     verbose = FALSE)
      } else if (method == "atmle_tmle") {
        res <- atmle_tmle(data,
                          S_node = 1,
                          W_node = c(2, 3, 4, 5),
                          A_node = 6,
                          Y_node = 7,
                          nuisance_method=nuisance_method,
                          working_model=working_model,
                          p_rct = pRCT,
                          verbose = FALSE)
      } else if (method == "atmle psi_tilde") {
        res <- psi_tilde_only_atmle(data,
                                    S_node = 1,
                                    W_node = c(2, 3, 4, 5),
                                    A_node = 6,
                                    Y_node = 7,
                                    nuisance_method=nuisance_method,
                                    working_model=working_model,
                                    p_rct = pRCT,
                                    verbose = FALSE)
      } else if (method == "tmle psi_tilde") {
        res <- psi_tilde_only_tmle(data,
                                   S_node = 1,
                                   W_node = c(2, 3, 4, 5),
                                   A_node = 6,
                                   Y_node = 7,
                                   nuisance_method=nuisance_method,
                                   working_model=working_model,
                                   p_rct = pRCT,
                                   verbose = FALSE)
      } else if (method == "escvtmle") {
        data$S <- 1 - data$S
        tmp <- ES.cvtmle(txinrwd = TRUE,
                         data = data,
                         study = "S",
                         covariates = c("W1", "W2", "W3", "W4"),
                         treatment_var = "A",
                         treatment = 1,
                         outcome = "Y",
                         pRCT = pRCT,
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
                             p_rct = pRCT,
                             verbose = FALSE)
      } else if (method == "rct_only") {
        res <- rct_only(data,
                        S_node = 1,
                        W_node = c(2, 3, 4, 5),
                        A_node = 6,
                        Y_node = 7,
                        nuisance_method = nuisance_method,
                        p_rct = pRCT,
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
