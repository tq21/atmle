library(ggplot2)
library(devtools)
load_all()

`%+%` <- function(a, b) paste0(a, b)

generate_data <- function(N, bA, bias){
  # U
  UY <- rnorm(N, 0, 1)

  # S
  S <- rbinom(N, 1, 0.2)

  # W
  W1 <- rnorm(N)
  W2 <- rnorm(N)
  W3 <- rnorm(N)
  W4 <- rnorm(N)

  # A
  A <- vector(length = N)
  A[S == 0] <- rbinom(N - sum(S), 1, plogis(-0.5*W1+0.9*W2-1.2*W3+0.3*W4))
  #A[S == 0] <- rbinom(N - sum(S), 1, 0.5)
  A[S == 1] <- rbinom(sum(S), 1, 0.5)

  # Y
  Y <- vector(length = N)
  if (is.numeric(bias)) {
    # Y_S0 <- 0.3+bA*A+1.2*as.numeric(W1 > 0)+0.7*as.numeric(W2 < 0)+0.3*as.numeric(W3 > 0)-0.5*as.numeric(W4 < 0)+UY+bias # RWD has bias
    # Y_S1 <- 0.3+bA*A+1.2*as.numeric(W1 > 0)+0.7*as.numeric(W2 < 0)+0.3*as.numeric(W3 > 0)-0.5*as.numeric(W4 < 0)+UY
    Y_S0 <- 0.3+bA*A+0.5*W1+0.1*W2+0.3*W3-0.5*W4+bias+UY # RWD has bias
    Y_S1 <- 0.3+bA*A+0.5*W1+0.1*W2+0.3*W3-0.5*W4+UY
    Y[S == 0] <- Y_S0[S == 0]
    Y[S == 1] <- Y_S1[S == 1]
  } else if (bias == "parametric") {
    Y_S0 <- 0.3+bA*A+0.5*W1+0.1*W2+0.3*W3-0.5*W4+0.2*W1+0.5*W2-0.4*W4+UY
    Y_S1 <- 0.3+bA*A+0.5*W1+0.1*W2+0.3*W3-0.5*W4+UY
    Y[S == 0] <- Y_S0[S == 0]
    Y[S == 1] <- Y_S1[S == 1]
  } else if (bias == "A") {
    bias <- 0.5
    Y_S0 <- 0.3+bA*A+0.5*W1+0.1*W2+0.3*W3-0.5*W4+bias*A+UY # RWD has bias
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
                    verbose=TRUE,
                    atmle_both=TRUE) {
  # results
  psi_coverage <- rep(NA, B)
  psi_est <- rep(NA, B)
  psi_ci_lower <- rep(NA, B)
  psi_ci_upper <- rep(NA, B)

  for (i in 1:B) {
    if (verbose) print(i)

    # simulate data
    data <- generate_data(n, bA, bias)

    # fit
    res <- NULL
    if (atmle_both) {
      res <- atmle(data,
                   S_node = 1,
                   W_node = c(2, 3, 4, 5),
                   A_node = 6,
                   Y_node = 7,
                   nuisance_method=nuisance_method,
                   working_model=working_model,
                   verbose = FALSE)
    } else {
      res <- atmle_tmle(data,
                        S_node = 1,
                        W_node = c(2, 3, 4, 5),
                        A_node = 6,
                        Y_node = 7,
                        nuisance_method=nuisance_method,
                        working_model=working_model,
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
              psi_ci_upper = psi_ci_upper))
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
                               verbose=TRUE,
                               atmle_both=TRUE) {

  n_seq <- seq(n_min, n_max, n_step)

  # results
  all_psi_est <- vector(mode = "list", length = length(n_seq))
  all_psi_coverage <- vector(mode = "list", length = length(n_seq))
  all_psi_ci_lower <- vector(mode = "list", length = length(n_seq))
  all_psi_ci_upper <- vector(mode = "list", length = length(n_seq))

  for (i in 1:length(n_seq)) {
    n <- n_seq[i]
    print(n)

    psi_est <- vector(length = B)
    psi_coverage <- vector(length = B)
    psi_ci_lower <- vector(length = B)
    psi_ci_upper <- vector(length = B)

    for (j in 1:B) {
      # simulate data
      data <- generate_data(n, bA, bias)

      # fit
      res <- NULL
      if (atmle_both) {
        res <- atmle(data,
                     S_node = 1,
                     W_node = c(2, 3, 4, 5),
                     A_node = 6,
                     Y_node = 7,
                     nuisance_method=nuisance_method,
                     working_model=working_model,
                     verbose = FALSE)
      } else {
        res <- atmle_tmle(data,
                          S_node = 1,
                          W_node = c(2, 3, 4, 5),
                          A_node = 6,
                          Y_node = 7,
                          nuisance_method=nuisance_method,
                          working_model=working_model,
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
  }

  return(list(all_psi_est = all_psi_est,
              all_psi_coverage = all_psi_coverage,
              all_psi_ci_lower = all_psi_ci_lower,
              all_psi_ci_upper = all_psi_ci_upper))
}
