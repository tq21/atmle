library(ggplot2)
library(devtools)
library(EScvtmle)
load_all()

`%+%` <- function(a, b) paste0(a, b)

generate_data <- function(N, bA, bias, pRCT){
  # U
  UY <- rnorm(N, 0, 1)

  # S
  S <- rbinom(N, 1, 0.2)
  #S <- sample(c(rep(1, 200), rep(0, N - 200)))

  # W
  # W1 <- rnorm(N)
  # W2 <- rnorm(N)
  # W3 <- rnorm(N)
  # W4 <- rnorm(N)
  W1 <- runif(N)
  W2 <- runif(N)
  W3 <- runif(N)
  W4 <- runif(N)

  # A
  A <- vector(length = N)
  #A[S == 0] <- rbinom(N - sum(S), 1, 0.6*W1)
  #A[S == 0] <- rbinom(N - sum(S), 1, plogis(-0.5*W1+0.9*W2-1.2*W3+0.3*W4))
  A[S == 0] <- rbinom(N - sum(S), 1, 0.67)
  #A[S == 0] <- rbinom(N - sum(S), 1, plogis(-2*W1-2*W2+1.2*W3-0.3*W4))
  A[S == 1] <- rbinom(sum(S), 1, pRCT)

  # Y
  Y <- vector(length = N)
  if (is.numeric(bias)) {
    #b <- rnorm(N, bias, 0.1)
    b <- bias
    Y_S0 <- 0.3+bA*A+0.5*W1+0.1*W2+0.3*W3-0.5*W4+UY-b*A
    Y_S1 <- 0.3+bA*A+0.5*W1+0.1*W2+0.3*W3-0.5*W4+UY
    Y[S == 0] <- Y_S0[S == 0]
    Y[S == 1] <- Y_S1[S == 1]
  } else if (bias == "param_simple") {
    b <- 2.7*W1+0.1+rnorm(N, 0.01)
    Y_S0 <- 0.3+bA*A+0.5*W1+0.1*W2+0.3*W3-0.5*W4+UY-b*A
    Y_S1 <- 0.3+bA*A+0.5*W1+0.1*W2+0.3*W3-0.5*W4+UY
    Y[S == 0] <- Y_S0[S == 0]
    Y[S == 1] <- Y_S1[S == 1]
  } else if (bias == "param_complex") {
    b <- -0.1+1.5*W1+1.2*W2+2.5*W3+0.6*W4
    Y_S0 <- 0.3+bA*A+0.5*W1+0.1*W2+0.3*W3-0.5*W4+UY-b*A
    Y_S1 <- 0.3+bA*A+0.5*W1+0.1*W2+0.3*W3-0.5*W4+UY
    Y[S == 0] <- Y_S0[S == 0]
    Y[S == 1] <- Y_S1[S == 1]
  } else if (bias == "HAL") {
    Y_S0 <- 0.3+bA*A+0.5*W1+0.1*W2+0.3*W3-0.5*W4+UY+
      100*as.numeric(W1 > 0.5)+5*as.numeric(W2 > 0.5)+5*as.numeric(W3 > 0.3)+5*as.numeric(W4 > 0.7)-0.5
    Y_S1 <- 0.3+bA*A+0.5*W1+0.1*W2+0.3*W3-0.5*W4+UY
    Y[S == 0] <- Y_S0[S == 0]
    Y[S == 1] <- Y_S1[S == 1]
  } else if (bias == "hard") {
    Y_S0 <- 0.3+bA*A+0.5*W1+0.1*W2+0.3*W3-0.5*W4+UY
    1.3*W1+0.9*W2+0.7*W3^3+1.5*W4^2+1.5*W1*W2+0.9*W1*W4^2-0.6
    Y_S1 <- 0.3+bA*A+0.5*W1+0.1*W2+0.3*W3-0.5*W4+UY
    Y[S == 0] <- Y_S0[S == 0]
    Y[S == 1] <- Y_S1[S == 1]
  } else if (bias == "test") {
    Y_S0 <- 0.3+bA*A+0.5*W1+0.1*W2+0.3*W3-0.5*W4+1.2*W1+0.8*W2+0.7*W3-0.4*W4+0.8+UY
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

# generate_data <- function(N, bA, bias, pRCT, prop_rct=0.5){
#   n_rct <- round(N * prop_rct)
#   n_rwd <- N - n_rct
#
#   # RCT data
#   # U
#   UY_rct <- rnorm(n_rct, 0, 1)
#
#   # W
#   W1_rct <- runif(n_rct)
#   W2_rct <- runif(n_rct)
#   W3_rct <- runif(n_rct)
#   W4_rct <- runif(n_rct)
#
#   # A
#   A_rct <- rbinom(n_rct, 1, pRCT)
#
#   # RWD data
#   # U
#   UY_rwd <- rnorm(n_rwd, 0, 1)
#
#   # W
#   W1_rwd <- runif(n_rwd)
#   W2_rwd <- runif(n_rwd)
#   W3_rwd <- runif(n_rwd)
#   W4_rwd <- runif(n_rwd)
#
#   # A
#   A_rwd <- rbinom(n_rwd, 1, pRCT)
#
#   # combined
#   UY <- c(UY_rct, UY_rwd)
#   W1 <- c(W1_rct, W1_rwd)
#   W2 <- c(W2_rct, W2_rwd)
#   W3 <- c(W3_rct, W3_rwd)
#   W4 <- c(W4_rct, W4_rwd)
#   A <- c(A_rct, A_rwd)
#   S <- c(rep(1, n_rct), rep(0, n_rwd))
#
#   # Y
#   Y <- vector(length = N)
#   if (is.numeric(bias)) {
#     b <- bias
#     if (bias != 0) {
#       b <- rnorm(N, bias, 0.1)
#     } else {
#       b <- 0
#     }
#     Y_S0 <- 0.3+0.5*W1+0.1*W2+0.3*W3-0.5*W4+bA*A-b+UY
#     Y_S1 <- 0.3+0.5*W1+0.1*W2+0.3*W3-0.5*W4+bA*A+UY
#     Y[S == 0] <- Y_S0[S == 0]
#     Y[S == 1] <- Y_S1[S == 1]
#   } else if (bias == "param_simple") {
#     b <- 2.7*W1+0.1+rnorm(N, 0.05)
#     Y_S0 <- 0.3+bA*A+0.5*W1+0.1*W2+0.3*W3-0.5*W4-b+UY
#     Y_S1 <- 0.3+bA*A+0.5*W1+0.1*W2+0.3*W3-0.5*W4+UY
#     Y[S == 0] <- Y_S0[S == 0]
#     Y[S == 1] <- Y_S1[S == 1]
#   } else if (bias == "param_complex") {
#     b <- -0.1+1.5*W1+1.2*W2+2.5*W3+0.6*W4+rnorm(N, 0.05)
#     Y_S0 <- 0.3+0.5*W1+0.1*W2+0.3*W3-0.5*W4+A*bA-b+UY
#     Y_S1 <- 0.3+0.5*W1+0.1*W2+0.3*W3-0.5*W4+A*bA+UY
#     Y[S == 0] <- Y_S0[S == 0]
#     Y[S == 1] <- Y_S1[S == 1]
#   } else if (bias == "HAL") {
#     b <- 1.2+0.7*W1+2.1*W2+0.3*W3+1.5*W4+
#       1.5*as.numeric(W1 > 0.5)-0.7*as.numeric(W2 < 0.7)+1.8*as.numeric(W4 > 0.1)+rnorm(N, 0.01)
#     Y_S0 <- 0.3+0.5*W1+0.1*W2+0.3*W3-0.5*W4+A*bA-b*A+UY
#     Y_S1 <- 0.3+0.5*W1+0.1*W2+0.3*W3-0.5*W4+A*bA+UY
#     Y[S == 0] <- Y_S0[S == 0]
#     Y[S == 1] <- Y_S1[S == 1]
#   } else if (bias == "sinusoidal") {
#     b <- 0.3+0.5*W1+0.7*W2^2+0.1*W3+
#       1.2*W3*sin(pi/2*abs(W1))+
#       0.6*as.numeric(W2 > 0.7)*cos(pi/2*abs(W1))+
#       2.3*W3*cos(abs(W4-W3))+rnorm(N, 0.01)
#     Y_S0 <- 0.3+bA*A+0.5*W1+0.1*W2+0.3*W3-0.5*W4-b*A+UY
#     Y_S1 <- 0.3+bA*A+0.5*W1+0.1*W2+0.3*W3-0.5*W4+UY
#     Y[S == 0] <- Y_S0[S == 0]
#     Y[S == 1] <- Y_S1[S == 1]
#   } else if (bias == "test") {
#     b <- -1.2+1.5*W1^2+0.3*W1*W2-1.8*W3^3+0.5*W4+rnorm(N, 0.1)
#     Y_S0 <- 0.3+0.5*W1+0.1*W2+0.3*W3-0.5*W4+A*bA-b*A+UY
#     Y_S1 <- 0.3+0.5*W1+0.1*W2+0.3*W3-0.5*W4+A*bA+UY
#     Y[S == 0] <- Y_S0[S == 0]
#     Y[S == 1] <- Y_S1[S == 1]
#   }
#
#   # data
#   data <- data.frame(S, W1, W2, W3, W4, A, Y)
#
#   #A <- rbinom(N, 1, pRCT)
#   #A[S == 0] <- rbinom(N - sum(S), 1, 0.67)
#   #A[S == 0] <- rbinom(N - sum(S), 1, plogis(-0.5*W1+0.9*W2-1.2*W3+0.3*W4))
#   #A[S == 0] <- rbinom(N - sum(S), 1, plogis(-5*W1^2-5*W2^3+0.9*W2-1.2*W3+0.3*W4))
#   #A[S == 0] <- rbinom(N - sum(S), 1, 0.05)
#   #A[S == 1] <- rbinom(sum(S), 1, pRCT) # rct
#
#   return(data)
# }

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
      #data$S <- 1 - data$S
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
        #data$S <- 1 - data$S
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
