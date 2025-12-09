sim_data <- function(n,
                     prop_rct = 0.3, # this is ~ prop. of RCT patients
                     g_rct = 0.5, # prob. of receiving trt in RCT
                     bias = c("a", "b", "none"),
                     S_counter = NULL,
                     A_counter = NULL) {

  # baseline covariates
  W1 <- rnorm(n, 0, 1)
  W2 <- rnorm(n, 0, 1)
  W3 <- rnorm(n, 0, 1)

  # trial enrollment
  if (is.null(S_counter)) {
    logit_pS <- qlogis(prop_rct)+0.4*W1-0.2*W2-0.2*W3
    S <- rbinom(n, 1, plogis(logit_pS))
  } else {
    S <- rep(S_counter, n)
  }

  # treatment assignment
  if (is.null(A_counter)) {
    A <- numeric(n)
    A[S == 1] <- rbinom(sum(S == 1), 1, g_rct)
    A[S == 0] <- rbinom(sum(S == 0), 1, plogis(0.5*W1[S == 0]))
  } else {
    A <- rep(A_counter, n)
  }

  # outcome
  UY <- rnorm(n, 0, 1)
  tau_W <- 1.5+0.6*W1
  bias <- match.arg(bias)
  if (bias == "a") {
    b <- 0.2+0.1*W1*(1-A)
  } else if (bias == "b") {
    b <- 0.5+0.1*W1*(1-A)+0.8*W3
  } else {
    b <- 0
  }
  Y <- 2.5+0.9*W1+1.1*W2+2.7*W3+tau_W*A+UY+(1-S)*b

  return(data.frame(S, W1, W2, W3, A, Y))
}

#' @title Function to get true estimand value
#'
#' @param B Number of Monte-Carlo draws to approximate the true value. Default
#' is 10 million.
#'
#' @return A `numeric` vector of length 2. Index 1 should be the true ATE
#' averaged over pooled population, and index 2 should be the true ATE averaged
#' over the RCT-only population. For more details on the distinction of those
#' two estimands, see xxx.
get_truth_vdl25 <- function(B = 1e7) {
  set.seed(123)
  data_A1 <- sim_data_vdl25(n = B, bias = "none", A_counter = 1)
  data_A0 <- sim_data_vdl25(n = B, bias = "none", A_counter = 0)
  rct_PW <- mean(data_A1$Y[data_A1$S == 1]) - mean(data_A0$Y[data_A0$S == 1])
  return(c(1.5, rct_PW))
}
