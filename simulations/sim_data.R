
sim_four_covs <- function(ate, n, rct_prop, g_rct, bias, controls_only, outcome_type) {
  # error
  UY <- rnorm(n, 0, 1.5)
  U_bias <- rnorm(n, 0, 0.02)

  # baseline covariates
  W1 <- rnorm(n, 0, 1)
  W2 <- rnorm(n, 0, 1)
  W3 <- rnorm(n, 0, 1)
  W4 <- rnorm(n, 0, 1)

  # study indicator, S=1 for RCT, S=0 for RWD
  S <- rbinom(n, 1, rct_prop)

  # treatments (external data has both treated and controls)
  A <- numeric(length = n)
  A[S == 1] <- rbinom(sum(S), 1, g_rct)
  if (controls_only) {
    A[S == 0] <- rep(0, n - sum(S))
  } else {
    A[S == 0] <- rbinom(n - sum(S), 1, plogis(0.9*W1+0.4*W2))
    #A[S == 0] <- rbinom(n - sum(S), 1, g_rct)
  }

  # bias term for RWD data
  b <- NULL
  if (is.numeric(bias)) {
    b <- bias
  } else if (bias == "param_simple") {
    b <- 0.8*W1
  } else if (bias == "param_complex") {
    b <- 0.6*W1+1.5*W2+1.4*W3+0.8*W4
  } else if (bias == "HAL") {
    b <- 2.5*as.numeric(W1 > 0)+1.9*W2+0.7*as.numeric(W1 < 0)*W2
  } else if (bias == "sinusoidal") {
    b <- 3.8*W2*as.numeric(W1 < 0.5)*sin(pi/2*abs(W1))+
      4*as.numeric(W2 > 0.7)*cos(pi/2*abs(W2))
  }

  # outcome
  Y <- NULL
  if (outcome_type == "binomial") {
    Y <- rbinom(n, 1, plogis(0.5+0.8*W1+1.1*W2+0.9*W3+1.3*W4+ate*A+UY+(1-S)*(1-A)*b))
  } else {
    Y <- round(0.5+0.8*W1+1.1*W2+0.9*W3+1.3*W4+ate*A+UY+(1-S)*(1-A)*(b+U_bias), 2)
  }

  # data frames combining RCT and RWD
  data <- data.frame(S = S,
                     W1 = W1,
                     W2 = W2,
                     W3 = W3,
                     W4 = W4,
                     A = A,
                     Y = Y)

  return(data)
}
