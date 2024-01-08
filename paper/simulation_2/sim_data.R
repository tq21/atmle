sim_two_covs <- function(ate, n, rct_prop, g_rct, bias, controls_only) {
  # error
  #UY <- round(rnorm(n, 0, 1), 1)
  UY <- round(rnorm(n, 0, 1), 1)
  Ubias <- round(rnorm(n, 0, 0.01), 1)

  # baseline covariates
  W1 <- round(runif(n, -0.4, 0.4), 1)
  W2 <- round(runif(n, -0.4, 0.4), 1)

  # study indicator, S=1 for RCT, S=0 for RWD
  S <- rbinom(n, 1, rct_prop)

  # treatments (external data has both treated and controls)
  A <- numeric(length = n)
  A[S == 1] <- rbinom(sum(S), 1, g_rct)
  if (controls_only) {
    A[S == 0] <- rep(0, n - sum(S))
  } else {
    A[S == 0] <- rbinom(n - sum(S), 1, plogis(0.8*W1-0.9*W2))
  }

  # bias term for RWD data
  b <- NULL
  if (is.numeric(bias)) {
    b <- bias
  } else if (bias == "HAL_1") {
    b <- 3.8*W1*W2
  } else if (bias == "HAL_2") {
    b <- 5.9*W1+7.2*W2+8.3*W1*W2+0.9*sin(2*pi*abs(W1))
  } else if (bias == "HAL_3") {
    b <- 9.9*W1+3.2*W2+9.3*W1*W2+1.21*as.numeric(W1 >= 0)*W2
  }

  # outcome
  Y <- -2.4-1.9*W1-2.5*W2+ate*A+UY+(1-S)*(b)+(1-S)*Ubias

  # data frames combining RCT and RWD
  data <- data.frame(S = S,
                     W1 = W1,
                     W2 = W2,
                     A = A,
                     Y = Y)

  return(data)
}
