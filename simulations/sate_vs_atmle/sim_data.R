sim_four_covs <- function(ate, n, rct_prop, g_rct, bias, controls_only) {
  # error
  UY <- rnorm(n, 0, 1)
  U_bias <- rnorm(n, 0, 1)

  # baseline covariates
  W1 <- rnorm(n, 1, 1)# runif(n, 0, 1) #
  W2 <- rnorm(n, 1, 1)# runif(n, 0, 1) #
  W3 <- rnorm(n, 1, 1)# runif(n, 0, 1) #
  W4 <- rnorm(n, 1, 1)# runif(n, 0, 1) #

  # study indicator, S=1 for RCT, S=0 for RWD
  S <- rbinom(n, 1, rct_prop)

  # treatments (external data has both treated and controls)
  A <- numeric(length = n)
  A[S == 1] <- rbinom(sum(S), 1, g_rct)
  if (controls_only) {
    A[S == 0] <- rep(0, n - sum(S))
  } else {
    A[S == 0] <- rbinom(n - sum(S), 1, plogis(-0.1*W1+0.2*W2))
    #A[S == 0] <- rbinom(n - sum(S), 1, 0.5)
  }

  # bias term for RWD data
  b <- NULL
  if (is.numeric(bias)) {
    b <- bias
  } else if (bias == "param_simple") {
    b <- 3.2+2.3*W1
  } else if (bias == "param_complex") {
    b <- 10.5+3.5*W1+5.4*W2+3.4*W3+6*W4
  } else if (bias == "HAL") {
    b <- 0.3+2.3*W1+1.5*W1*W2+3.4*W1^2*W2
  } else if (bias == "sinusoidal") {
    b <- 3.8*W2*as.numeric(W1 < 0.5)*sin(pi/2*abs(W1))+
      4*as.numeric(W2 > 0.7)*cos(pi/2*abs(W2))
  }

  # outcome
  Y <- 0.5+0.8*W1+1.1*W2+0.9*W3+1.3*W4+ate*A+UY+(1-S)*b+(1-S)*U_bias
  #Y <- rbinom(n, 1, plogis(-0.5*W1-0.2*W2-0.9*W3-0.5*W4+ate*A+(1-S)*(b)))

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
