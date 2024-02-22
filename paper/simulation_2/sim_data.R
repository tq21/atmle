sim_data <- function(ate,
                     n,
                     rct,
                     g_rct,
                     bias,
                     controls_only) {
  # error
  UY <- rnorm(n, 0, 1)

  # baseline covariates
  # W1 <- round(runif(n, -1, 1), 1)
  # W2 <- round(runif(n, -1, 1), 1)
  W1 <- round(rnorm(n, 0, 1), 2)
  W2 <- rbinom(n, 1, 0.5)

  # study indicator S and treatment A
  if (rct) {
    S <- rep(1, n)
    A <- rbinom(n, 1, g_rct)
  } else {
    S <- rep(0, n)
    if (controls_only) {
      A <- rep(0, n)
    } else {
      A <- rbinom(n, 1, plogis(W1))
    }
  }

  # bias term for RWD data
  if (is.numeric(bias)) {
    b <- bias
  } else if (bias == "HAL_1") {
    b <- 0.3+0.1*W1*W2*(1-A)+0.6*(1-A)
  } else if (bias == "HAL_2") {
    b <- 0.5+0.4*sin(2*pi*abs(W1))*(1-A)+0.3*(1-A)
  } else if (bias == "HAL_3") {
    b <- 0.8+(0.5*as.numeric(W1 >= 0)+0.3*W1*W2)*(1-A)
  }

  # outcome
  Y <- 2.5+0.9*W1+1.1*W2+ate*A+UY+(1-S)*b

  # data frames combining RCT and RWD
  data <- data.frame(S = S,
                     W1 = W1,
                     W2 = W2,
                     A = A,
                     Y = Y)

  return(data)
}
