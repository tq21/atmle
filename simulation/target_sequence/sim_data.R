sim_data <- function(ate,
                     n,
                     rct,
                     g_rct,
                     bias,
                     controls_only) {
  # error
  UY <- rnorm(n, 0, 1)

  # baseline covariates
  W1 <- round(runif(n, 0, 1), 1)
  W2 <- round(runif(n, 0, 1), 1)
  W3 <- round(runif(n, 0, 1), 1)

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
  } else if (bias == "c") {
    b <- 0.3+0.9*W2*(1-A)+0.7*as.numeric(W2 > 0.5)*W3
  } else if (bias == "d") {
    b <- 0.3+1.1*W1*(1-A)+0.9*W2^2*W3
  }

  # outcome
  Y <- 1.9+ate*A+0.9*W1+1.4*W2+2.1*W3+(1-S)*b+UY

  # data frames combining RCT and RWD
  data <- data.frame(S = S,
                     W1 = W1,
                     W2 = W2,
                     W3 = W3,
                     A = A,
                     Y = Y)

  return(data)
}
