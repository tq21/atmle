sim_data <- function(ate,
                     n,
                     rct,
                     g_rct,
                     bias,
                     controls_only) {
  # error
  UY <- rnorm(n, 0, 1)
  Ubias <- rnorm(n, 0, 0.2)

  # baseline covariates
  W1 <- rnorm(n, 0, 1)
  W2 <- rnorm(n, 0, 1)
  W3 <- rnorm(n, 0, 1)
  W4 <- rnorm(n, 0, 1)

  # study indicator S and treatment A
  if (rct) {
    S <- rep(1, n)
    A <- rbinom(n, 1, g_rct)
  } else {
    S <- rep(0, n)
    if (controls_only) {
      A <- rep(0, n)
    } else {
      A <- rbinom(n, 1, plogis(0.5*W1+0.5*W2))
    }
  }

  # bias term for RWD data
  b <- NULL
  if (is.numeric(bias)) {
    b <- bias*(1-A)
  } else if (bias == "small_bias") {
    b <- 0.3+1.2*W1*(1-A)
  } else if (bias == "large_bias") {
    b <- 10.2+5.3*W1*(1-A)+2.7*W2
  }

  # outcome
  Y <- -2.4-1.9*W1-2.5*W2+1.2*W3+2.7*W4+ate*A+UY+(1-S)*b

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
