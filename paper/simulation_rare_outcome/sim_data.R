sim_data <- function(beta,
                     n,
                     rct,
                     g_rct,
                     bias,
                     A_counter = -1) {

  # baseline covariates
  W1 <- round(runif(n, 0, 1), 1)
  W2 <- round(runif(n, 0, 1), 1)
  W3 <- rbinom(n, 1, 0.5)

  # study indicator S and treatment A
  if (rct) {
    S <- rep(1, n)
    if (A_counter == 1) {
      A <- rep(1, n)
    } else if (A_counter == 0) {
      A <- rep(0, n)
    } else {
      A <- rbinom(n, 1, g_rct)
    }

  } else {
    S <- rep(0, n)
    A <- rep(0, n)
  }

  # bias term for RWD data
  if (bias == "a") {
    b <- 0.2+0.1*W1*(1-A)
  } else if (bias == "b") {
    b <- 0.5+1.1*W1*(1-A)+0.8*W3
  }

  # outcome
  Y <- rbinom(n = n, size = 1,
              prob = plogis(-3-0.9*W1+0.5*W2-0.7*W3+beta*A+(1-S)*b))

  # data frames combining RCT and RWD
  data <- data.frame(S = S,
                     W1 = W1,
                     W2 = W2,
                     W3 = W3,
                     A = A,
                     Y = Y)

  return(data)
}
