sim_data <- function(n, counter_A = NULL) {
  # error
  UY <- rnorm(n, 0, 1)

  # baseline covariates
  W1 <- round(runif(n, 0, 1), 1)
  W2 <- round(runif(n, 0, 1), 1)
  W3 <- round(runif(n, 0, 1), 1)

  # treatment
  if (is.null(counter_A)) {
    A <- rbinom(n, 1, plogis(-9*W1+9*W2))
  } else {
    A <- rep(counter_A, n)
  }

  # outcome
  Y <- 1.9+1.5*A+2.5*W1*A+0.7*W2*A+0.9*W1+1.4*W2+2.1*W3+UY

  data <- data.frame(W1 = W1,
                     W2 = W2,
                     W3 = W3,
                     A = A,
                     Y = Y)

  return(data)
}

get_truth <- function() {
  data_A1 <- sim_data(1e7, counter_A = 1)
  data_A0 <- sim_data(1e7, counter_A = 0)
  return(mean(data_A1$Y- data_A0$Y))
}
