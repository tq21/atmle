
generate_data <- function(N, pRCT, bA=-0.7){
  # U
  U_Y <- rnorm(N, 0, 1)

  # W
  W1 <- rnorm(N, 0, 1)
  W2 <- rnorm(N, 0, 1)

  # A
  A <- rbinom(N, 1, pRCT)

  # S
  S <- rbinom(N, 1, plogis(0.1 * W1^3 - 0.5 * W2^2 + 0.1 * A))

  # Y
  Y <- 0.1 * W1 - 0.5 * W2 + bA * A + 0.1 * W1^2 * S + U_Y

  data <- data.frame(S, W1, W2, A, Y)

  return(data)
}

set.seed(21411)

data <- generate_data(1000, 0.5)

S <- data$S
W <- as.matrix(data[, c("W1", "W2")])
A <- data$A
Y <- data$Y
