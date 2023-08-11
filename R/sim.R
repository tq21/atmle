
generate_data <- function(n_rct=200, n_rwd=1000, p_rct=0.67, bA=-0.6, B=0.21){
  # U
  UY_rct <- rnorm(n_rct, 0, 1.5)
  UY_rwd <- rnorm(n_rwd, 0, 1.5)

  # W
  W1_rct <- rnorm(n_rct, 0, 1)
  W1_rwd <- rnorm(n_rwd, 0, 1)
  W2_rct <- rnorm(n_rct, 0, 1)
  W2_rwd <- rnorm(n_rwd, 0, 1)

  # A
  A_rct <- rbinom(n_rct, 1, p_rct)
  A_rwd <- rbinom(n_rwd, 1, plogis(0.1 * W1_rwd - 0.5 * W2_rwd))

  # bias
  B1 <- rnorm(n_rwd, 3/4*B, 0.02)
  B2 <- rnorm(n_rwd, 1/4*B, 0.02)

  # Y
  Y_rct <- -3+2*W1_rct^2+W2_rct^3+W1_rct*W2_rct+bA*A_rct+UY_rct
  Y_rwd <- -3+2*W1_rwd^2+W2_rwd^3+W1_rwd*W2_rwd+bA*A_rwd+UY_rwd

  # data
  dt_rct <- data.frame(S=1, W1=W1_rct, W2=W2_rct, A=A_rct, Y=Y_rct)
  dt_rwd <- data.frame(S=0, W1=W1_rwd, W2=W2_rwd, A=A_rwd, Y=Y_rwd)
  data <- rbind(dt_rct, dt_rwd)

  return(data)
}

set.seed(21411)

data <- generate_data()

S <- data$S
W <- as.matrix(data[, c("W1", "W2")])
A <- data$A
Y <- data$Y
