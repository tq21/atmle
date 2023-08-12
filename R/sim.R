

generate_data <- function(n_rct=200, n_rwd=1000, p_rct=0.67, bA=-3.2, gamma=0.5, B=0.21){
  N <- n_rct + n_rwd

  # U
  UY <- rnorm(N, 0, 1)

  # S
  S <- c(rep(1, n_rct), rep(0, n_rwd))

  # W
  W1 <- runif(N, -1, 1)
  W2 <- runif(N, -1, 1)
  W3 <- runif(N, -1, 1)
  W4 <- runif(N, -1, 1)

  # A
  A_rct <- rbinom(n_rct, 1, p_rct)
  g_rwd <- as.numeric(W1 < -3)*W3+1.1*as.numeric(W1 > -2)-
    2*as.numeric(W1 > 0)+2.5*as.numeric(W1 > 2)*W3-2.5*as.numeric(W1 > 3) +
    as.numeric(W2 > -1)-4*as.numeric(W2 > 1)*W3+2*as.numeric(W2 > 3)
  A_rwd <- rbinom(n_rwd, 1, plogis(g_rwd))
  A <- c(A_rct, A_rwd)

  # Y
  Y <- 1.1+bA*A+3.8*W1*as.numeric(W2 < 0)*sin(pi/2*abs(W1))+
    4*as.numeric(W2 > 0)*cos(pi/2*abs(W1))+4*as.numeric(W1 < 0)*cos(pi/2*abs(W3))+3*S+
    0.1*W3*sin(pi*W4)+W3*cos(abs(W4-W3))+UY

  # data
  data <- data.frame(S, W1, W2, W3, W4, A, Y)

  return(data)
}

set.seed(21411)
data <- generate_data(n_rct=200, n_rwd=500, p_rct=0.67, bA=5.9, gamma=0.5)

S <- data$S
W <- as.matrix(data[, c("W1", "W2", "W3", "W4")])
A <- data$A
Y <- data$Y
