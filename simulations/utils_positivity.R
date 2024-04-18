generate_data <- function(N, bA, bias, pRCT){
  # U
  UY <- rnorm(N, 0, 1)

  # S
  S <- rbinom(N, 1, 0)

  # W
  W1 <- runif(N)
  W2 <- runif(N)
  W3 <- runif(N)
  W4 <- runif(N)

  # A
  A <- vector(length = N)
  #A[S == 0] <- rbinom(N - sum(S), 1, 0.01)
  #A[S == 0] <- rbinom(N - sum(S), 1, plogis(-2*W1-2*W2+0.9*W2-1.2*W3+0.3*W4))
  A[S == 0] <- rbinom(N - sum(S), 1, plogis(0.1*W1-0.2*W2+1.2*W3-0.3*W4))
  A[S == 1] <- rbinom(sum(S), 1, pRCT) # rct

  # Y
  Y <- vector(length = N)
  if (is.numeric(bias)) {
    Y_S0 <- 0.3+bA*A+0.5*W1+0.1*W2+0.3*W3-0.5*W4+bias+UY # RWD has bias
    Y_S1 <- 0.3+bA*A+0.5*W1+0.1*W2+0.3*W3-0.5*W4+UY
    Y[S == 0] <- Y_S0[S == 0]
    Y[S == 1] <- Y_S1[S == 1]
  } else if (bias == "param_simple") {
    Y_S0 <- 0.3+bA*A+0.5*W1+0.1*W2+0.3*W3-0.5*W4+1.3*W1-0.6+UY
    Y_S1 <- 0.3+bA*A+0.5*W1+0.1*W2+0.3*W3-0.5*W4+UY
    Y[S == 0] <- Y_S0[S == 0]
    Y[S == 1] <- Y_S1[S == 1]
  } else if (bias == "param_complex") {
    Y_S0 <- 0.3+bA*A+0.5*W1+0.1*W2+0.3*W3-0.5*W4+1.3*W1+1.1*W2+0.9*W3-2.4*W4+0.6+UY
    Y_S1 <- 0.3+bA*A+0.5*W1+0.1*W2+0.3*W3-0.5*W4+UY
    Y[S == 0] <- Y_S0[S == 0]
    Y[S == 1] <- Y_S1[S == 1]
  } else if (bias == "HAL") {
    Y_S0 <- 0.3+bA*A+0.5*W1+0.1*W2+0.3*W3-0.5*W4+UY+
      -0.5
    Y_S1 <- 0.3+bA*A+0.5*W1+0.1*W2+0.3*W3-0.5*W4+UY
    Y[S == 0] <- Y_S0[S == 0]
    Y[S == 1] <- Y_S1[S == 1]
  } else if (bias == "hard") {
    Y_S0 <- 0.3+bA*A+0.5*W1+0.1*W2+0.3*W3-0.5*W4+UY
    1.3*W1+0.9*W2^2+0.7*W3^3+1.5*W4^2+1.5*W1*W2+1.9*W1*W2*W4^2-0.6
    Y_S1 <- 0.3+bA*A+0.5*W1+0.1*W2+0.3*W3-0.5*W4+UY
    Y[S == 0] <- Y_S0[S == 0]
    Y[S == 1] <- Y_S1[S == 1]
  } else if (bias == "sinusoidal") {
    Y_S0 <- 0.3+bA*A+0.5*W1+0.1*W2+0.3*W3-0.5*W4+UY+
      3.8*W3*as.numeric(W1 < 0.5)* sin(pi/2*abs(W1))+
      4*as.numeric(W2 > 0.7)*cos(pi/2*abs(W1))+
      1.9*W4*sin(pi*W4)+3.2*W3*cos(abs(W4-W3))
    Y_S1 <- 0.3+bA*A+0.5*W1+0.1*W2+0.3*W3-0.5*W4+UY
    Y[S == 0] <- Y_S0[S == 0]
    Y[S == 1] <- Y_S1[S == 1]
  } else if (bias == "test 2") {
    Y_S0 <- 0.3+bA*A+0.5*W1+0.1*W2+0.3*W3-0.5*W4+UY+100
    Y_S1 <- 0.3+bA*A+0.5*W1+0.1*W2+0.3*W3-0.5*W4+UY
    Y[S == 0] <- Y_S0[S == 0]
    Y[S == 1] <- Y_S1[S == 1]
  }

  # data
  data <- data.frame(S, W1, W2, W3, W4, A, Y)

  return(data)
}
