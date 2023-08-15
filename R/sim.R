generate_data <- function(N, p_rct, bA){
  # U
  UY <- rnorm(N, 0, 1)

  # W
  W1 <- runif(N, -1, 1)
  W2 <- runif(N, -1, 1)
  W3 <- runif(N, -1, 1)
  W4 <- runif(N, -1, 1)

  # S
  S <- rbinom(N, 1, plogis(-3.8*as.numeric(W1 > 0)*as.numeric(W2 > 0)))

  # A
  A <- vector(length = N)
  A[S == 1] <- rbinom(sum(S), 1, p_rct)
  g_rwd <- 1.2*as.numeric(W1 < 0)+1.1*as.numeric(W1 > 0.4)+
    2*as.numeric(W3 > 0)-2.5*as.numeric(W4 > -0.5)

  # g_rwd <- as.numeric(W1 < 0)*W3+1.1*as.numeric(W1 > 0.4)-
  #   2*as.numeric(W1 > 0)+2.5*as.numeric(W1 > 0.5)*W3-2.5*as.numeric(W1 > 0.4) +
  #   as.numeric(W2 > -0.4)-4*as.numeric(W2 > 0)*W3+2*as.numeric(W2 > -0.3)
  A[S == 0] <- rbinom(N-sum(S), 1, plogis(g_rwd[S == 0]))

  # bias
  B <- rnorm(N, 0.2, 0.02)

  # Y
  # Y <- 1.1+bA*A+3.8*W1*as.numeric(W2 < 0)*sin(pi/2*abs(W1))+
  #   4*as.numeric(W2 > 0)*cos(pi/2*abs(W1))+4*as.numeric(W1 < 0)*cos(pi/2*abs(W3))+6*S+
  #   0.1*W3*sin(pi*W4)+W3*cos(abs(W4-W3))+UY
  Y <- 1.1+bA*A+
    3.8*as.numeric(W2 < 0)+4*as.numeric(W2 > 0)+4*as.numeric(W1 < 0)+
    1.5*as.numeric(W3 < -0.5)+1.2*as.numeric(W4 > 0.7)+B*(1-S)+UY

  # data
  data <- data.frame(S, W1, W2, W3, W4, A, Y)

  return(data)
}

# data <- generate_data(N=500, p_rct=0.67, bA=0.7)
