generate_data <- function(N, bA, B, rwd_rand){
  # U
  UY <- rnorm(N, 0, 1)

  # W
  W1 <- runif(N, -1, 1)
  W2 <- runif(N, -1, 1)
  W3 <- runif(N, -1, 1)
  W4 <- runif(N, -1, 1)

  # S
  S <- rbinom(N, 1, 0.2)

  # A
  A <- vector(length = N)
  A[S == 1] <- rbinom(sum(S), 1, 0.5)
  g_rwd <- rep(0, N)
  if (!rwd_rand) {
    g_rwd <- 1.2*as.numeric(W1 >= 0)+1.1*as.numeric(W2 >= 0.4)+
      2*as.numeric(W3 >= 0)-2.5*as.numeric(W4 >= -0.5)
  }
  A[S == 0] <- rbinom(N-sum(S), 1, plogis(g_rwd[S == 0]))

  # Y
  Y <- 1.1+bA*A+3.8*as.numeric(W2 >= 0)+4*as.numeric(W2 >= 0)+
    1.5*as.numeric(W3 >= -0.5)+1.2*as.numeric(W4 >= 0.7)+
    B+UY

  # data
  data <- data.frame(S, W1, W2, W3, W4, A, Y)

  return(data)
}

# RWD also randomized
data <- generate_data(N=500, bA=0.7, B=0, rwd_rand=TRUE)
res <- atmle(data,
             S_node = 1,
             W_node = c(2, 3, 4, 5),
             A_node = 6,
             Y_node = 7,
             target_Pi = FALSE,
             g_rct=0.5)


data
S_node = 1
W_node = c(2, 3, 4, 5)
A_node = 6
Y_node = 7
target_Pi = FALSE
g_rct=0.5
