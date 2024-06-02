sim_data <- function(n) {

  W1 <- round(runif(n, 0, 1), 2)
  W2 <- round(runif(n, 0, 1), 2)
  A <- rbinom(n, 1, prob = plogis(0.5 * W1 - 0.25 * W2))

  # failure time
  lambda_T <- exp(W1 + W2 + 0.5 * A)
  T <- rpois(n, lambda = lambda_T) + 1

  # censoring time
  lambda_C <- exp(W1 + 2 * W2)
  C <- rpois(n, lambda = lambda_C) + 1

  # compute T_tilde=min(T,C) and Delta=I(T<=C)
  T_tilde <- pmin(T, C)
  Delta <- as.numeric(T <= C)

  return(data.frame(W1 = W1,
                    W2 = W2,
                    A = A,
                    T_tilde = T_tilde,
                    Delta = Delta))
}

library(devtools)
load_all()
library(hal9001)


set.seed(123)
data <- sim_data(100)
S_node <- 1
W_node <- c(1, 2)
A_node <- 3
T_tilde_node <- 4
Delta_node <- 5
t0 <- 2
g_rct <- 0.5
controls_only <- FALSE
g_method <- "glm"
G_bar_method <- "glm"
lambda_method <- "HAL"
v_folds <- 5
g_bounds <- c(0.01, 0.99)
G_bar_bounds <- c(0.01, 0.99)


