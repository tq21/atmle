library(atmle)
library(EScvtmle)
library(devtools)
load_all()

`%+%` <- function(a, b) paste0(a, b)

# DGP: SCM VERSION
sim_data <- function(n, gamma, rct) {
  W1 <- rnorm(n, 0, 1)
  W2 <- rbinom(n, 1, 0.5)
  W3 <- rnorm(n, 0.1*W2, 1)
  U1 <- rnorm(n, 0, 1)
  U2 <- rnorm(n, 0, 0.35)
  A <- numeric(length = n)

  if (rct) {
    S <- rep(1, n)
    A <- rbinom(n, 1, 0.5)
  } else {
    S <- rep(0, n)
    A <- rbinom(n - sum(S), 1, plogis(W1+W2+gamma*U1+gamma*U2))
  }

  #Y <- rnorm(n, (0.2+0.1*W3)*A+U1+U2, 1)
  Y <- rnorm(n, 0.2*A+0.1*W3+U1+U2, 1)

  data <- data.frame(S = S,
                     W1 = W1,
                     W2 = W2,
                     W3 = W3,
                     A = A,
                     Y = Y)

  return(data)
}

# SCM SIMULATIONS
bias_seq <- c(0,0.5)#seq(0, 5, 0.1)
g_rct <- 0.5
n_rct <- 5000
n_rwd <- 5000
controls_only <- FALSE
S_node <- 1
W_node <- 2:4
A_node <- 5
Y_node <- 6
atmle_glmnet <- list()

# simulate data
rct_data_list <- lapply(bias_seq, function(bias) {
  return(sim_data(n_rct, bias, TRUE))
})
rwd_data_list <- lapply(bias_seq, function(bias) {
  return(sim_data(n_rwd, bias, FALSE))
})

for (i in 1:length(bias_seq)) {
  data <- rbind(rct_data_list[[i]], rwd_data_list[[i]])

  # atmle with glmnet
  atmle_glmnet_res <- atmle(data = data,
                            S_node = S_node,
                            W_node = W_node,
                            A_node = A_node,
                            Y_node = Y_node,
                            atmle_pooled = TRUE,
                            controls_only = FALSE,
                            theta_method = "glm",
                            Pi_method = "glm",
                            g_method = "glm",
                            theta_tilde_method = "glm",
                            Q_method = "glm",
                            bias_working_model = "glmnet",
                            pooled_working_model = "glmnet",
                            g_rct = g_rct,
                            family = "gaussian",
                            verbose = FALSE)

  atmle_glmnet[[i]] <- atmle_glmnet_res
}
