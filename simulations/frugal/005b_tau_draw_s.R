library(atmle)
library(EScvtmle)
library(devtools)
load_all()

`%+%` <- function(a, b) paste0(a, b)

# DGP: SCM VERSION
sim_data <- function(n, gamma, rct_prop, A_counter = -1) {
  W1 <- rnorm(n, 0, 1)
  W2 <- rbinom(n, 1, 0.5)
  W3 <- rnorm(n, 0.1*W2, 1)
  U1 <- rnorm(n, 0, 1)
  U2 <- rnorm(n, 0, 0.35)

  S <- rbinom(n, 1, rct_prop)

  A <- numeric(length = n)
  A[S == 1] <- rbinom(sum(S), 1, 0.5)
  A[S == 0] <- rbinom(n - sum(S), 1, plogis(W1+W2+gamma*U1+gamma*U2))

  if (A_counter == 1) {
    A <- rep(1, n)
  } else if (A_counter == 0) {
    A <- rep(0, n)
  }

  Y <- rnorm(n, (0.2+0.1*W3)*A+U1+U2, 1)

  data <- data.frame(S = S,
                     W1 = W1,
                     W2 = W2,
                     W3 = W3,
                     A = A,
                     Y = Y)

  return(data)
}

# SCM SIMULATIONS
set.seed(124)
B <- 100
n <- 1500
n_large <- 10^7
bias_seq <- 0.5#seq(0, 2, 0.1)
g_rct <- 0.5
rct_prop <- 0.2
controls_only <- FALSE
S_node <- 1
W_node <- 2:4
A_node <- 5
Y_node <- 6

true_tau_A1 <- numeric(length = length(bias_seq))
true_tau_A0 <- numeric(length = length(bias_seq))
atmle_tau_A1 <- matrix(nrow = B, ncol = length(bias_seq))
atmle_tau_A0 <- matrix(nrow = B, ncol = length(bias_seq))

# use monte-carlo simulation to evaluate the true tau
for (i in 1:length(bias_seq)) {
  # simulate large data
  rct_A1 <- sim_data(n_large, bias_seq[i], 1, A_counter = 1)
  rct_A0 <- sim_data(n_large, bias_seq[i], 1, A_counter = 0)
  rwd_A1 <- sim_data(n_large, bias_seq[i], 0, A_counter = 1)
  rwd_A0 <- sim_data(n_large, bias_seq[i], 0, A_counter = 0)

  # estimate true tau(1,W)
  true_tau_A1[i] <- mean(rct_A0$Y) - mean(rwd_A0$Y)
  true_tau_A0[i] <- mean(rct_A1$Y) - mean(rwd_A1$Y)
}

for (b in 1:B) {
  print("run: " %+% b)

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

    # collect results
    atmle_tau_A1[b, i] <- mean(atmle_glmnet_res$tau_A1)
    atmle_tau_A0[b, i] <- mean(atmle_glmnet_res$tau_A0)
  }
}

mean(atmle_tau_A1)-true_tau_A1
mean(atmle_tau_A0)-true_tau_A0

save(list = c("true_tau_A1",
              "true_tau_A0",
              "atmle_tau_A1",
              "atmle_tau_A0"), file = "out/tau_draw_s.RData")
