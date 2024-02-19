#library(atmle)
library(EScvtmle)
library(dplyr)
devtools::load_all()

`%+%` <- function(x, y) paste0(x, y)

sim_data <- function(n, gamma, rct, A_counter = -1) {
  W1 <- rnorm(n, 0, 1)
  W2 <- rbinom(n, 1, 0.5)
  W3 <- 0.1*W2+rnorm(n, 0, 1)
  U1 <- rnorm(n, 0, 1)
  U2 <- rnorm(n, 0, 0.35)
  A <- numeric(length = n)

  if (rct) {
    S <- rep(1, n)
    A <- rbinom(n, 1, 0.5)
  } else {
    S <- rep(0, n)
    A <- rbinom(n, 1, plogis(W1+W2+gamma*U1+gamma*U2))
  }

  if (A_counter == 1) {
    A <- rep(1, n)
  } else if (A_counter == 0) {
    A <- rep(0, n)
  }

  Y <- (0.2+0.1*W3)*A+U1+U2+rnorm(n, 0, 1)

  data <- data.frame(S = S,
                     W1 = W1,
                     W2 = W2,
                     W3 = W3,
                     A = A,
                     Y = Y)

  return(data)
}

sim_data_v2 <- function(n, rct, g_rct, controls_only, A_counter = -1) {
  W1 <- rnorm(n, 0, 1)
  W2 <- rnorm(n, 0, 1)
  W3 <- rnorm(n, 0, 1)
  A <- numeric(length = n)

  if (rct) {
    S <- rep(1, n)
    A <- rbinom(n, 1, g_rct)
  } else {
    S <- rep(0, n)
    if (controls_only) {
      A <- rep(0, n)
    } else {
      A <- rbinom(n, 1, plogis(W1+W2))
    }
  }

  if (A_counter == 1) {
    A <- rep(1, n)
  } else if (A_counter == 0) {
    A <- rep(0, n)
  }

  # no bias
  # Y <- 1.9+1.5*A+4.2*W1+3.9*W2+rnorm(n, 0, 1)

  # constant bias
  # Y <- 1.9+1.5*A+4.2*W1+3.9*W2+(1-S)*(6)+rnorm(n, 0, 1)

  # parametric bias
  Y <- 1.9+1.5*A+4.2*W1+3.9*W2+(1-S)*(2+3*W1)+rnorm(n, 0, 1)

  data <- data.frame(S = S,
                     W1 = W1,
                     W2 = W2,
                     W3 = W3,
                     A = A,
                     Y = Y)

  return(data)
}

# n_large <- 10^6
# rct_large_A1 <- sim_data(n_large, 1, rct = TRUE, A_counter = 1)
# rct_large_A0 <- sim_data(n_large, 1, rct = TRUE, A_counter = 0)
# rwd_large_A1 <- sim_data(n_large, 1, rct = FALSE, A_counter = 1)
# rwd_large_A0 <- sim_data(n_large, 1, rct = FALSE, A_counter = 0)
# mean(rct_large_A1$Y)-mean(rwd_large_A1$Y)
# mean(rct_large_A0$Y)-mean(rwd_large_A0$Y)
#
# rct_large_A1_v2 <- sim_data_v2(n_large, rct = TRUE, A_counter = 1)
# rct_large_A0_v2 <- sim_data_v2(n_large, rct = TRUE, A_counter = 0)
# rwd_large_A1_v2 <- sim_data_v2(n_large, rct = FALSE, A_counter = 1)
# rwd_large_A0_v2 <- sim_data_v2(n_large, rct = FALSE, A_counter = 0)
# mean(rct_large_A1_v2$Y)-mean(rwd_large_A1_v2$Y)
# mean(rct_large_A0_v2$Y)-mean(rwd_large_A0_v2$Y)

# SIMULATION PARAMETERS
controls_only <- TRUE
B <- 500
ate <- 1.5
g_rct <- 0.67
atmle_cover <- numeric(B)
atmle_est <- numeric(B)
escvtmle_cover <- numeric(B)
escvtmle_est <- numeric(B)

for (i in 1:B) {
  print("run: " %+% i)

  # rct_data <- sim_data_v2(300, rct = TRUE)
  # rwd_large <- sim_data_v2(50000, rct = FALSE)
  # rwd_treated <- sample_n(rwd_large[rwd_large$A == 1, ], 600, replace = FALSE)
  # rwd_control <- sample_n(rwd_large[rwd_large$A == 0, ], 600, replace = FALSE)
  # data <- rbind(rct_data, rwd_treated, rwd_control)

  rct_data <- sim_data_v2(300, g_rct, controls_only, rct = TRUE)
  rwd_data <- sim_data_v2(200, g_rct, controls_only, rct = FALSE)
  data <- rbind(rct_data, rwd_data)

  S_node <- 1
  W_node <- 2:4
  A_node <- 5
  Y_node <- 6

  atmle_res <- atmle(data = data,
                     S_node = S_node,
                     W_node = W_node,
                     A_node = A_node,
                     Y_node = Y_node,
                     atmle_pooled = TRUE,
                     controls_only = controls_only,
                     theta_method = "glm",
                     Pi_method = "glm",
                     g_method = "glm",
                     theta_tilde_method = "glm",
                     Q_method = "glm",
                     bias_working_model = "glmnet",
                     pooled_working_model = "glmnet",
                     g_rct = g_rct,
                     family = "gaussian",
                     verbose = FALSE,
                     target_gwt = FALSE)
  atmle_cover[i] <- (atmle_res$lower <= ate) & (atmle_res$upper >= ate)
  atmle_est[i] <- atmle_res$est

  if (atmle_cover[i] == 1) {
    print("covered")
  } else {
    print("not covered")
  }

  escvtmle_res <- ES.cvtmle(txinrwd = !controls_only,
                            data = data,
                            study = "S",
                            covariates = c("W1", "W2", "W3"),
                            treatment_var = "A",
                            treatment = 1,
                            outcome = "Y",
                            pRCT = g_rct,
                            family = "gaussian",
                            Q.SL.library = c("SL.glm"),
                            g.SL.library = c("SL.glm"),
                            Q.discreteSL = TRUE,
                            g.discreteSL = TRUE,
                            V = 5)
  escvtmle_cover[i] <- (escvtmle_res$CI$b2v[1] <= ate) & (escvtmle_res$CI$b2v[2] >= ate)
  escvtmle_est[i] <- escvtmle_res$ATE$b2v
}

print("atmle coverage: " %+% mean(atmle_cover))
print("escvtmle coverage: " %+% mean(escvtmle_cover))
print("atmle mse: " %+% (round((mean(atmle_est)-ate)^2+var(atmle_est),5)))
print("escvtmle mse: " %+% (round((mean(escvtmle_est)-ate)^2+var(escvtmle_est),5)))
