library(atmle)
library(EScvtmle)
library(devtools)
load_all()

`%+%` <- function(a, b) paste0(a, b)

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

  # Y <- (0.2+0.1*W3)*A+(1-S)*(10+rnorm(n, 0, 0.1))+U1+U2+rnorm(n, 0, 1)
  Y <- (0.2+0.1*W3)*A+(1-S)*(-0.5*A+2.1*W1+rnorm(n, 0, 0.1))+rnorm(n, 0, 1)

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
ate <- 0.205
g_rct <- 0.5
n_rct <- 1000
n_rwd <- 5000
controls_only <- FALSE
S_node <- 1
W_node <- 2:4
A_node <- 5
Y_node <- 6

est_atmle <- numeric(length = B)
est_escvtmle <- numeric(length = B)
est_rct <- numeric(length = B)
est_pooled <- numeric(length = B)
psi_pound_atmle <- numeric(length = B)
psi_tilde_atmle <- numeric(length = B)

# truths
ate <- 0.205
n_large <- 10^6
rct_large_A1 <- sim_data(n_large, 0, rct = TRUE, A_counter = 1)
rct_large_A0 <- sim_data(n_large, 0, rct = TRUE, A_counter = 0)
rwd_large_A1 <- sim_data(n_large, 0, rct = FALSE, A_counter = 1)
rwd_large_A0 <- sim_data(n_large, 0, rct = FALSE, A_counter = 0)
mean(rct_large_A1$Y)-mean(rwd_large_A1$Y)
mean(rct_large_A0$Y)-mean(rwd_large_A0$Y)

psi_pound <-

for (b in 1:B) {
  print("run: " %+% b)

  # simulate data
  data <- rbind(sim_data(n_rct, g_rct, TRUE),
                sim_data(n_rwd, g_rct, FALSE))

  # tmle on rct
  tmle_rct_res <- tmle(Y = data[data$S == 1, Y_node],
                       A = data[data$S == 1, A_node],
                       W = data[data$S == 1, W_node],
                       Qform = "Y ~ A + A*W3",
                       g.SL.library = c("SL.glm"),
                       Q.SL.library = c("SL.glm"),
                       family = "gaussian")

  # tmle on pooled
  tmle_pooled_res <- tmle(Y = data[data$S == 0, Y_node],
                          A = data[data$S == 0, A_node],
                          W = data[data$S == 0, W_node],
                          Qform = "Y ~ A + A*W3",
                          g.SL.library = c("SL.glm"),
                          Q.SL.library = c("SL.glm"),
                          family = "gaussian")

  # atmle with glmnet
  atmle_res <- atmle(data = data,
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

  # ES-CVTMLE
  escvtmle_res <- ES.cvtmle(txinrwd = TRUE,
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

  # collect results
  est_rct[i] <- tmle_rct_res$estimates$ATE$psi
  est_pooled[i] <- tmle_pooled_res$estimates$ATE$psi
  est_atmle[i] <- atmle_res$est
  psi_pound_atmle[i] <- atmle_res$psi_pound
  psi_tilde_atmle[i] <- atmle_res$psi_tilde
  est_escvtmle[i] <- escvtmle_res$ATE$b2v
}

save(list = c("est_atmle_glmnet",
              "est_escvtmle",
              "est_rct_only",
              "est_rwd"), file = "scm_replicate_larger_n_rct.RData")
