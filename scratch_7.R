#library(atmle)
library(EScvtmle)
library(dplyr)
library(hal9001)
devtools::load_all()

`%+%` <- function(x, y) paste0(x, y)

sim_data <- function(n, gamma, rct) {
  W1 <- round(rnorm(n, 0, 1), 2)
  W2 <- rbinom(n, 1, 0.5)
  W3 <- round(rnorm(n, 0.1*W2, 1), 2)
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

  Y <- rnorm(n, (0.2+0.1*W3)*A+U1+U2, 1)

  data <- data.frame(S = S,
                     W1 = W1,
                     W2 = W2,
                     W3 = W3,
                     A = A,
                     Y = Y)

  return(data)
}

# SIMULATION PARAMETERS
controls_only <- FALSE
B <- 100
ate <- 0.205
g_rct <- 0.5
atmle_cover <- numeric(B)
atmle_est <- numeric(B)
atmle_pooled_est <- numeric(B)
atmle_bias_est <- numeric(B)
escvtmle_cover <- numeric(B)
escvtmle_est <- numeric(B)
gamma <- 0.5

for (i in 1:B) {
  print("run: " %+% i)

  rct_data <- sim_data(300, gamma, TRUE)
  rwd_data <- sim_data(1200, gamma, FALSE)
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
                     bias_working_model = "HAL",
                     pooled_working_model = "glmnet",
                     g_rct = g_rct,
                     family = "gaussian",
                     verbose = FALSE,
                     target_gwt = TRUE,
                     min_working_model = FALSE,
                     undersmooth = FALSE)

  atmle_cover[i] <- (atmle_res$lower <= ate) & (atmle_res$upper >= ate)
  atmle_est[i] <- atmle_res$est
  atmle_pooled_est[i] <- atmle_res$psi_tilde_est
  atmle_bias_est[i] <- atmle_res$psi_pound_est

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

library(purrr)
oracle_lower <- atmle_est - 1.96*sd(atmle_est)
oracle_upper <- atmle_est + 1.96*sd(atmle_est)
oracle_coverage <- mean(map2_vec(oracle_lower, oracle_upper, function(.x, .y) {
  return(ate >= .x & ate <= .y)
})) # 0.925

# save results
# save(atmle_cover,
#      atmle_est,
#      atmle_pooled_est,
#      atmle_bias_est,
#      escvtmle_cover,
#      escvtmle_est,
#      file = "out/HAL_slow.RData")
