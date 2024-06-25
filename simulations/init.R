library(devtools)
library(data.table)
library(hal9001)
library(glmnet)
library(purrr)
load_all()
set.seed(213)
source("sim_data.R")

t0 <- 3
data_A1 <- sim_data(100000, A_counter = 1)
data_A0 <- sim_data(100000, A_counter = 0)
truth <- mean(data_A1$T_tilde > t0)-mean(data_A0$T_tilde > t0)

n_seq <- c(1000, 2000, 3000, 4000)
B <- 200

all_psi_tilde_r_learner <- vector(mode = "list", length = length(n_seq))
all_psi_tilde_no_tmle_lambda <- vector(mode = "list", length = length(n_seq))
all_psi_tilde_tmle_lambda <- vector(mode = "list", length = length(n_seq))
all_psi_tilde_r_learner_tmle_beta <- vector(mode = "list", length = length(n_seq))

for (n_idx in 1:length(n_seq)) {
  n <- n_seq[n_idx]
  print(n)
  for (b in 1:B) {
    print(b)
    data <- sim_data(n)
    res <- atmle_surv(data = data,
                      S = "T_tilde",
                      W = c("W1", "W2"),
                      A = "A",
                      T_tilde = "T_tilde",
                      Delta = "Delta",
                      t0 = t0,
                      g_rct = 0.5,
                      controls_only = FALSE,
                      g_method = "glm",
                      lambda_method = "glm")
    all_psi_tilde_r_learner[[n_idx]][b] <- res$psi_tilde_r_learner
    all_psi_tilde_no_tmle_lambda[[n_idx]][b] <- res$psi_tilde_no_tmle_lambda
    all_psi_tilde_tmle_lambda[[n_idx]][b] <- res$psi_tilde_tmle_lambda
    all_psi_tilde_r_learner_tmle_beta[[n_idx]][b] <- res$psi_tilde_r_learner_tmle_beta
  }
}

save(list = c("all_psi_tilde_r_learner",
              "all_psi_tilde_no_tmle_lambda",
              "all_psi_tilde_tmle_lambda",
              "all_psi_tilde_r_learner_tmle_beta"), file = "out/results_0625.RData")
