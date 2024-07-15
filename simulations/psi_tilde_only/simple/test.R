library(devtools)
library(purrr)
source("sim_data.R")
load_all()

t0 <- 3
data_A1 <- sim_data(100000, A_counter = 1)
data_A0 <- sim_data(100000, A_counter = 0)
truth <- mean(data_A1$T_tilde > t0)-mean(data_A0$T_tilde > t0)

n <- 2000
B <- 500

all_psi_tilde_r_learner <- numeric(B)
all_psi_tilde_no_tmle_lambda <- numeric(B)
all_psi_tilde_tmle_lambda <- numeric(B)
all_psi_tilde_r_learner_tmle_beta <- numeric(B)
all_psi_tilde_tmle_lambda_cover <- numeric(B)
all_psi_tilde_r_learner_tmle_beta_cover <- numeric(B)

for (b in 1:B) {
  print(b)
  data <- sim_data(n)
  res <- atmle_surv(data = data,
                    W = c("W1", "W2"),
                    A = "A",
                    T_tilde = "T_tilde",
                    tau = 5,
                    Delta = "Delta",
                    t0 = t0,
                    g_method = "glm",
                    lambda_method = "glm")
  all_psi_tilde_r_learner[b] <- res$psi_tilde_r_learner
  all_psi_tilde_no_tmle_lambda[b] <- res$psi_tilde_no_tmle_lambda
  all_psi_tilde_tmle_lambda[b] <- res$psi_tilde_tmle_lambda
  all_psi_tilde_r_learner_tmle_beta[b] <- res$psi_tilde_r_learner_tmle_beta

  if (res$psi_tilde_tmle_lambda_lower <= truth & res$psi_tilde_tmle_lambda_upper >= truth) {
    all_psi_tilde_tmle_lambda_cover[b] <- 1
  } else {
    all_psi_tilde_tmle_lambda_cover[b] <- 0
  }

  if (res$psi_tilde_r_learner_tmle_beta_lower <= truth & res$psi_tilde_r_learner_tmle_beta_upper >= truth) {
    all_psi_tilde_r_learner_tmle_beta_cover[b] <- 1
  } else {
    all_psi_tilde_r_learner_tmle_beta_cover[b] <- 0
  }

  print("TMLE: " %+% round(mean(all_psi_tilde_tmle_lambda_cover[1:b]), 2) %+%
          ", TMLE beta: " %+% round(mean(all_psi_tilde_r_learner_tmle_beta_cover[1:b]), 2))
}

par(mfrow = c(2, 2))
hist(all_psi_tilde_r_learner, main = "IPCW R-learner")
abline(v = truth, col = "red")
hist(all_psi_tilde_no_tmle_lambda, main = "No TMLE lambda")
abline(v = truth, col = "red")
hist(all_psi_tilde_tmle_lambda, main = "TMLE lambda")
abline(v = truth, col = "red")
hist(all_psi_tilde_r_learner_tmle_beta, main = "IPCW R-learner TMLE beta")
abline(v = truth, col = "red")

tmp <- get_true_cate(data, t0)
mean((res$ipcw_r_loss_cate-tmp)^2)
mean((res$tmle_cate-tmp)^2)
