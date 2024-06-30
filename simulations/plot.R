library(ggplot2)
library(data.table)
library(purrr)
load("out/results_0630_fixed.RData")
source("sim_data.R")

set.seed(238974)
t0 <- 4
data_A1 <- sim_data(100000, A_counter = 1)
data_A0 <- sim_data(100000, A_counter = 0)
truth <- mean(data_A1$T_tilde > t0)-mean(data_A0$T_tilde > t0)

# diagnostic plots
par(mfrow = c(1,3))
hist(all_psi_tilde_r_learner[[1]], main = "IPCW R-learner", breaks = 20, xlab = "")
abline(v = truth, col = "red")
hist(all_psi_tilde_no_tmle_lambda[[1]], main = "Initial lambda", breaks = 20, xlab = "")
abline(v = truth, col = "red")
hist(all_psi_tilde_r_learner_tmle_beta[[1]], main = "IPCW R-learner TMLE beta", breaks = 20, xlab = "")
abline(v = truth, col = "red")

df_plot <- data.frame(n = seq(500, 2000, 500),
                      r_learner_mse = map_dbl(all_psi_tilde_r_learner, ~mean((.x-truth)^2)),
                      r_learner_tmle_mse = map_dbl(all_psi_tilde_tmle_lambda, ~mean((.x-truth)^2)))

ggplot(data = df_plot, aes(x = n)) +
  geom_line(aes(y = r_learner_mse, color = "IPCW R-learner")) +
  geom_line(aes(y = r_learner_tmle_mse, color = "IPCW R-learner + TMLE")) +
  labs(title = "MSE of R-learner and R-learner + TMLE",
       y = "MSE",
       x = "Sample size") +
  theme_minimal() +
  scale_color_manual(values = c("IPCW R-learner" = "blue", "IPCW R-learner + TMLE" = "red"))
