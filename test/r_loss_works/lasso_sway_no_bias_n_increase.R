library(ggpubr)
library(purrr)
source("utils.R")

# n increase -------------------------------------------------------------------
B <- 500
n_min <- 200
n_max <- 3000
n_step <- 200

no_bias <- run_sim_n_increase(B = B,
                              n_min = n_min,
                              n_max = n_max,
                              n_step = n_step,
                              bA = 1.5,
                              bias = 0,
                              nuisance_method = "lasso",
                              working_model = "lasso",
                              verbose = TRUE)

constant_bias <- run_sim_n_increase(B = B,
                                    n_min = n_min,
                                    n_max = n_max,
                                    n_step = n_step,
                                    bA = 1.5,
                                    bias = 1.8,
                                    nuisance_method = "lasso",
                                    working_model = "lasso",
                                    verbose = TRUE)

n_seq <- seq(n_min, n_max, n_step)
truth <- 1.5

data <- data.frame(n = n_seq,
                   no_bias_bias = unlist(map(no_bias$all_psi_est, ~ (mean(.x) - truth)^2)),
                   no_bias_se = unlist(map(no_bias$all_psi_est, ~ sd(.x))),
                   no_bias_coverage = unlist(map(no_bias$all_psi_coverage, ~ mean(.x))),
                   constant_bias_bias = unlist(map(constant_bias$all_psi_est, ~ (mean(.x) - truth)^2)),
                   constant_bias_se = unlist(map(constant_bias$all_psi_est, ~ sd(.x))),
                   constant_bias_coverage = unlist(map(constant_bias$all_psi_coverage, ~ mean(.x))))

# no bias ----------------------------------------------------------------------
p_no_bias_bias <- ggplot(data, aes(x = n, y = no_bias_bias)) +
  geom_point(color = "blue") +
  geom_line(color = "blue") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  #scale_y_continuous(breaks = seq(1.2, 1.8, 0.1), limits = c(1.2, 1.8)) +
  labs(title = "",
       x = "n",
       y = "bias") +
  theme_minimal() +
  theme(text = element_text(size = 16))

p_no_bias_se <- ggplot(data, aes(x = n, y = no_bias_se)) +
  geom_point(color = "blue") +
  geom_line(color = "blue") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  #scale_y_continuous(breaks = seq(1.2, 1.8, 0.1), limits = c(1.2, 1.8)) +
  labs(title = "",
       x = "n",
       y = "standard error") +
  theme_minimal() +
  theme(text = element_text(size = 16))

p_no_bias_coverage <- ggplot(data, aes(x = n, y = no_bias_coverage)) +
  geom_point(color = "blue") +
  geom_line(color = "blue") +
  geom_hline(yintercept = 0.95, color = "red", linetype = "dashed", linewidth = 1) +
  scale_y_continuous(breaks = seq(0.8, 1.0, 0.05), limits = c(0.8, 1.0)) +
  labs(title = "",
       x = "n",
       y = "coverage") +
  theme_minimal() +
  theme(text = element_text(size = 16))

# constant bias ----------------------------------------------------------------
p_constant_bias_bias <- ggplot(data, aes(x = n, y = constant_bias_bias)) +
  geom_point(color = "blue") +
  geom_line(color = "blue") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "",
       x = "n",
       y = "bias") +
  theme_minimal() +
  theme(text = element_text(size = 16))

p_constant_bias_se <- ggplot(data, aes(x = n, y = constant_bias_se)) +
  geom_point(color = "blue") +
  geom_line(color = "blue") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "",
       x = "n",
       y = "standard error") +
  theme_minimal() +
  theme(text = element_text(size = 16))

p_constant_bias_coverage <- ggplot(data, aes(x = n, y = constant_bias_coverage)) +
  geom_point(color = "blue") +
  geom_line(color = "blue") +
  geom_hline(yintercept = 0.95, color = "red", linetype = "dashed", linewidth = 1) +
  scale_y_continuous(breaks = seq(0.8, 1.0, 0.05), limits = c(0.8, 1.0)) +
  labs(title = "",
       x = "n",
       y = "coverage") +
  theme_minimal() +
  theme(text = element_text(size = 16))

plt <- ggarrange(p_no_bias_bias, p_no_bias_se, p_no_bias_coverage,
                 p_constant_bias_bias, p_constant_bias_se, p_constant_bias_coverage,
                 nrow = 2, ncol = 3)

save(list = c("no_bias", "constant_bias", "n_seq"),
     file = "out/no_bias_and_constant_bias_8_29.RData")

no_bias_n_increase <- run_sim_n_increase(B = 1,
                                         n_min = 200,
                                         n_max = 5000,
                                         n_step = 50,
                                         bA = 1.5,
                                         bias = 0,
                                         nuisance_method = "lasso",
                                         working_model = "lasso",
                                         verbose=TRUE)

constant_bias_n_increase <- run_sim_n_increase(B = 1,
                                               n_min = 200,
                                               n_max = 5000,
                                               n_step = 50,
                                               bA = 1.5,
                                               bias = 1.8,
                                               nuisance_method = "lasso",
                                               working_model = "lasso",
                                               verbose=TRUE)

data_n <- data.frame(n = seq(200, 5000, 50),
                     no_bias = unlist(no_bias_n_increase$all_psi_est),
                     constant_bias = unlist(constant_bias_n_increase$all_psi_est))

p_no_bias_bias_converge <- ggplot(data_n, aes(x = n, y = no_bias)) +
  geom_point(color = "blue") +
  geom_line(color = "blue") +
  geom_hline(yintercept = truth, color = "red", linetype = "dashed", linewidth = 1) +
  #scale_y_continuous(breaks = seq(1.2, 1.8, 0.1), limits = c(1.2, 1.8)) +
  labs(title = "",
       x = "n",
       y = "bias") +
  theme_minimal() +
  theme(text = element_text(size = 16))

p_constant_bias_bias_converge <- ggplot(data_n, aes(x = n, y = constant_bias)) +
  geom_point(color = "blue") +
  geom_line(color = "blue") +
  geom_hline(yintercept = truth, color = "red", linetype = "dashed", linewidth = 1) +
  #scale_y_continuous(breaks = seq(1.2, 1.8, 0.1), limits = c(1.2, 1.8)) +
  labs(title = "",
       x = "n",
       y = "bias") +
  theme_minimal() +
  theme(text = element_text(size = 16))


