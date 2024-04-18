library(ggpubr)

data <- data.frame(n = n_seq,
                   #no_bias_bias = unlist(map(no_bias$all_psi_est, ~ (mean(.x) - truth)^2)),
                   #no_bias_se = unlist(map(no_bias$all_psi_est, ~ sd(.x))),
                   #no_bias_coverage = unlist(map(no_bias$all_psi_coverage, ~ mean(.x))),
                   #small_bias_bias = unlist(map(small_bias$all_psi_est, ~ (mean(.x) - truth)^2)),
                   #small_bias_se = unlist(map(small_bias$all_psi_est, ~ sd(.x))),
                   #small_bias_coverage = unlist(map(small_bias$all_psi_coverage, ~ mean(.x))),
                   #large_bias_bias = unlist(map(large_bias$all_psi_est, ~ (mean(.x) - truth)^2)),
                   #large_bias_se = unlist(map(large_bias$all_psi_est, ~ sd(.x))),
                   #large_bias_coverage = unlist(map(large_bias$all_psi_coverage, ~ mean(.x))),
                   param_bias_bias = unlist(map(param_bias$all_psi_est, ~ (mean(.x) - truth)^2)),
                   param_bias_se = unlist(map(param_bias$all_psi_est, ~ sd(.x))),
                   param_bias_coverage = unlist(map(param_bias$all_psi_coverage, ~ mean(.x))))

# no bias ----------------------------------------------------------------------
p_no_bias_bias <- ggplot(data, aes(x = n, y = no_bias_bias)) +
  geom_point(color = "blue") +
  geom_line(color = "blue") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "",
       x = "n",
       y = "bias") +
  theme_minimal() +
  theme(text = element_text(size = 16))

p_no_bias_se <- ggplot(data, aes(x = n, y = no_bias_se)) +
  geom_point(color = "blue") +
  geom_line(color = "blue") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
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

# small bias -------------------------------------------------------------------
p_small_bias_bias <- ggplot(data, aes(x = n, y = small_bias_bias)) +
  geom_point(color = "blue") +
  geom_line(color = "blue") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "",
       x = "n",
       y = "bias") +
  theme_minimal() +
  theme(text = element_text(size = 16))

p_small_bias_se <- ggplot(data, aes(x = n, y = small_bias_se)) +
  geom_point(color = "blue") +
  geom_line(color = "blue") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "",
       x = "n",
       y = "standard error") +
  theme_minimal() +
  theme(text = element_text(size = 16))

p_small_bias_coverage <- ggplot(data, aes(x = n, y = small_bias_coverage)) +
  geom_point(color = "blue") +
  geom_line(color = "blue") +
  geom_hline(yintercept = 0.95, color = "red", linetype = "dashed", linewidth = 1) +
  scale_y_continuous(breaks = seq(0.8, 1.0, 0.05), limits = c(0.8, 1.0)) +
  labs(title = "",
       x = "n",
       y = "coverage") +
  theme_minimal() +
  theme(text = element_text(size = 16))

# large bias -------------------------------------------------------------------
p_large_bias_bias <- ggplot(data, aes(x = n, y = large_bias_bias)) +
  geom_point(color = "blue") +
  geom_line(color = "blue") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "",
       x = "n",
       y = "bias") +
  theme_minimal() +
  theme(text = element_text(size = 16))

p_large_bias_se <- ggplot(data, aes(x = n, y = large_bias_se)) +
  geom_point(color = "blue") +
  geom_line(color = "blue") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "",
       x = "n",
       y = "standard error") +
  theme_minimal() +
  theme(text = element_text(size = 16))

p_large_bias_coverage <- ggplot(data, aes(x = n, y = large_bias_coverage)) +
  geom_point(color = "blue") +
  geom_line(color = "blue") +
  geom_hline(yintercept = 0.95, color = "red", linetype = "dashed", linewidth = 1) +
  scale_y_continuous(breaks = seq(0.8, 1.0, 0.05), limits = c(0.8, 1.0)) +
  labs(title = "",
       x = "n",
       y = "coverage") +
  theme_minimal() +
  theme(text = element_text(size = 16))

# param bias -------------------------------------------------------------------
p_param_bias_bias <- ggplot(data, aes(x = n, y = param_bias_bias)) +
  geom_point(color = "blue") +
  geom_line(color = "blue") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "",
       x = "n",
       y = "bias") +
  theme_minimal() +
  theme(text = element_text(size = 16))

p_param_bias_se <- ggplot(data, aes(x = n, y = param_bias_se)) +
  geom_point(color = "blue") +
  geom_line(color = "blue") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "",
       x = "n",
       y = "standard error") +
  theme_minimal() +
  theme(text = element_text(size = 16))

p_param_bias_coverage <- ggplot(data, aes(x = n, y = param_bias_coverage)) +
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
                 p_small_bias_bias, p_small_bias_se, p_small_bias_coverage,
                 p_large_bias_bias, p_large_bias_se, p_large_bias_coverage,
                 p_param_bias_bias, p_param_bias_se, p_param_bias_coverage,
                 nrow = 4, ncol = 3)
