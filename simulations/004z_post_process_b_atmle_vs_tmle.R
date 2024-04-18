library(ggpubr)
library(purrr)
source("utils_plot.R")

# simulation parameters --------------------------------------------------------
B <- 500 # number of runs for each sample size n
n_min <- 500 # smallest sample size
n_max <- 5000 # largest sample size
n_step <- 500 # sample size increment
bA <- 1.5

# atmle for psi_tilde
atmle_res <- readRDS("out/atmle_psi_tilde_positivity_922.RDS")
bias_atmle <- unlist(map(atmle_res$all_psi_est, function(.x) mean(.x-bA)))
se_atmle <- unlist(map(atmle_res$all_psi_est, function(.x) sd(.x)))
var_atmle <- unlist(map(atmle_res$all_psi_est, function(.x) var(.x)))
mse_atmle <- bias_atmle^2 + var_atmle
cover_atmle <- unlist(map(atmle_res$all_psi_coverage, function(.x) mean(.x)))

# tmle for psi_tilde
tmle_res <- readRDS("out/tmle_psi_tilde_positivity_922.RDS")
bias_tmle <- unlist(map(tmle_res$all_psi_est, function(.x) mean(.x-bA)))
se_tmle <- unlist(map(tmle_res$all_psi_est, function(.x) sd(.x)))
var_tmle <- unlist(map(tmle_res$all_psi_est, function(.x) var(.x)))
mse_tmle <- bias_tmle^2 + var_tmle
cover_tmle <- unlist(map(tmle_res$all_psi_coverage, function(.x) mean(.x)))

dt <- data.frame(n = seq(n_min, n_max, n_step),
                 estimator = rep(c("A-TMLE", "TMLE"), each = length(bias_atmle)),
                 bias = c(bias_atmle, bias_tmle),
                 se = c(se_atmle, se_tmle),
                 mse = c(mse_atmle, mse_tmle),
                 coverage = c(cover_atmle, cover_tmle))
dt <- dt[order(dt$n), ]

p_bias <- ggplot(dt, aes(x = n, y = abs(bias), color = estimator)) +
  geom_point() +
  geom_line() +
  labs(title = "",
       x = "n",
       y = "bias") +
  theme_minimal() +
  theme(text = element_text(size = 16))

p_se <- ggplot(dt, aes(x = n, y = se, color = estimator)) +
  geom_point() +
  geom_line() +
  labs(title = "",
       x = "n",
       y = "standard error") +
  theme_minimal() +
  theme(text = element_text(size = 16))

p_mse <- ggplot(dt, aes(x = n, y = mse, color = estimator)) +
  geom_point() +
  geom_line() +
  labs(title = "",
       x = "n",
       y = "mse") +
  theme_minimal() +
  theme(text = element_text(size = 16))

p_cover <- ggplot(dt, aes(x = n, y = coverage, color = estimator)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0.95, color = "red", linetype = "dashed", linewidth = 1) +
  scale_y_continuous(breaks = seq(0.8, 1.0, 0.05), limits = c(0.8, 1.0)) +
  labs(title = "",
       x = "n",
       y = "coverage") +
  theme_minimal() +
  theme(text = element_text(size = 16))

# relative mse plot
dt_relative <- data.frame(n = seq(n_min, n_max, n_step),
                          names = rep(c("mse(TMLE)/mse(A-TMLE)"),
                                      each = length(seq(n_min, n_max, n_step))),
                          ratio = c(mse_tmle / mse_atmle))
relative_mse_plot <- ggplot(dt_relative, aes(x = n, y = ratio, color = names)) +
  geom_point(size = 1.5) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 1, color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "",
       x = "n",
       y = "relative mse") +
  theme_minimal() +
  theme(text = element_text(size = 16),
        plot.title = element_text(hjust = 0.5))

plt <- ggarrange(p_bias, p_se, p_mse, p_cover, relative_mse_plot,
                 nrow = 2, ncol = 3, common.legend = TRUE)
plt <- annotate_figure(plt, top = text_grob("pooled-ATE, positivity violation",
                                            face = "bold", size = 16))

ggsave(filename = "pooled_ate_atmle_vs_tmle_positivity.pdf", plot = plt, device = "pdf",
       path = "plot", width = 12, height = 8, dpi = 300)
