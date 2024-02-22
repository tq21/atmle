library(ggpubr)
library(purrr)
source("utils_plot.R")

# simulation parameters --------------------------------------------------------
B <- 500
n_rct_seq <- c(400, 500, 600, 700, 800)
n_rwd_seq <- n_rct_seq * 3
total_sample_sizes <- n_rct_seq + n_rwd_seq
ate <- 1.5

# HAL 1 ------------------------------------------------------------------------
atmle_both_res <- readRDS("out/atmle_both_HAL_1_bias_0218.RDS")
escvtmle_res <- readRDS("out/escvtmle_HAL_1_bias_0218.RDS")
tmle_res <- readRDS("out/tmle_HAL_1_bias_0218.RDS")
rct_only_res <- readRDS("out/rct_only_HAL_1_bias_0218.RDS")

# MSE
plt_mse_a <- get_mse_plot("No bias",
                          c("A-TMLE", "ES-CVTMLE", "TMLE"),
                          atmle_both_res, escvtmle_res, tmle_res)

# relative MSE
plt_relative_mse_a <- get_relative_mse_plot(
  "",
  c("A-TMLE", "ES-CVTMLE", "TMLE", "RCT ONLY"),
  list(c("A-TMLE", "RCT ONLY"), c("ES-CVTMLE", "RCT ONLY"), c("TMLE", "RCT ONLY")),
  atmle_both_res, escvtmle_res, tmle_res, rct_only_res)

# coverage
plt_coverage_a <- get_cover_plot(
  "",
  c("A-TMLE", "ES-CVTMLE", "TMLE", "RCT ONLY"),
  atmle_both_res, escvtmle_res, tmle_res, rct_only_res)

# proportion of selected
plt_prop_a <- get_plot_prop_selected(escvtmle_res, "")

# bias b -----------------------------------------------------------------------
atmle_both_res <- readRDS("out/atmle_both_HAL_2_bias_0218.RDS")
escvtmle_res <- readRDS("out/escvtmle_HAL_2_bias_0218.RDS")
tmle_res <- readRDS("out/tmle_HAL_2_bias_0218.RDS")
rct_only_res <- readRDS("out/rct_only_HAL_2_bias_0218.RDS")

# MSE
plt_mse_b <- get_mse_plot("No bias",
                          c("A-TMLE", "ES-CVTMLE", "TMLE"),
                          atmle_both_res, escvtmle_res, tmle_res)

# relative MSE
plt_relative_mse_b <- get_relative_mse_plot(
  "",
  c("A-TMLE", "ES-CVTMLE", "TMLE", "RCT ONLY"),
  list(c("A-TMLE", "RCT ONLY"), c("ES-CVTMLE", "RCT ONLY"), c("TMLE", "RCT ONLY")),
  atmle_both_res, escvtmle_res, tmle_res, rct_only_res)

# coverage
plt_coverage_b <- get_cover_plot(
  "",
  c("A-TMLE", "ES-CVTMLE", "TMLE", "RCT ONLY"),
  atmle_both_res, escvtmle_res, tmle_res, rct_only_res)

# proportion of selected
plt_prop_b <- get_plot_prop_selected(escvtmle_res, "")

plt <- ggarrange(plt_mse_a, plt_relative_mse_a, plt_coverage_a, plt_prop_a,
                 plt_mse_b, plt_relative_mse_b, plt_coverage_b, plt_prop_b,
                 nrow = 2, ncol = 4, common.legend = TRUE)
ggsave(filename = "HAL.pdf", plot = plt, device = "pdf",
       path = "figs", width = 16, height = 8, dpi = 300)
