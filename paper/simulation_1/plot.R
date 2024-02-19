library(ggpubr)
library(purrr)
source("utils_plot.R")

# simulation parameters --------------------------------------------------------
B <- 500
n_rct_seq <- c(400, 500, 600, 700, 800)
n_rwd_seq <- n_rct_seq * 3
total_sample_sizes <- n_rct_seq + n_rwd_seq
ate <- 1.5

# bias a -----------------------------------------------------------------------
atmle_both_res <- readRDS("out/atmle_both_a_bias_0218.RDS")
escvtmle_res <- readRDS("out/escvtmle_a_bias_0218.RDS")
tmle_res <- readRDS("out/tmle_a_bias_0218.RDS")
rct_only_res <- readRDS("out/rct_only_a_bias_0218.RDS")
plt_bias_a <- get_mse_plot("No bias",
                           c("A-TMLE", "ES-CVTMLE", "TMLE"),
                           atmle_both_res, escvtmle_res, tmle_res)

# relative MSE
plt_relative_mse_a <- get_relative_mse_plot(
  "A bias",
  c("A-TMLE", "ES-CVTMLE", "TMLE", "RCT ONLY"),
  list(c("A-TMLE", "RCT ONLY"), c("ES-CVTMLE", "RCT ONLY"), c("TMLE", "RCT ONLY")),
  atmle_both_res, escvtmle_res, tmle_res, rct_only_res)

# bias b -----------------------------------------------------------------------
atmle_both_res <- readRDS("out/atmle_both_b_bias_0218.RDS")
escvtmle_res <- readRDS("out/escvtmle_b_bias_0218.RDS")
tmle_res <- readRDS("out/tmle_b_bias_0218.RDS")
rct_only_res <- readRDS("out/rct_only_b_bias_0218.RDS")
plt_bias_b <- get_mse_plot("No bias",
                           c("A-TMLE", "ES-CVTMLE", "TMLE"),
                           atmle_both_res, escvtmle_res, tmle_res)

# relative MSE
plt_relative_mse_b <- get_relative_mse_plot(
  "B bias",
  c("A-TMLE", "ES-CVTMLE", "TMLE", "RCT ONLY"),
  list(c("A-TMLE", "RCT ONLY"), c("ES-CVTMLE", "RCT ONLY"), c("TMLE", "RCT ONLY")),
  atmle_both_res, escvtmle_res, tmle_res, rct_only_res)



# ES-CVTMLE proportion of folds selected, TODO
plt_no_bias_prop <- get_plot_selected(escvtmle_res, "No bias")
ggsave(filename = "no_bias_prop.pdf", plot = plt_no_bias_prop, device = "pdf",
       path = "plot", width = 10, height = 8, dpi = 300)

