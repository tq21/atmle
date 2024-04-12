library(ggpubr)
library(purrr)
source("utils_plot.R")

# simulation parameters --------------------------------------------------------
B <- 500
n_rct_seq <- c(400, 500, 600, 700, 800)
n_rwd_seq <- n_rct_seq * 5
total_sample_sizes <- n_rct_seq + n_rwd_seq
ate <- 4.2

# HAL a ------------------------------------------------------------------------
atmle_both_res <- readRDS("out/atmle_both_HAL_1_bias_0328.RDS")
escvtmle_res <- readRDS("out/escvtmle_HAL_1_bias_0328.RDS")
tmle_res <- readRDS("out/tmle_HAL_1_bias_0328.RDS")
rct_only_res <- readRDS("out/rct_only_HAL_1_bias_0328.RDS")

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

# average ci length
atmle_ci_length <- get_avg_ci_length(atmle_both_res)
escvtmle_ci_length <- get_avg_ci_length(escvtmle_res)
tmle_ci_length <- get_avg_ci_length(tmle_res)
rct_only_ci_length <- get_avg_ci_length(rct_only_res)

print("A-TMLE vs. ES-CVTMLE = " %+% (round(atmle_ci_length / escvtmle_ci_length, 3)))
print("A-TMLE vs. TMLE = " %+% (round(atmle_ci_length / tmle_ci_length, 3)))
print("A-TMLE vs. RCT ONLY = " %+% (round(atmle_ci_length / rct_only_ci_length, 3)))

# HAL b ------------------------------------------------------------------------
atmle_both_res <- readRDS("out/atmle_both_HAL_2_bias_0328.RDS")
escvtmle_res <- readRDS("out/escvtmle_HAL_2_bias_0328.RDS")
tmle_res <- readRDS("out/tmle_HAL_2_bias_0328.RDS")
rct_only_res <- readRDS("out/rct_only_HAL_2_bias_0328.RDS")

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

# average ci length
atmle_ci_length <- get_avg_ci_length(atmle_both_res)
escvtmle_ci_length <- get_avg_ci_length(escvtmle_res)
tmle_ci_length <- get_avg_ci_length(tmle_res)
rct_only_ci_length <- get_avg_ci_length(rct_only_res)

print("A-TMLE vs. ES-CVTMLE = " %+% (round(atmle_ci_length / escvtmle_ci_length, 3)))
print("A-TMLE vs. TMLE = " %+% (round(atmle_ci_length / tmle_ci_length, 3)))
print("A-TMLE vs. RCT ONLY = " %+% (round(atmle_ci_length / rct_only_ci_length, 3)))

# HAL c -----------------------------------------------------------------------
atmle_both_res <- readRDS("out/atmle_both_HAL_3_bias_0328.RDS")
escvtmle_res <- readRDS("out/escvtmle_HAL_3_bias_0328.RDS")
tmle_res <- readRDS("out/tmle_HAL_3_bias_0328.RDS")
rct_only_res <- readRDS("out/rct_only_HAL_3_bias_0328.RDS")

# MSE
plt_mse_c <- get_mse_plot("",
                          c("A-TMLE", "ES-CVTMLE", "TMLE"),
                          atmle_both_res, escvtmle_res, tmle_res)

# relative MSE
plt_relative_mse_c <- get_relative_mse_plot(
  "",
  c("A-TMLE", "ES-CVTMLE", "TMLE", "RCT ONLY"),
  list(c("A-TMLE", "RCT ONLY"), c("ES-CVTMLE", "RCT ONLY"), c("TMLE", "RCT ONLY")),
  atmle_both_res, escvtmle_res, tmle_res, rct_only_res)

# coverage
plt_coverage_c <- get_cover_plot(
  "",
  c("A-TMLE", "ES-CVTMLE", "TMLE", "RCT ONLY"),
  atmle_both_res, escvtmle_res, tmle_res, rct_only_res)

# proportion of selected
plt_prop_c <- get_plot_prop_selected(escvtmle_res, "")

# average ci length
atmle_ci_length <- get_avg_ci_length(atmle_both_res)
escvtmle_ci_length <- get_avg_ci_length(escvtmle_res)
tmle_ci_length <- get_avg_ci_length(tmle_res)
rct_only_ci_length <- get_avg_ci_length(rct_only_res)

print("A-TMLE vs. ES-CVTMLE = " %+% (round(atmle_ci_length / escvtmle_ci_length, 3)))
print("A-TMLE vs. TMLE = " %+% (round(atmle_ci_length / tmle_ci_length, 3)))
print("A-TMLE vs. RCT ONLY = " %+% (round(atmle_ci_length / rct_only_ci_length, 3)))

plt <- ggarrange(plt_mse_a, plt_relative_mse_a, plt_coverage_a, plt_prop_a,
                 plt_mse_c, plt_relative_mse_c, plt_coverage_c, plt_prop_c,
                 nrow = 2, ncol = 4, common.legend = TRUE)
ggsave(filename = "complex.pdf", plot = plt, device = "pdf",
       path = "figs", width = 16, height = 8, dpi = 300)
