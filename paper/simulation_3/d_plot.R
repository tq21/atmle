library(ggpubr)
library(purrr)
source("utils_plot.R")

# simulation parameters --------------------------------------------------------
B <- 1000 # number of runs for each sample size n
n_min <- 0.05 # smallest sample size
n_max <- 0.3 # largest sample size
n_step <- 0.05 # sample size increment
ate <- 1.5

# a. no bias -------------------------------------------------------------------
atmle_both_res <- readRDS("out/atmle_both_no_bias_1217.RDS")
escvtmle_res <- readRDS("out/escvtmle_no_bias_1217.RDS")
tmle_res <- readRDS("out/tmle_no_bias_1217.RDS")
plt_no_bias <- get_plot(
  "No bias",
  c("A-TMLE", "A-TMLE*", "ES-CVTMLE", "TMLE", "RCT ONLY"),
  atmle_both_res, atmle_tmle_res, escvtmle_res, tmle_res, rct_only_res)
ggsave(filename = "no_bias.pdf", plot = plt_no_bias, device = "pdf",
       path = "plot", width = 10, height = 8, dpi = 300)

# all comparisons relative MSE
plt_relative_mse_no_bias <- get_relative_mse_plot(
  "No bias",
  c("A-TMLE", "A-TMLE*", "ES-CVTMLE", "TMLE", "RCT ONLY"),
  list(c("A-TMLE", "RCT ONLY"), c("A-TMLE*", "RCT ONLY"), c("ES-CVTMLE", "RCT ONLY"), c("TMLE", "RCT ONLY")),
  atmle_both_res, atmle_tmle_res, escvtmle_res, tmle_res, rct_only_res)
ggsave(filename = "no_bias_relative.pdf", plot = plt_relative_mse_no_bias, device = "pdf",
       path = "plot", width = 10, height = 8, dpi = 300)

plt_no_bias_prop <- get_plot_selected(escvtmle_res, "No bias")
ggsave(filename = "no_bias_prop.pdf", plot = plt_no_bias_prop, device = "pdf",
       path = "plot", width = 10, height = 8, dpi = 300)

# b. small bias -------------------------------------------------------------
atmle_both_res <- readRDS("out/atmle_both_small_bias_1217.RDS")
atmle_tmle_res <- readRDS("out/atmle_tmle_small_bias_1217.RDS")
escvtmle_res <- readRDS("out/escvtmle_small_bias_1217.RDS")
tmle_res <- readRDS("out/tmle_small_bias_1217.RDS")
rct_only_res <- readRDS("out/rct_only_small_bias_1217.RDS")
plt_param_simple_bias <- get_plot(
  "Small bias",
  c("A-TMLE", "A-TMLE*", "ES-CVTMLE", "TMLE", "RCT ONLY"),
  atmle_both_res, atmle_tmle_res, escvtmle_res, tmle_res, rct_only_res)
ggsave(filename = "param_simple_bias.pdf", plot = plt_param_simple_bias, device = "pdf",
       path = "plot", width = 10, height = 8, dpi = 300)

# all comparisons relative MSE
plt_relative_mse_param_simple_bias <- get_relative_mse_plot(
  "Small bias",
  c("A-TMLE", "A-TMLE*", "ES-CVTMLE", "TMLE", "RCT ONLY"),
  list(c("A-TMLE", "RCT ONLY"), c("A-TMLE*", "RCT ONLY"), c("ES-CVTMLE", "RCT ONLY"), c("TMLE", "RCT ONLY")),
  atmle_both_res, atmle_tmle_res, escvtmle_res, tmle_res, rct_only_res)
ggsave(filename = "param_simple_bias_relative.pdf", plot = plt_relative_mse_param_simple_bias, device = "pdf",
       path = "plot", width = 10, height = 8, dpi = 300)

plt_param_simple_bias_prop <- get_plot_selected(escvtmle_res, "Simple parametric bias")
ggsave(filename = "param_simple_bias_prop.pdf", plot = plt_param_simple_bias_prop, device = "pdf",
       path = "plot", width = 10, height = 8, dpi = 300)

# c. large bias -------------------------------------------------------------
atmle_both_res <- readRDS("out/atmle_both_large_bias_1217.RDS")
atmle_tmle_res <- readRDS("out/atmle_tmle_large_bias_1217.RDS")
escvtmle_res <- readRDS("out/escvtmle_large_bias_1217.RDS")
tmle_res <- readRDS("out/tmle_large_bias_1217.RDS")
rct_only_res <- readRDS("out/rct_only_large_bias_1217.RDS")
plt_param_complex_bias <- get_plot(
  "Large bias",
  c("A-TMLE", "A-TMLE*", "ES-CVTMLE", "TMLE", "RCT ONLY"),
  atmle_both_res, atmle_tmle_res, escvtmle_res, tmle_res, rct_only_res)
ggsave(filename = "param_complex_bias.pdf", plot = plt_param_complex_bias, device = "pdf",
       path = "plot", width = 10, height = 8, dpi = 300)

# all comparisons relative MSE
plt_relative_mse_param_complex_bias <- get_relative_mse_plot(
  "Large bias",
  c("A-TMLE", "A-TMLE*", "ES-CVTMLE", "TMLE", "RCT ONLY"),
  list(c("A-TMLE", "RCT ONLY"), c("A-TMLE*", "RCT ONLY"), c("ES-CVTMLE", "RCT ONLY"), c("TMLE", "RCT ONLY")),
  atmle_both_res, atmle_tmle_res, escvtmle_res, tmle_res, rct_only_res)
ggsave(filename = "param_complex_bias_relative.pdf", plot = plt_relative_mse_param_complex_bias, device = "pdf",
       path = "plot", width = 10, height = 8, dpi = 300)

plt_param_complex_bias_prop <- get_plot_selected(escvtmle_res, "Complex parametric bias")
ggsave(filename = "param_complex_bias_prop.pdf", plot = plt_param_complex_bias_prop, device = "pdf",
       path = "plot", width = 10, height = 8, dpi = 300)

# combine all oracle relative MSE plots
# relative_mse_plt_oracle <- ggarrange(plt_no_bias_oracle_relative,
#                                      plt_param_simple_bias_oracle_relative,
#                                      plt_param_complex_oracle_relative,
#                                      nrow = 2, ncol = 2, common.legend = TRUE)
# ggsave(filename = "relative_mse_oracle.pdf", plot = relative_mse_plt_oracle, device = "pdf",
#        path = "plot", width = 10, height = 8, dpi = 300)

# combine all comparisons relative MSE plots
relative_mse_plt <- ggarrange(plt_relative_mse_no_bias,
                              plt_relative_mse_param_simple_bias,
                              plt_relative_mse_param_complex_bias,
                              nrow = 1, ncol = 3, common.legend = TRUE)
ggsave(filename = "relative_mse_all.pdf", plot = relative_mse_plt, device = "pdf",
       path = "plot", width = 12, height = 4, dpi = 300)

# e. HAL bias ------------------------------------------------------------------
# all comparisons
# atmle_both_res <- readRDS("out/atmle_both_HAL_glm_hal_1016.RDS")
# escvtmle_res <- readRDS("out/escvtmle_HAL_glm_hal_1016.RDS")
# tmle_res <- readRDS("out/tmle_HAL_glm_hal_1016.RDS")
# rct_only_res <- readRDS("out/rct_only_HAL_glm_hal_1016.RDS")
# plt_HAL_bias <- get_plot(
#   "Complex bias",
#   c("A-TMLE", "ES-CVTMLE", "TMLE", "RCT ONLY"),
#   atmle_both_res, escvtmle_res, tmle_res, rct_only_res)
# ggsave(filename = "HAL_bias.pdf", plot = plt_HAL_bias, device = "pdf",
#        path = "plot", width = 10, height = 8, dpi = 300)
#
# # all comparisons relative MSE
# plt_relative_mse_HAL_bias <- get_relative_mse_plot(
#   "Complex bias",
#   c("A-TMLE", "ES-CVTMLE", "TMLE", "RCT ONLY"),
#   list(c("A-TMLE", "RCT ONLY"), c("ES-CVTMLE", "RCT ONLY"), c("TMLE", "RCT ONLY")),
#   atmle_both_res, escvtmle_res, tmle_res, rct_only_res)
# ggsave(filename = "HAL_bias_relative.pdf", plot = plt_relative_mse_HAL_bias, device = "pdf",
#        path = "plot", width = 10, height = 8, dpi = 300)
#
# plt_HAL_bias_prop <- get_plot_selected(escvtmle_res, "HAL bias")
