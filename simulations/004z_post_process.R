library(ggpubr)
library(purrr)
source("utils_plot.R")

# simulation parameters --------------------------------------------------------
B <- 500 # number of runs for each sample size n
n_min <- 500 # smallest sample size
n_max <- 5000 # largest sample size
n_step <- 500 # sample size increment
bA <- 1.5

# a. no bias -------------------------------------------------------------------
atmle_both_res <- readRDS("out/atmle_both_0_bias_glm_lasso_20230916_084333.RDS")
atmle_tmle_res <- readRDS("out/atmle_tmle_0_bias_glm_lasso_20230916_031404.RDS")
escvtmle_res <- readRDS("out/escvtmle_0_bias_glm_lasso_20230916_064412.RDS")
rct_only_res <- readRDS("out/rct_only_0_bias_glm_lasso_20230916_071239.RDS")
plt_no_bias <- get_plot(atmle_both_res, atmle_tmle_res, escvtmle_res, rct_only_res, "No bias")
plt_relative_mse_no_bias <- get_relative_mse_plot(atmle_both_res, atmle_tmle_res, escvtmle_res, "No bias")
plt_no_bias_prop <- get_plot_selected(escvtmle_res, "No bias")
# ggsave(filename = "no_bias.pdf", plot = plt_no_bias, device = "pdf",
#        path = "plot", width = 10, height = 8, dpi = 300)

# b. constant bias -------------------------------------------------------------
atmle_both_res <- readRDS("out/atmle_both_1.8_bias_glm_lasso_20230915_221142.RDS")
atmle_tmle_res <- readRDS("out/atmle_tmle_1.8_bias_glm_lasso_20230916_025210.RDS")
escvtmle_res <- readRDS("out/escvtmle_1.8_bias_glm_lasso_20230916_062259.RDS")
rct_only_res <- readRDS("out/rct_only_1.8_bias_glm_lasso_20230916_065130.RDS")
plt_constant_bias <- get_plot(atmle_both_res, atmle_tmle_res, escvtmle_res, rct_only_res, "Constant bias")
plt_relative_mse_constant_bias <- get_relative_mse_plot(atmle_both_res, atmle_tmle_res, escvtmle_res, "Constant bias")
plt_constant_bias_prop <- get_plot_selected(escvtmle_res, "Constant bias")
# ggsave(filename = "constant_bias.pdf", plot = plt_constant_bias, device = "pdf",
#        path = "plot", width = 10, height = 8, dpi = 300)

# c. simple parametric bias -------------------------------------------------------------
atmle_both_res <- readRDS("out/atmle_both_param_simple_glm_lasso_20230915_221730.RDS")
atmle_tmle_res <- readRDS("out/atmle_tmle_param_simple_glm_lasso_20230916_030246.RDS")
escvtmle_res <- readRDS("out/escvtmle_param_simple_glm_lasso_20230916_063327.RDS")
rct_only_res <- readRDS("out/rct_only_param_simple_glm_lasso_20230916_070123.RDS")
plt_param_simple_bias <- get_plot(atmle_both_res, atmle_tmle_res, escvtmle_res, rct_only_res, "Simple parametric bias")
plt_relative_mse_param_simple_bias <- get_relative_mse_plot(atmle_both_res, atmle_tmle_res, escvtmle_res, "Simple parametric bias")
plt_param_simple_bias_prop <- get_plot_selected(escvtmle_res, "Simple parametric bias")
# ggsave(filename = "param_simple_bias.pdf", plot = plt_param_simple_bias, device = "pdf",
#        path = "plot", width = 10, height = 8, dpi = 300)

# d. complex parametric bias -------------------------------------------------------------
atmle_both_res <- readRDS("out/atmle_both_param_complex_glm_lasso_20230913_143429.RDS")
atmle_tmle_res <- readRDS("out/atmle_tmle_param_complex_glm_lasso_913.RDS")
escvtmle_res <- readRDS("out/escvtmle_param_complex_glm_lasso_20230915_053026.RDS")
plt_param_complex_bias <- get_plot(atmle_both_res, atmle_tmle_res, escvtmle_res, "Complex parametric bias")
plt_relative_mse_param_complex_bias <- get_relative_mse_plot(atmle_both_res, atmle_tmle_res, escvtmle_res, "Complex parametric bias")
plt_param_complex_bias_prop <- get_plot_selected(escvtmle_res, "Complex parametric bias")
# ggsave(filename = "param_complex_bias.pdf", plot = plt_param_complex_bias, device = "pdf",
#        path = "plot", width = 10, height = 8, dpi = 300)

# e. misspecified bias ---------------------------------------------------------
# atmle_both_res <- readRDS("out/atmle_both_misspecify_glm_lasso_20230907_152456.RDS")
# atmle_tmle_res <- readRDS("out/atmle_tmle_misspecify_glm_lasso_20230907_172917.RDS")
# escvtmle_res <- readRDS("out/escvtmle_misspecify_glm_lasso_20230907_204604.RDS")
# get_plot(atmle_both_res, atmle_tmle_res, escvtmle_res)
# ggsave(filename = "misspecified_bias.pdf", plot = plt_param_complex_bias, device = "pdf",
#        path = "plot", width = 10, height = 8, dpi = 300)

# prop selected
prop_plt <- ggarrange(plt_no_bias_prop,
                      plt_constant_bias_prop,
                      plt_param_simple_bias_prop,
                      plt_param_complex_bias_prop,
                      nrow = 2, ncol = 2)

# mse ratios
relative_mse_plt <- ggarrange(plt_relative_mse_no_bias,
                              plt_relative_mse_constant_bias,
                              plt_relative_mse_param_simple_bias,
                              plt_relative_mse_param_complex_bias,
                              nrow = 2, ncol = 2, common.legend = TRUE)

ggsave(filename = "prop.pdf", plot = prop_plt, device = "pdf",
       path = "plot", width = 10, height = 8, dpi = 300)
