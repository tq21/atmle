library(ggpubr)
library(purrr)
source("utils_plot.R")

# simulation parameters --------------------------------------------------------
B <- 500 # number of runs for each sample size n
n_min <- 1000 # smallest sample size
n_max <- 3000 # largest sample size
n_step <- 500 # sample size increment
bA <- 1.5

# a. no bias -------------------------------------------------------------------
# oracle comparisons
atmle_both_res <- readRDS("out/atmle_both_0_bias_glm_glmnet_1009.RDS")
oracle_atmle_both_res <- readRDS("out/oracle_atmle_both_0_bias_glm_glmnet_1009.RDS")
atmle_tmle_res <- readRDS("out/atmle_tmle_0_bias_glm_glmnet_1009.RDS")
oracle_atmle_tmle_res <- readRDS("out/oracle_atmle_tmle_0_bias_glm_glmnet_1009.RDS")
plt_no_bias_oracle <- get_plot(
  "No bias",
  c("A-TMLE", "Oracle A-TMLE", "A-TMLE*", "Oracle A-TMLE*"),
  atmle_both_res, oracle_atmle_both_res, atmle_tmle_res, oracle_atmle_tmle_res)
ggsave(filename = "no_bias_oracle.pdf", plot = plt_no_bias_oracle, device = "pdf",
       path = "plot", width = 10, height = 8, dpi = 300)

# oracle relative MSE
plt_no_bias_oracle_relative <- get_relative_mse_plot(
  "No bias",
  c("A-TMLE", "Oracle A-TMLE", "A-TMLE*", "Oracle A-TMLE*"),
  list(c("A-TMLE", "Oracle A-TMLE"), c("A-TMLE*", "Oracle A-TMLE*")),
  atmle_both_res, oracle_atmle_both_res, atmle_tmle_res, oracle_atmle_tmle_res)



escvtmle_res <- readRDS("out/escvtmle_0_bias_glm_glmnet_1007.RDS")
rct_only_res <- readRDS("out/rct_only_0_bias_glm_glmnet_1007.RDS")
tmle_res <- readRDS("out/tmle_0_bias_glm_glmnet_1007.RDS")
plt_no_bias <- get_plot(atmle_both_res, atmle_tmle_res, escvtmle_res, tmle_res, "No bias")
plt_relative_mse_no_bias <- get_relative_mse_plot(atmle_both_res, atmle_tmle_res, escvtmle_res, rct_only_res, "No bias")
plt_no_bias_prop <- get_plot_selected(escvtmle_res, "No bias")
ggsave(filename = "no_bias.pdf", plot = plt_no_bias, device = "pdf",
       path = "plot", width = 10, height = 8, dpi = 300)

# c. simple parametric bias -------------------------------------------------------------
# oracle comparison
atmle_both_res <- readRDS("out/atmle_both_param_simple_glm_glmnet_1009.RDS")
oracle_atmle_both_res <- readRDS("out/oracle_atmle_both_param_simple_glm_glmnet_1009.RDS")
atmle_tmle_res <- readRDS("out/atmle_tmle_param_simple_glm_glmnet_1009.RDS")
oracle_atmle_tmle_res <- readRDS("out/oracle_atmle_tmle_param_simple_glm_glmnet_1009.RDS")
plt_param_simple_bias_oracle <- get_plot(
  "Parametric bias 1",
  c("A-TMLE", "Oracle A-TMLE", "A-TMLE*", "Oracle A-TMLE*"),
  atmle_both_res, oracle_atmle_both_res, atmle_tmle_res, oracle_atmle_tmle_res)
ggsave(filename = "param_simple_oracle.pdf", plot = plt_param_simple_bias_oracle, device = "pdf",
       path = "plot", width = 10, height = 8, dpi = 300)

# oracle relative MSE
plt_param_simple_oracle_relative <- get_relative_mse_plot(
  "Parametric bias 1",
  c("A-TMLE", "Oracle A-TMLE", "A-TMLE*", "Oracle A-TMLE*"),
  list(c("A-TMLE", "Oracle A-TMLE"), c("A-TMLE*", "Oracle A-TMLE*")),
  atmle_both_res, oracle_atmle_both_res, atmle_tmle_res, oracle_atmle_tmle_res)

plt_param_simple_bias <- get_plot(atmle_both_res, atmle_tmle_res, escvtmle_res, tmle_res, "Simple parametric bias")
plt_relative_mse_param_simple_bias <- get_relative_mse_plot(atmle_both_res, atmle_tmle_res, escvtmle_res, tmle_res, "Simple parametric bias")
plt_param_simple_bias_prop <- get_plot_selected(escvtmle_res, "Simple parametric bias")
ggsave(filename = "param_simple_bias.pdf", plot = plt_param_simple_bias, device = "pdf",
       path = "plot", width = 10, height = 8, dpi = 300)

# d. complex parametric bias -------------------------------------------------------------
# oracle comparison
atmle_both_res <- readRDS("out/atmle_both_param_complex_glm_glmnet_1009.RDS")
oracle_atmle_both_res <- readRDS("out/oracle_atmle_both_param_complex_glm_glmnet_1009.RDS")
atmle_tmle_res <- readRDS("out/atmle_tmle_param_complex_glm_glmnet_1009.RDS")
oracle_atmle_tmle_res <- readRDS("out/oracle_atmle_tmle_param_complex_glm_glmnet_1009.RDS")
plt_param_complex_bias_oracle <- get_plot(
  "Parametric bias 2",
  c("A-TMLE", "Oracle A-TMLE", "A-TMLE*", "Oracle A-TMLE*"),
  atmle_both_res, oracle_atmle_both_res, atmle_tmle_res, oracle_atmle_tmle_res)
ggsave(filename = "param_complex_oracle.pdf", plot = plt_param_complex_bias_oracle, device = "pdf",
       path = "plot", width = 10, height = 8, dpi = 300)

# oracle relative MSE
plt_param_complex_oracle_relative <- get_relative_mse_plot(
  "Parametric bias 2",
  c("A-TMLE", "Oracle A-TMLE", "A-TMLE*", "Oracle A-TMLE*"),
  list(c("A-TMLE", "Oracle A-TMLE"), c("A-TMLE*", "Oracle A-TMLE*")),
  atmle_both_res, oracle_atmle_both_res, atmle_tmle_res, oracle_atmle_tmle_res)


atmle_both_res <- readRDS("out/atmle_both_param_complex_glm_glmnet_1007.RDS")
atmle_tmle_res <- readRDS("out/atmle_tmle_param_complex_glm_glmnet_1007.RDS")
escvtmle_res <- readRDS("out/escvtmle_param_complex_glm_glmnet_1007.RDS")
rct_only_res <- readRDS("out/rct_only_param_complex_glm_glmnet_1007.RDS")
tmle_res <- readRDS("out/tmle_param_complex_glm_glmnet_1007.RDS")
plt_param_complex_bias <- get_plot(atmle_both_res, atmle_tmle_res, escvtmle_res, tmle_res, "Parametric bias")
plt_relative_mse_param_complex_bias <- get_relative_mse_plot(atmle_both_res, atmle_tmle_res, escvtmle_res, rct_only_res, "Complex parametric bias")
plt_param_complex_bias_prop <- get_plot_selected(escvtmle_res, "Complex parametric bias")
ggsave(filename = "param_complex_bias.pdf", plot = plt_param_complex_bias, device = "pdf",
       path = "plot", width = 10, height = 8, dpi = 300)

# e. HAL bias ------------------------------------------------------------------
atmle_both_res <- readRDS("out/atmle_both_HAL_glm_HAL_1008.RDS")
escvtmle_res <- readRDS("out/escvtmle_HAL_glm_HAL_1008.RDS")
rct_only_res <- readRDS("out/rct_only_HAL_glm_HAL_1008.RDS")
tmle_res <- readRDS("out/tmle_HAL_glm_HAL_1008.RDS")
plt_HAL_bias <- get_plot("HAL bias",
                         c("A-TMLE", "ESCVTMLE", "TMLE", "RCT ONLY"),
                         atmle_both_res, escvtmle_res, tmle_res, rct_only_res)


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
ggsave(filename = "relative_mse_200_rct.pdf", plot = relative_mse_plt, device = "pdf",
       path = "plot", width = 10, height = 8, dpi = 300)

ggsave(filename = "prop_200_rct.pdf", plot = prop_plt, device = "pdf",
       path = "plot", width = 10, height = 8, dpi = 300)
