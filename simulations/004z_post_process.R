library(ggpubr)
library(purrr)

# simulation parameters --------------------------------------------------------
B <- 500 # number of runs for each sample size n
n_min <- 500 # smallest sample size
n_max <- 5000 # largest sample size
n_step <- 500 # sample size increment

get_res <- function(atmle_both_res, atmle_tmle_res, escvtmle_res) {
  bias_atmle_both <- unlist(map(atmle_both_res$all_psi_est, function(.x) mean(.x-bA)))
  bias_atmle_tmle <- unlist(map(atmle_tmle_res$all_psi_est, function(.x) mean(.x-bA)))
  bias_escvtmle <- unlist(map(escvtmle_res$all_psi_est, function(.x) mean(.x-bA)))

  se_atmle_both <- unlist(map(atmle_both_res$all_psi_est, function(.x) sd(.x)))
  se_atmle_tmle <- unlist(map(atmle_tmle_res$all_psi_est, function(.x) sd(.x)))
  se_escvtmle <- unlist(map(escvtmle_res$all_psi_est, function(.x) sd(.x)))

  var_atmle_both <- unlist(map(atmle_both_res$all_psi_est, function(.x) var(.x)))
  var_atmle_tmle <- unlist(map(atmle_tmle_res$all_psi_est, function(.x) var(.x)))
  var_escvtmle <- unlist(map(escvtmle_res$all_psi_est, function(.x) var(.x)))

  mse_atmle_both <- bias_atmle_both^2 + var_atmle_both
  mse_atmle_tmle <- bias_atmle_tmle^2 + var_atmle_tmle
  mse_escvtmle <- bias_escvtmle^2 + var_escvtmle

  cover_atmle_both <- unlist(map(atmle_both_res$all_psi_coverage, function(.x) mean(.x)))
  cover_atmle_tmle <- unlist(map(atmle_tmle_res$all_psi_coverage, function(.x) mean(.x)))
  cover_escvtmle <- unlist(map(escvtmle_res$all_psi_coverage, function(.x) mean(.x)))

  return(data.frame(n = seq(n_min, n_max, n_step),
                    estimator = rep(c("A-TMLE", "A-TMLE*", "ESCVTMLE"), each = length(bias_atmle_both)),
                    bias = c(bias_atmle_both, bias_atmle_tmle, bias_escvtmle),
                    se = c(se_atmle_both, se_atmle_tmle, se_escvtmle),
                    mse = c(mse_atmle_both, mse_atmle_tmle, mse_escvtmle),
                    coverage = c(cover_atmle_both, cover_atmle_tmle, cover_escvtmle)))
}


get_plot <- function(atmle_both_res, atmle_tmle_res, escvtmle_res, title) {
  dt_res <- get_res(atmle_both_res, atmle_tmle_res, escvtmle_res)
  dt_res <- dt_res[order(dt_res$n), ]

  p_bias <- ggplot(dt_res, aes(x = n, y = abs(bias), color = estimator)) +
    geom_point() +
    geom_line() +
    labs(title = "",
         x = "n",
         y = "bias") +
    theme_minimal() +
    theme(text = element_text(size = 16))

  p_se <- ggplot(dt_res, aes(x = n, y = se, color = estimator)) +
    geom_point() +
    geom_line() +
    labs(title = "",
         x = "n",
         y = "standard error") +
    theme_minimal() +
    theme(text = element_text(size = 16))

  p_mse <- ggplot(dt_res, aes(x = n, y = mse, color = estimator)) +
    geom_point() +
    geom_line() +
    labs(title = "",
         x = "n",
         y = "mse") +
    theme_minimal() +
    theme(text = element_text(size = 16))

  p_cover <- ggplot(dt_res, aes(x = n, y = coverage, color = estimator)) +
    geom_point() +
    geom_line() +
    geom_hline(yintercept = 0.95, color = "red", linetype = "dashed", linewidth = 1) +
    scale_y_continuous(breaks = seq(0.8, 1.0, 0.05), limits = c(0.8, 1.0)) +
    labs(title = "",
         x = "n",
         y = "coverage") +
    theme_minimal() +
    theme(text = element_text(size = 16))

  plt <- ggarrange(p_bias, p_se, p_mse, p_cover,
                   nrow = 2, ncol = 2, common.legend = TRUE)
  plt <- annotate_figure(plt, top = text_grob(title, face = "bold", size = 16))

  return(plt)
}

# a. no bias -------------------------------------------------------------------
atmle_both_res <- readRDS("out/atmle_both_0_bias_glm_lasso_20230907_152348.RDS")
atmle_tmle_res <- readRDS("out/atmle_tmle_0_bias_glm_lasso_20230907_200301.RDS")
escvtmle_res <- readRDS("out/escvtmle_0_bias_glm_lasso_20230907_232009.RDS")
plt_no_bias <- get_plot(atmle_both_res, atmle_tmle_res, escvtmle_res, "No bias")
ggsave(filename = "no_bias.pdf", plot = plt_no_bias, device = "pdf",
       path = "plot", width = 10, height = 8, dpi = 300)

# b. constant bias -------------------------------------------------------------
atmle_both_res <- readRDS("out/atmle_both_1.8_bias_glm_lasso_20230907_152229.RDS")
atmle_tmle_res <- readRDS("out/atmle_tmle_1.8_bias_glm_lasso_20230907_193918.RDS")
escvtmle_res <- readRDS("out/escvtmle_1.8_bias_glm_lasso_20230907_230126.RDS")
plt_constant_bias <- get_plot(atmle_both_res, atmle_tmle_res, escvtmle_res, "Constant bias")
ggsave(filename = "constant_bias.pdf", plot = plt_constant_bias, device = "pdf",
       path = "plot", width = 10, height = 8, dpi = 300)

# c. simple parametric bias -------------------------------------------------------------
atmle_both_res <- readRDS("out/atmle_both_param_simple_bias_glm_lasso_20230907_152342.RDS")
atmle_tmle_res <- readRDS("out/atmle_tmle_param_simple_glm_lasso_20230907_195015.RDS")
escvtmle_res <- readRDS("out/escvtmle_param_simple_glm_lasso_20230907_231001.RDS")
plt_param_simple_bias <- get_plot(atmle_both_res, atmle_tmle_res, escvtmle_res, "Simple parametric bias")
ggsave(filename = "param_simple_bias.pdf", plot = plt_param_simple_bias, device = "pdf",
       path = "plot", width = 10, height = 8, dpi = 300)

# d. complex parametric bias -------------------------------------------------------------
atmle_both_res <- readRDS("out/atmle_both_param_complex_glm_lasso_20230907_152214.RDS")
atmle_tmle_res <- readRDS("out/atmle_tmle_param_complex_glm_lasso_20230907_191426.RDS")
escvtmle_res <- readRDS("out/escvtmle_param_complex_glm_lasso_20230907_222628.RDS")
plt_param_complex_bias <- get_plot(atmle_both_res, atmle_tmle_res, escvtmle_res, "Complex parametric bias")
ggsave(filename = "param_complex_bias.pdf", plot = plt_param_complex_bias, device = "pdf",
       path = "plot", width = 10, height = 8, dpi = 300)

# e. misspecified bias ---------------------------------------------------------
atmle_both_res <- readRDS("out/atmle_both_misspecify_glm_lasso_20230907_152456.RDS")
atmle_tmle_res <- readRDS("out/atmle_tmle_misspecify_glm_lasso_20230907_172917.RDS")
escvtmle_res <- readRDS("out/escvtmle_misspecify_glm_lasso_20230907_204604.RDS")
get_plot(atmle_both_res, atmle_tmle_res, escvtmle_res)
ggsave(filename = "misspecified_bias.pdf", plot = plt_param_complex_bias, device = "pdf",
       path = "plot", width = 10, height = 8, dpi = 300)
