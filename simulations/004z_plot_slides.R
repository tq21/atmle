library(ggpubr)
library(purrr)
library(grid)
source("utils_plot.R")

# simulation parameters --------------------------------------------------------
B <- 1000 # number of runs for each sample size n
n_min <- 1000 # smallest sample size
n_max <- 3000 # largest sample size
n_step <- 500 # sample size increment
bA <- 1.5
names <- c("A-TMLE", "ES-CVTMLE", "TMLE")

# no bias ----------------------------------------------------------------------
no_bias_atmle_res <- readRDS("out/atmle_both_0_bias_glm_glmnet_1028.RDS")
no_bias_escvtmle_res <- readRDS("out/escvtmle_0_bias_glm_glmnet_1028.RDS")
no_bias_tmle_res <- readRDS("out/tmle_0_bias_glm_glmnet_1028.RDS")
no_bias_rct_only_res <- readRDS("out/rct_only_0_bias_glm_glmnet_1028.RDS")
no_bias_res_list <- list(no_bias_atmle_res, no_bias_escvtmle_res, no_bias_tmle_res)
no_bias_res_dt <- map_dfr(1:length(no_bias_res_list), function(i) {
  return(get_res(no_bias_res_list[[i]], bA, names[i]))
})
no_bias_res_dt <- no_bias_res_dt[order(no_bias_res_dt$n), ]
no_bias_res_dt$scenario <- "No bias"

# small bias -------------------------------------------------------------------
small_bias_atmle_res <- readRDS("out/atmle_both_param_simple_glm_glmnet_1028.RDS")
small_bias_escvtmle_res <- readRDS("out/escvtmle_param_simple_glm_glmnet_1028.RDS")
small_bias_tmle_res <- readRDS("out/tmle_param_simple_glm_glmnet_1028.RDS")
small_bias_rct_only_res <- readRDS("out/rct_only_param_simple_glm_glmnet_1028.RDS")
small_bias_res_list <- list(small_bias_atmle_res, small_bias_escvtmle_res, small_bias_tmle_res)
small_bias_res_dt <- map_dfr(1:length(small_bias_res_list), function(i) {
  return(get_res(small_bias_res_list[[i]], bA, names[i]))
})
small_bias_res_dt <- small_bias_res_dt[order(small_bias_res_dt$n), ]
small_bias_res_dt$scenario <- "Small bias"

# large bias -------------------------------------------------------------------
large_bias_atmle_res <- readRDS("out/atmle_both_param_complex_glm_glmnet_1028.RDS")
large_bias_escvtmle_res <- readRDS("out/escvtmle_param_complex_glm_glmnet_1028.RDS")
large_bias_tmle_res <- readRDS("out/tmle_param_complex_glm_glmnet_1028.RDS")
large_bias_rct_only_res <- readRDS("out/rct_only_param_complex_glm_glmnet_1028.RDS")
large_bias_res_list <- list(large_bias_atmle_res, large_bias_escvtmle_res, large_bias_tmle_res)
large_bias_res_dt <- map_dfr(1:length(large_bias_res_list), function(i) {
  return(get_res(large_bias_res_list[[i]], bA, names[i]))
})
large_bias_res_dt <- large_bias_res_dt[order(large_bias_res_dt$n), ]
large_bias_res_dt$scenario <- "Large bias"

# plot 1: MSE comparisons ------------------------------------------------------
no_bias_mse_plt <- ggplot(no_bias_res_dt, aes(x = n, y = mse, color = Estimator)) +
  geom_point(size = 1.5) +
  geom_line(linewidth = 1.5) +
  labs(title = "No bias",
       x = "n",
       y = "MSE") +
  theme_minimal() +
  theme(text = element_text(size = 20),
        plot.title = element_text(face = "bold", size = 20, hjust = 0.5))

small_bias_mse_plt <- ggplot(small_bias_res_dt, aes(x = n, y = mse, color = Estimator)) +
  geom_point(size = 1.5) +
  geom_line(linewidth = 1.5) +
  labs(title = "Small bias",
       x = "n",
       y = "MSE") +
  theme_minimal() +
  theme(text = element_text(size = 20),
        plot.title = element_text(face = "bold", size = 20, hjust = 0.5))

large_bias_mse_plt <- ggplot(large_bias_res_dt, aes(x = n, y = mse, color = Estimator)) +
  geom_point(size = 1.5) +
  geom_line(linewidth = 1.5) +
  labs(title = "Large bias",
       x = "n",
       y = "MSE") +
  theme_minimal() +
  theme(text = element_text(size = 20),
        plot.title = element_text(face = "bold", size = 20, hjust = 0.5))

mse_plt <- ggarrange(no_bias_mse_plt, small_bias_mse_plt, large_bias_mse_plt,
                     nrow = 1, ncol = 3, common.legend = TRUE)
mse_plt <- annotate_figure(mse_plt, top = text_grob("MSE comparisons\n", face = "bold", size = 30))
ggsave(filename = "mse.pdf", plot = mse_plt, device = "pdf",
       path = "plot", width = 18, height = 6, dpi = 300)

# plot 2: relative efficiency comparisons --------------------------------------
no_bias_relative_plt <- get_relative_mse_plot(
  "No bias",
  c("A-TMLE", "ES-CVTMLE", "TMLE", "RCT ONLY"),
  list(c("A-TMLE", "RCT ONLY"), c("ES-CVTMLE", "RCT ONLY"), c("TMLE", "RCT ONLY")),
  no_bias_atmle_res, no_bias_escvtmle_res, no_bias_tmle_res, no_bias_rct_only_res)

small_bias_relative_plt <- get_relative_mse_plot(
  "Small bias",
  c("A-TMLE", "ES-CVTMLE", "TMLE", "RCT ONLY"),
  list(c("A-TMLE", "RCT ONLY"), c("ES-CVTMLE", "RCT ONLY"), c("TMLE", "RCT ONLY")),
  small_bias_atmle_res, small_bias_escvtmle_res, small_bias_tmle_res, small_bias_rct_only_res)

large_bias_relative_plt <- get_relative_mse_plot(
  "Large bias",
  c("A-TMLE", "ES-CVTMLE", "TMLE", "RCT ONLY"),
  list(c("A-TMLE", "RCT ONLY"), c("ES-CVTMLE", "RCT ONLY"), c("TMLE", "RCT ONLY")),
  large_bias_atmle_res, large_bias_escvtmle_res, large_bias_tmle_res, large_bias_rct_only_res)

relative_plt <- ggarrange(no_bias_relative_plt, small_bias_relative_plt, large_bias_relative_plt,
                          nrow = 1, ncol = 3, common.legend = TRUE)
relative_plt <- annotate_figure(relative_plt, top = text_grob("Relative efficiency comparisons\n", face = "bold", size = 30))
ggsave(filename = "relative.pdf", plot = relative_plt, device = "pdf",
       path = "plot", width = 18, height = 6, dpi = 300)

# plot 3: coverage comparisons -------------------------------------------------
no_bias_cover_plt <- ggplot(no_bias_res_dt, aes(x = n, y = cover, color = Estimator)) +
  geom_point(size = 1.5) +
  geom_line(linewidth = 1.5) +
  geom_hline(yintercept = 0.95, color = "red", linetype = "dashed", linewidth = 1) +
  scale_y_continuous(breaks = seq(0.8, 1.0, 0.05), limits = c(0.8, 1.0)) +
  labs(title = "No bias",
       x = "n",
       y = "Coverage") +
  theme_minimal() +
  theme(text = element_text(size = 20),
        plot.title = element_text(face = "bold", size = 20, hjust = 0.5))

small_bias_cover_plt <- ggplot(small_bias_res_dt, aes(x = n, y = cover, color = Estimator)) +
  geom_point(size = 1.5) +
  geom_line(linewidth = 1.5) +
  geom_hline(yintercept = 0.95, color = "red", linetype = "dashed", linewidth = 1) +
  scale_y_continuous(breaks = seq(0.8, 1.0, 0.05), limits = c(0.8, 1.0)) +
  labs(title = "Small bias",
       x = "n",
       y = "Coverage") +
  theme_minimal() +
  theme(text = element_text(size = 20),
        plot.title = element_text(face = "bold", size = 20, hjust = 0.5))

large_bias_cover_plt <- ggplot(large_bias_res_dt, aes(x = n, y = cover, color = Estimator)) +
  geom_point(size = 1.5) +
  geom_line(linewidth = 1.5) +
  geom_hline(yintercept = 0.95, color = "red", linetype = "dashed", linewidth = 1) +
  scale_y_continuous(breaks = seq(0.8, 1.0, 0.05), limits = c(0.8, 1.0)) +
  labs(title = "Large bias",
       x = "n",
       y = "Coverage") +
  theme_minimal() +
  theme(text = element_text(size = 20),
        plot.title = element_text(face = "bold", size = 20, hjust = 0.5))

cover_plt <- ggarrange(no_bias_cover_plt, small_bias_cover_plt, large_bias_cover_plt,
                       nrow = 1, ncol = 3, common.legend = TRUE)
cover_plt <- annotate_figure(cover_plt, top = text_grob("Coverage comparisons\n", face = "bold", size = 30))
ggsave(filename = "cover.pdf", plot = cover_plt, device = "pdf",
       path = "plot", width = 18, height = 6, dpi = 300)

B <- 500 # number of runs for each sample size n

# d. HAL 1 bias ----------------------------------------------------------------
hal_1_bias_atmle_res <- readRDS("out/atmle_both_HAL_1_glm_HAL_1030.RDS")
hal_1_bias_escvtmle_res <- readRDS("out/escvtmle_HAL_1_glm_HAL_1030.RDS")
hal_1_bias_tmle_res <- readRDS("out/tmle_HAL_1_glm_HAL_1030.RDS")
hal_1_bias_rct_only_res <- readRDS("out/rct_only_HAL_1_glm_HAL_1030.RDS")
hal_1_bias_res_list <- list(hal_1_bias_atmle_res, hal_1_bias_escvtmle_res, hal_1_bias_tmle_res)
hal_1_bias_res_dt <- map_dfr(1:length(hal_1_bias_res_list), function(i) {
  return(get_res(hal_1_bias_res_list[[i]], bA, names[i]))
})
hal_1_bias_res_dt <- hal_1_bias_res_dt[order(hal_1_bias_res_dt$n), ]
hal_1_bias_res_dt$scenario <- "HAL 1 bias"

# e. HAL 2 bias ----------------------------------------------------------------
hal_2_bias_atmle_res <- readRDS("out/atmle_both_HAL_2_glm_HAL_1030.RDS")
hal_2_bias_escvtmle_res <- readRDS("out/escvtmle_HAL_2_glm_HAL_1030.RDS")
hal_2_bias_tmle_res <- readRDS("out/tmle_HAL_2_glm_HAL_1030.RDS")
hal_2_bias_rct_only_res <- readRDS("out/rct_only_HAL_2_glm_HAL_1030.RDS")
hal_2_bias_res_list <- list(hal_2_bias_atmle_res, hal_2_bias_escvtmle_res, hal_2_bias_tmle_res)
hal_2_bias_res_dt <- map_dfr(1:length(hal_2_bias_res_list), function(i) {
  return(get_res(hal_2_bias_res_list[[i]], bA, names[i]))
})
hal_2_bias_res_dt <- hal_2_bias_res_dt[order(hal_2_bias_res_dt$n), ]
hal_2_bias_res_dt$scenario <- "HAL 2 bias"

# f. HAL 3 bias ----------------------------------------------------------------
hal_3_bias_atmle_res <- readRDS("out/atmle_both_HAL_3_glm_HAL_1030.RDS")
hal_3_bias_escvtmle_res <- readRDS("out/escvtmle_HAL_3_glm_HAL_1030.RDS")
hal_3_bias_tmle_res <- readRDS("out/tmle_HAL_3_glm_HAL_1030.RDS")
hal_3_bias_rct_only_res <- readRDS("out/rct_only_HAL_3_glm_HAL_1030.RDS")
hal_3_bias_res_list <- list(hal_3_bias_atmle_res, hal_3_bias_escvtmle_res, hal_3_bias_tmle_res)
hal_3_bias_res_dt <- map_dfr(1:length(hal_3_bias_res_list), function(i) {
  return(get_res(hal_3_bias_res_list[[i]], bA, names[i]))
})
hal_3_bias_res_dt <- hal_3_bias_res_dt[order(hal_3_bias_res_dt$n), ]
hal_3_bias_res_dt$scenario <- "HAL 3 bias"

# plot 4: MSE comparisons, HAL -------------------------------------------------
HAL_1_bias_mse_plt <- ggplot(hal_1_bias_res_dt, aes(x = n, y = mse, color = Estimator)) +
  geom_point(size = 1.5) +
  geom_line(linewidth = 1.5) +
  labs(title = "Complex bias 1",
       x = "n",
       y = "MSE") +
  theme_minimal() +
  theme(text = element_text(size = 20),
        plot.title = element_text(face = "bold", size = 20, hjust = 0.5))

HAL_2_bias_mse_plt <- ggplot(hal_2_bias_res_dt, aes(x = n, y = mse, color = Estimator)) +
  geom_point(size = 1.5) +
  geom_line(linewidth = 1.5) +
  labs(title = "Complex bias 2",
       x = "n",
       y = "MSE") +
  theme_minimal() +
  theme(text = element_text(size = 20),
        plot.title = element_text(face = "bold", size = 20, hjust = 0.5))

HAL_3_bias_mse_plt <- ggplot(hal_3_bias_res_dt, aes(x = n, y = mse, color = Estimator)) +
  geom_point(size = 1.5) +
  geom_line(linewidth = 1.5) +
  labs(title = "Complex bias 3",
       x = "n",
       y = "MSE") +
  theme_minimal() +
  theme(text = element_text(size = 20),
        plot.title = element_text(face = "bold", size = 20, hjust = 0.5))

HAL_mse_plt <- ggarrange(HAL_1_bias_mse_plt, HAL_2_bias_mse_plt, HAL_3_bias_mse_plt,
                         nrow = 1, ncol = 3, common.legend = TRUE)
HAL_mse_plt <- annotate_figure(HAL_mse_plt, top = text_grob("MSE comparisons\n", face = "bold", size = 30))
ggsave(filename = "HAL_mse.pdf", plot = HAL_mse_plt, device = "pdf",
       path = "plot", width = 18, height = 6, dpi = 300)

# plot 5: relative MSE comparisons, HAL ----------------------------------------
HAL_1_bias_relative_plt <- get_relative_mse_plot(
  "Complex bias 1",
  c("A-TMLE", "ES-CVTMLE", "TMLE", "RCT ONLY"),
  list(c("A-TMLE", "RCT ONLY"), c("ES-CVTMLE", "RCT ONLY"), c("TMLE", "RCT ONLY")),
  hal_1_bias_atmle_res, hal_1_bias_escvtmle_res, hal_1_bias_tmle_res, hal_1_bias_rct_only_res)

HAL_2_bias_relative_plt <- get_relative_mse_plot(
  "Complex bias 2",
  c("A-TMLE", "ES-CVTMLE", "TMLE", "RCT ONLY"),
  list(c("A-TMLE", "RCT ONLY"), c("ES-CVTMLE", "RCT ONLY"), c("TMLE", "RCT ONLY")),
  hal_2_bias_atmle_res, hal_2_bias_escvtmle_res, hal_2_bias_tmle_res, hal_2_bias_rct_only_res)

HAL_3_bias_relative_plt <- get_relative_mse_plot(
  "Complex bias 3",
  c("A-TMLE", "ES-CVTMLE", "TMLE", "RCT ONLY"),
  list(c("A-TMLE", "RCT ONLY"), c("ES-CVTMLE", "RCT ONLY"), c("TMLE", "RCT ONLY")),
  hal_3_bias_atmle_res, hal_3_bias_escvtmle_res, hal_3_bias_tmle_res, hal_3_bias_rct_only_res)

HAL_relative_plt <- ggarrange(HAL_1_bias_relative_plt, HAL_2_bias_relative_plt, HAL_3_bias_relative_plt,
                              nrow = 1, ncol = 3, common.legend = TRUE)
HAL_relative_plt <- annotate_figure(HAL_relative_plt, top = text_grob("Relative efficiency comparisons\n", face = "bold", size = 30))
ggsave(filename = "HAL_relative.pdf", plot = HAL_relative_plt, device = "pdf",
       path = "plot", width = 18, height = 6, dpi = 300)

# plot 5: Coverage comparisons, HAL --------------------------------------------
HAL_1_bias_cover_plt <- ggplot(hal_1_bias_res_dt, aes(x = n, y = cover, color = Estimator)) +
  geom_point(size = 1.5) +
  geom_line(linewidth = 1.5) +
  geom_hline(yintercept = 0.95, color = "red", linetype = "dashed", linewidth = 1) +
  scale_y_continuous(breaks = seq(0.8, 1.0, 0.05), limits = c(0.8, 1.0)) +
  labs(title = "Complex bias 1",
       x = "n",
       y = "Coverage") +
  theme_minimal() +
  theme(text = element_text(size = 20),
        plot.title = element_text(face = "bold", size = 20, hjust = 0.5))

HAL_2_bias_cover_plt <- ggplot(hal_2_bias_res_dt, aes(x = n, y = cover, color = Estimator)) +
  geom_point(size = 1.5) +
  geom_line(linewidth = 1.5) +
  geom_hline(yintercept = 0.95, color = "red", linetype = "dashed", linewidth = 1) +
  scale_y_continuous(breaks = seq(0.8, 1.0, 0.05), limits = c(0.8, 1.0)) +
  labs(title = "Complex bias 2",
       x = "n",
       y = "Coverage") +
  theme_minimal() +
  theme(text = element_text(size = 20),
        plot.title = element_text(face = "bold", size = 20, hjust = 0.5))

HAL_3_bias_cover_plt <- ggplot(hal_3_bias_res_dt, aes(x = n, y = cover, color = Estimator)) +
  geom_point(size = 1.5) +
  geom_line(linewidth = 1.5) +
  geom_hline(yintercept = 0.95, color = "red", linetype = "dashed", linewidth = 1) +
  scale_y_continuous(breaks = seq(0.8, 1.0, 0.05), limits = c(0.8, 1.0)) +
  labs(title = "Complex bias 3",
       x = "n",
       y = "Coverage") +
  theme_minimal() +
  theme(text = element_text(size = 20),
        plot.title = element_text(face = "bold", size = 20, hjust = 0.5))

HAL_cover_plt <- ggarrange(HAL_1_bias_cover_plt, HAL_2_bias_cover_plt, HAL_3_bias_cover_plt,
                           nrow = 1, ncol = 3, common.legend = TRUE)
HAL_cover_plt <- annotate_figure(HAL_cover_plt, top = text_grob("Coverage comparisons\n", face = "bold", size = 30))
ggsave(filename = "HAL_cover.pdf", plot = HAL_cover_plt, device = "pdf",
       path = "plot", width = 18, height = 6, dpi = 300)
