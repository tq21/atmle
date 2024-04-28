library(ggpubr)
library(purrr)
source("utils_plot.R")

# simulation parameters --------------------------------------------------------
B <- 500
n_rct_seq <- c(400, 500, 600, 700, 800)
n_rwd_seq <- n_rct_seq * 3
total_sample_sizes <- n_rct_seq + n_rwd_seq
ate <- 1.5

atmle_both_res <- readRDS("out/atmle_both_surrogate_0421.RDS")
escvtmle_res <- readRDS("out/escvtmle_surrogate_0421.RDS")
tmle_res <- readRDS("out/tmle_surrogate_0421.RDS")
rct_only_res <- readRDS("out/rct_only_surrogate_0421.RDS")

# MSE
plt_mse <- get_mse_plot("Differential Outcomes",
                        c("A-TMLE", "ES-CVTMLE", "TMLE"),
                        atmle_both_res, escvtmle_res, tmle_res)

# relative MSE
plt_relative_mse <- get_relative_mse_plot(
  "",
  c("A-TMLE", "ES-CVTMLE", "TMLE", "RCT ONLY"),
  list(c("A-TMLE", "RCT ONLY"), c("ES-CVTMLE", "RCT ONLY"), c("TMLE", "RCT ONLY")),
  atmle_both_res, escvtmle_res, tmle_res, rct_only_res)

# coverage
plt_coverage <- get_cover_plot(
  "",
  c("A-TMLE", "ES-CVTMLE", "TMLE", "RCT ONLY"),
  atmle_both_res, escvtmle_res, tmle_res, rct_only_res)

# proportion of selected
plt_prop <- get_plot_prop_selected(escvtmle_res, "")

plt <- ggarrange(plt_mse, plt_relative_mse, plt_coverage, plt_prop,
                 nrow = 1, ncol = 4, common.legend = TRUE)
ggsave(filename = "param.pdf", plot = plt, device = "pdf",
       path = "figs", width = 16, height = 4, dpi = 300)
