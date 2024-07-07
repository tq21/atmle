library(purrr)
library(ggplot2)
library(ggpubr)

n_rct_seq <- c(400, 500, 600, 700, 800)
n_rwd_seq <- 3*n_rct_seq
n_seq <- n_rct_seq + n_rwd_seq
ate <- 1.5

# bias a -----------------------------------------------------------------------
load("out/bias_a.RData")
res_list <- list(all_naive_atmle, all_atmle, all_tmle)
names(res_list) <- c("Naive A-TMLE", "A-TMLE", "TMLE")
df_bias_a <- map_dfr(names(res_list), function(name) {
  map_dfr(1:length(n_seq), function(.i) {
    .x <- res_list[[name]][[.i]]
    return(data.frame(Estimator = name,
                      n = n_seq[.i],
                      bias = mean(abs(.x - ate)),
                      var = var(.x),
                      mse = mean((.x - ate)^2)))
  })
})

plt_a <- get_mse_plot(df_bias_a)

# bias a -----------------------------------------------------------------------
load("out/bias_b.RData")
res_list <- list(all_naive_atmle, all_atmle, all_tmle)
names(res_list) <- c("Naive A-TMLE", "A-TMLE", "TMLE")
df_bias_b <- map_dfr(names(res_list), function(name) {
  map_dfr(1:length(n_seq), function(.i) {
    .x <- res_list[[name]][[.i]]
    return(data.frame(Estimator = name,
                      n = n_seq[.i],
                      bias = mean(abs(.x - ate)),
                      var = var(.x),
                      mse = mean((.x - ate)^2)))
  })
})

plt_b <- get_mse_plot(df_bias_b)

plt <- ggarrange(plt_a, plt_b, nrow = 1, ncol = 2, common.legend = TRUE)
ggsave(filename = "plt.pdf", plot = plt, device = "pdf",
       path = "figs", width = 8, height = 4, dpi = 300)
