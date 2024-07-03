source("utils.R")
source("utils_plot.R")
source("sim_data_002.R")
res_list <- readRDS("out/res_list_002.RDS")
set.seed(123)

n_seq <- seq(500, 2000, 500)
t0 <- 3
truth <- get_truth(t0)

# tmp quick result table
df_r_learner <- map_dfr(1:length(n_seq), function(.i) {
  .x <- res_list$all_r_learner[[.i]]
  .y <- res_list$all_r_learner_cover[[.i]]
  return(data.frame(estimator = "IPCW R-learner",
                    n = n_seq[.i],
                    bias = mean(abs(.x - truth)),
                    var = var(.x),
                    mse = mean((.x - truth)^2),
                    cover = mean(.y)))
})

df_tmle <- map_dfr(1:length(n_seq), function(.i) {
  .x <- res_list$all_tmle[[.i]]
  .y <- res_list$all_tmle_cover[[.i]]
  return(data.frame(estimator = "TMLE survival",
                    n = n_seq[.i],
                    bias = mean(abs(.x - truth)),
                    var = var(.x),
                    mse = mean((.x - truth)^2),
                    cover = mean(.y)))
})

df_r_learner_tmle_lambda <- map_dfr(1:length(n_seq), function(.i) {
  .x <- res_list$all_r_learner_tmle_lambda[[.i]]
  .y <- res_list$all_r_learner_tmle_lambda_cover[[.i]]
  return(data.frame(estimator = "IPCW R-learner + TMLE lambda",
                    n = n_seq[.i],
                    bias = mean(abs(.x - truth)),
                    var = var(.x),
                    mse = mean((.x - truth)^2),
                    cover = mean(.y)))
})

df_r_learner_tmle_proj <- map_dfr(1:length(n_seq), function(.i) {
  .x <- res_list$all_r_learner_tmle_proj[[.i]]
  .y <- res_list$all_r_learner_tmle_proj_cover[[.i]]
  return(data.frame(estimator = "IPCW R-learner + TMLE projection",
                    n = n_seq[.i],
                    bias = mean(abs(.x - truth)),
                    var = var(.x),
                    mse = mean((.x - truth)^2),
                    cover = mean(.y)))
})

df_list <- list(df_r_learner, df_tmle, df_r_learner_tmle_lambda, df_r_learner_tmle_proj)
df <- do.call(rbind, df_list)
df$base_mse <- subset(df, estimator == "IPCW R-learner")$mse
df$relative_mse <- df$base_mse / df$mse

# plot
bias_plt <- get_bias_plot(df)
var_plt <- get_var_plot(df)
mse_plt <- get_mse_plot(df)
relative_mse_plt <- get_relative_mse_plot(df)
cover_plt <- get_cover_plot(df)
plt <- ggarrange(bias_plt, var_plt, mse_plt, relative_mse_plt, cover_plt,
                 nrow = 1, ncol = 5, common.legend = TRUE)
ggsave(filename = "simple_002.pdf", plot = plt, device = "pdf",
       path = "figs", width = 16, height = 4, dpi = 300)
