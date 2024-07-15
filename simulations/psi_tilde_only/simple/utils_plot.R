library(ggplot2)
library(ggpubr)

get_bias_plot <- function(df) {
  p_bias <- ggplot(df, aes(x = n, y = bias, color = Estimator)) +
    geom_point(size = 1.5) +
    geom_line(linewidth = 1) +
    labs(title = "",
         x = "n",
         y = "Absolute bias") +
    theme_minimal() +
    theme(text = element_text(size = 16),
          legend.position = "none")

  return(p_bias)
}

get_var_plot <- function(df) {
  p_var <- ggplot(df, aes(x = n, y = var, color = Estimator)) +
    geom_point(size = 1.5) +
    geom_line(linewidth = 1) +
    labs(title = "",
         x = "n",
         y = "Variance") +
    theme_minimal() +
    theme(text = element_text(size = 16),
          legend.position = "none")

  return(p_var)
}

get_mse_plot <- function(df) {
  p_mse <- ggplot(df, aes(x = n, y = mse, color = Estimator)) +
    geom_point(size = 1.5) +
    geom_line(linewidth = 1) +
    labs(title = "",
         x = "n",
         y = "ATE MSE") +
    theme_minimal() +
    theme(text = element_text(size = 16),
          legend.position = "none")

  return(p_mse)
}

get_cover_plot <- function(df) {
  p_cover <- ggplot(df, aes(x = n, y = cover, color = Estimator)) +
    geom_point(size = 1.5) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = 0.95, color = "red", linetype = "dashed", linewidth = 1) +
    scale_y_continuous(breaks = seq(0.8, 1.0, 0.05), limits = c(0.8, 1.0)) +
    labs(title = "",
         x = "n",
         y = "Coverage") +
    theme_minimal() +
    theme(text = element_text(size = 16),
          legend.position = "none")

  return(p_cover)
}

get_relative_mse_plot <- function(df) {
  relative_mse_plot <- ggplot(df, aes(x = n, y = relative_mse, color = Estimator)) +
    geom_point(size = 1.5) +
    geom_line(linewidth = 1) +
    labs(title = "",
         x = "n",
         y = "Relative MSE") +
    theme_minimal() +
    theme(text = element_text(size = 16),
          legend.position = "none")

  return(relative_mse_plot)
}

get_cate_mse_plot <- function(df) {
  cate_mse_plot <- ggplot(df, aes(x = n, y = cate_mse, color = Estimator)) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.5) +
    labs(title = "",
         x = "n",
         y = "CATE MSE") +
    theme_minimal() +
    theme(text = element_text(size = 16),
          legend.position = "none")

  return(cate_mse_plot)
}

plt_fun <- function() {
  truth <- get_truth(t0)

  # tmp quick result table
  df_r_learner <- map_dfr(1:length(n_seq), function(.i) {
    .x <- res_list$all_r_learner[[.i]]
    .y <- res_list$all_r_learner_cover[[.i]]
    return(data.frame(Estimator = "IPCW-R-loss",
                      n = n_seq[.i],
                      bias = mean(abs(.x - truth)),
                      var = var(.x),
                      mse = mean((.x - truth)^2),
                      cover = mean(.y)))
  })

  df_r_learner_tmle_lambda <- map_dfr(1:length(n_seq), function(.i) {
    .x <- res_list$all_r_learner_tmle_lambda[[.i]]
    .y <- res_list$all_r_learner_tmle_lambda_cover[[.i]]
    return(data.frame(Estimator = "IPCW-R-loss + TMLE hazard",
                      n = n_seq[.i],
                      bias = mean(abs(.x - truth)),
                      var = var(.x),
                      mse = mean((.x - truth)^2),
                      cover = mean(.y)))
  })

  df_r_learner_tmle_proj <- map_dfr(1:length(n_seq), function(.i) {
    .x <- res_list$all_r_learner_tmle_proj[[.i]]
    .y <- res_list$all_r_learner_tmle_proj_cover[[.i]]
    return(data.frame(Estimator = "IPCW-R-loss + TMLE hazard + proj. param.",
                      n = n_seq[.i],
                      bias = mean(abs(.x - truth)),
                      var = var(.x),
                      mse = mean((.x - truth)^2),
                      cover = mean(.y)))
  })

  # cate mse
  df_ipcw_r_loss_cate_mse <- map_dfr(1:length(n_seq), function(.i) {
    .x <- res_list$all_r_learner_cate_mse[[.i]]
    return(data.frame(Estimator = "IPCW R-learner",
                      n = n_seq[.i],
                      cate_mse = mean(.x)))
  })

  df_tmle_cate_mse <- map_dfr(1:length(n_seq), function(.i) {
    .x <- res_list$all_r_learner_tmle_proj_cate_mse[[.i]]
    return(data.frame(Estimator = "IPCW-R-loss + TMLE hazard + proj. param.",
                      n = n_seq[.i],
                      cate_mse = mean(.x)))
  })

  df_list <- list(df_r_learner, df_r_learner_tmle_lambda, df_r_learner_tmle_proj)
  df <- do.call(rbind, df_list)
  df$base_mse <- subset(df, Estimator == "IPCW-R-loss")$mse
  df$relative_mse <- df$base_mse / df$mse
  df_r_learner_tmle_lambda_cate_mse <- df_ipcw_r_loss_cate_mse
  df_r_learner_tmle_lambda_cate_mse$cate_mse <- NA
  df_r_learner_tmle_lambda_cate_mse$Estimator <- "IPCW-R-loss + TMLE hazard"
  df_cate <- rbind(df_ipcw_r_loss_cate_mse, df_r_learner_tmle_lambda_cate_mse, df_tmle_cate_mse)

  # plot
  mse_plt <- get_mse_plot(df)
  relative_mse_plt <- get_relative_mse_plot(df)
  cover_plt <- get_cover_plot(df)
  cate_mse_plt <- get_cate_mse_plot(df_cate)

  return(list(mse_plt,
              relative_mse_plt,
              cover_plt,
              cate_mse_plt))
}
