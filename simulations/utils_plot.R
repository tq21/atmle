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

get_plot_selected <- function(escvtmle_res, name) {
  df <- data.frame(prop = escvtmle_res$escvtmle_prop_selected)
  p <- ggplot(df, aes(x = prop)) +
    geom_bar(stat = "count") +
    scale_x_continuous(breaks = seq(0, 1, 0.2)) +
    labs(title = name,
         x = "prop. of folds selected",
         y = "frequency") +
    theme_minimal() +
    theme(text = element_text(size = 16),
          plot.title = element_text(hjust = 0.5))

  return(p)
}

get_relative_mse_plot <- function(atmle_both_res, atmle_tmle_res, escvtmle_res, name) {
  dt_res <- get_res(atmle_both_res, atmle_tmle_res, escvtmle_res)
  mse_atmle <- dt_res[dt_res$estimator == "A-TMLE", "mse"]
  mse_atmle_star <- dt_res[dt_res$estimator == "A-TMLE*", "mse"]
  mse_escvtmle <- dt_res[dt_res$estimator == "ESCVTMLE", "mse"]
  dt_relative <- data.frame(n = rep(seq(n_min, n_max, n_step), 2),
                            names = rep(c("mse(ESCVTMLE)/mse(A-TMLE)", "mse(ESCVTMLE)/mse(A-TMLE*)"),
                                        each = length(seq(n_min, n_max, n_step))),
                            ratio = c(mse_escvtmle / mse_atmle, mse_escvtmle / mse_atmle_star))
  relative_mse_plot <- ggplot(dt_relative, aes(x = n, y = ratio, color = names)) +
    geom_point() +
    geom_line() +
    labs(title = name,
         x = "n",
         y = "relative mse") +
    theme_minimal() +
    theme(text = element_text(size = 16))
}
