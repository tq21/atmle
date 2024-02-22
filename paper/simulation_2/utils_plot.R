get_bias <- function(res, ate) {
  return(unlist(map(res$all_psi_est, function(.x) mean(.x-ate))))
}

get_se <- function(res) {
  return(unlist(map(res$all_psi_est, function(.x) sd(.x))))
}

get_var <- function(res) {
  return(unlist(map(res$all_psi_est, function(.x) var(.x))))
}

get_cover <- function(res) {
  return(unlist(map(res$all_psi_coverage, function(.x) mean(.x))))
}

get_res <- function(res, ate, estimator_name) {

  bias <- get_bias(res, ate)
  se <- get_se(res)
  var <- get_var(res)
  mse <- bias^2+var
  cover <- get_cover(res)

  return(data.frame(n = total_sample_sizes,
                    estimator = estimator_name,
                    bias = bias,
                    se = se,
                    mse = mse,
                    cover = cover))
}

get_mse_plot <- function(title, names, ...) {
  res_list <- list(...)
  dt_res <- map_dfr(1:length(res_list), function(i) {
    return(get_res(res_list[[i]], ate, names[i]))
  })
  dt_res <- dt_res[order(dt_res$n), ]

  p_mse <- ggplot(dt_res, aes(x = n, y = mse, color = estimator)) +
    geom_point() +
    geom_line() +
    labs(title = "",
         x = "n",
         y = "mse") +
    theme_minimal() +
    theme(text = element_text(size = 16))

  plt <- ggarrange(p_mse, nrow = 1, ncol = 1)
  plt <- annotate_figure(plt, top = text_grob(title, face = "bold", size = 16))

  return(plt)
}

get_cover_plot <- function(title, names, ...) {
  res_list <- list(...)
  dt_res <- map_dfr(1:length(res_list), function(i) {
    return(get_res(res_list[[i]], ate, names[i]))
  })
  dt_res <- dt_res[order(dt_res$n), ]

  p_cover <- ggplot(dt_res, aes(x = n, y = cover, color = estimator)) +
    geom_point() +
    geom_line() +
    geom_hline(yintercept = 0.95, color = "red", linetype = "dashed", linewidth = 1) +
    scale_y_continuous(breaks = seq(0.8, 1.0, 0.05), limits = c(0.8, 1.0)) +
    labs(title = "",
         x = "n",
         y = "coverage") +
    theme_minimal() +
    theme(text = element_text(size = 16))

  plt <- ggarrange(p_cover, nrow = 1, ncol = 1)
  plt <- annotate_figure(plt, top = text_grob(title, face = "bold", size = 16))

  return(plt)
}

get_plot <- function(title, names, ...) {
  res_list <- list(...)
  dt_res <- map_dfr(1:length(res_list), function(i) {
    return(get_res(res_list[[i]], ate, names[i]))
  })
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

  p_cover <- ggplot(dt_res, aes(x = n, y = cover, color = estimator)) +
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

# average proportion of folds selected plot for each n
get_plot_prop_selected <- function(escvtmle_res, name) {
  df <- data.frame(n = total_sample_sizes,
                   prop = map_vec(escvtmle_res$all_escvtmle_prop_selected, mean))
  p <- ggplot(df, aes(x = n, y = prop)) +
    geom_point() +
    geom_line() +
    scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
    labs(title = "",
         x = "n",
         y = "Average proportion of folds that pooled by ES-CVTMLE") +
    theme_minimal() +
    theme(text = element_text(size = 16))

  return(p)
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

get_relative_mse_plot <- function(title, estimator_names, comparisons, ...) {
  res_list <- list(...)
  dt_res <- map_dfr(1:length(res_list), function(i) {
    return(get_res(res_list[[i]], ate, estimator_names[i]))
  })

  # get second estimator name of all comparisons
  comparator_names <- unlist(map(comparisons, function(.x) .x[1]))

  # compute relative MSE, for each comparison, first is comparator, second is reference
  dt_relative <- data.frame(n = dt_res$n,
                            names = rep(comparator_names, each = length(dt_res$n)),
                            ratio = rep(NA, length(dt_res$n) * length(comparisons)))
  for (i in 1:length(comparisons)) {
    comparator <- comparisons[[i]][1]
    reference <- comparisons[[i]][2]
    dt_relative$ratio[dt_relative$names == comparator] <- dt_res$mse[dt_res$estimator == reference] / dt_res$mse[dt_res$estimator == comparator]
  }

  relative_mse_plot <- ggplot(dt_relative, aes(x = n, y = ratio, color = names)) +
    geom_point(size = 1.5) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = 1, color = "red", linetype = "dashed", linewidth = 1) +
    scale_y_continuous(breaks = seq(1, 2, 0.5), limits = c(1, 2)) +
    labs(title = title,
         x = "n",
         y = "Relative MSE") +
    theme_minimal() +
    theme(text = element_text(size = 16),
          plot.title = element_text(face = "bold", size = 16, hjust = 0.5))
}
