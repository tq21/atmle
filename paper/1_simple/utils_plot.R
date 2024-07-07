`%+%` <- function(a, b) paste0(a, b)

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
                    Estimator = estimator_name,
                    bias = bias,
                    se = se,
                    mse = mse,
                    cover = cover))
}

get_mse_plot <- function(title, label, names, ...) {
  res_list <- list(...)
  dt_res <- map_dfr(1:length(res_list), function(i) {
    return(get_res(res_list[[i]], ate, names[i]))
  })
  dt_res <- dt_res[order(dt_res$n), ]

  p_mse <- ggplot(dt_res, aes(x = n, y = mse, color = Estimator)) +
    geom_point(size = 1.5) +
    geom_line(linewidth = 1) +
    scale_x_continuous(breaks = total_sample_sizes) +
    labs(title = "",
         x = "n",
         y = (label %+% "\nMSE")) +
    theme_minimal() +
    theme(text = element_text(size = 16),
          legend.position = "none")

  return(p_mse)
}

get_cover_plot <- function(title, names, ...) {
  res_list <- list(...)
  dt_res <- map_dfr(1:length(res_list), function(i) {
    return(get_res(res_list[[i]], ate, names[i]))
  })
  dt_res <- dt_res[order(dt_res$n), ]

  p_cover <- ggplot(dt_res, aes(x = n, y = cover, color = Estimator)) +
    geom_point(size = 1.5) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = 0.95, color = "red", linetype = "dashed", linewidth = 1) +
    scale_x_continuous(breaks = total_sample_sizes) +
    scale_y_continuous(breaks = seq(0, 1.0, 0.1), limits = c(0, 1.0)) +
    labs(title = "",
         x = "n",
         y = "Coverage") +
    theme_minimal() +
    theme(text = element_text(size = 16),
          legend.position = "none")

  return(p_cover)
}

get_plot <- function(title, names, ...) {
  res_list <- list(...)
  dt_res <- map_dfr(1:length(res_list), function(i) {
    return(get_res(res_list[[i]], ate, names[i]))
  })
  dt_res <- dt_res[order(dt_res$n), ]

  p_bias <- ggplot(dt_res, aes(x = n, y = abs(bias), color = Estimator)) +
    geom_point() +
    geom_line() +
    labs(title = "",
         x = "n",
         y = "bias") +
    theme_minimal() +
    theme(text = element_text(size = 16))

  p_se <- ggplot(dt_res, aes(x = n, y = se, color = Estimator)) +
    geom_point() +
    geom_line() +
    labs(title = "",
         x = "n",
         y = "standard error") +
    theme_minimal() +
    theme(text = element_text(size = 16))

  p_mse <- ggplot(dt_res, aes(x = n, y = mse, color = Estimator)) +
    geom_point() +
    geom_line() +
    labs(title = "",
         x = "n",
         y = "mse") +
    theme_minimal() +
    theme(text = element_text(size = 16))

  p_cover <- ggplot(dt_res, aes(x = n, y = cover, color = Estimator)) +
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
    geom_point(size = 1.5, color = "#7CAE00") +
    geom_line(linewidth = 1, color = "#7CAE00") +
    #scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
    scale_x_continuous(breaks = total_sample_sizes) +
    labs(title = "",
         x = "n",
         y = "Avg. prop. folds selected") +
    theme_minimal() +
    theme(text = element_text(size = 16),
          legend.position = "none")

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
    dt_relative$ratio[dt_relative$names == comparator] <- dt_res$mse[dt_res$Estimator == reference] / dt_res$mse[dt_res$Estimator == comparator]
  }

  relative_mse_plot <- ggplot(dt_relative, aes(x = n, y = ratio, color = names)) +
    geom_point(size = 1.5) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = 1, color = "red", linetype = "dashed", linewidth = 1) +
    scale_x_continuous(breaks = total_sample_sizes) +
    scale_y_continuous(breaks = seq(0.8, 2, 0.2), limits = c(0.8, 2)) +
    labs(title = title,
         x = "n",
         y = "Relative MSE") +
    theme_minimal() +
    theme(text = element_text(size = 16),
          legend.position = "none")
}

# gg_color_hue <- function(n) {
#   hues = seq(15, 375, length = n + 1)
#   hcl(h = hues, l = 65, c = 100)[1:n]
# }
# n = 4
# cols = gg_color_hue(n)

get_avg_ci_length <- function(obj) {
  return(mean(unlist(obj$all_psi_ci_upper) - unlist(obj$all_psi_ci_lower)))
}
