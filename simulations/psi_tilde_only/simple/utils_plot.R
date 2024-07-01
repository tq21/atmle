library(ggplot2)
library(ggpubr)

get_bias_plot <- function(df) {
  p_bias <- ggplot(df, aes(x = n, y = bias, color = estimator)) +
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
  p_var <- ggplot(df, aes(x = n, y = var, color = estimator)) +
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
  p_mse <- ggplot(df, aes(x = n, y = mse, color = estimator)) +
    geom_point(size = 1.5) +
    geom_line(linewidth = 1) +
    labs(title = "",
         x = "n",
         y = "MSE") +
    theme_minimal() +
    theme(text = element_text(size = 16),
          legend.position = "none")

  return(p_mse)
}

get_cover_plot <- function(df) {
  p_cover <- ggplot(df, aes(x = n, y = cover, color = estimator)) +
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
  relative_mse_plot <- ggplot(df, aes(x = n, y = relative_mse, color = estimator)) +
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
