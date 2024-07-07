get_mse_plot <- function(df) {
  p_mse <- ggplot(df, aes(x = n, y = mse, color = Estimator)) +
    geom_point(size = 1.5) +
    geom_line(linewidth = 1) +
    scale_x_continuous(breaks = n_seq) +
    labs(title = "",
         x = "n",
         y = "MSE") +
    theme_minimal() +
    theme(text = element_text(size = 16),
          legend.position = "none")

  return(p_mse)
}
