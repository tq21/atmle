res <- data.frame(
  name = c("atmle", "atmle*", "escvtmle"),
  est = c(-0.00852, -0.00797, -0.0155),
  lower = c(-0.0229, -0.0233, -0.0288),
  upper = c(0.00586, 0.00737, -0.00388)
)

ggplot(res, aes(x = name, y = est)) +
  geom_point(position = position_dodge(width = 0.5), size = 3.5, color = "darkorange") +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                position = position_dodge(width = 0.5), width = 0.3, size = 1.5, color = "darkorange") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.5) +
  labs(title = "", x = "", y = "Estimate", color = "") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 14))
