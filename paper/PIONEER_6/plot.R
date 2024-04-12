# res_df <- data.frame(
#   Estimator = c("ES-CVTMLE", "A-TMLE", "A-TMLE (psi tilde)", "A-TMLE (psi pound)"),
#   est = (c(-0.0153, -0.0107374, -0.0121455, -0.0014081)*100),
#   lower = (c(-0.0275, -0.0199295, -0.0236712, -0.0097927)*100),
#   upper = (c(-0.003, -0.0015454, -0.0006199, 0.0069765)*100)
# )

# res_df <- data.frame(
#   Estimator = c("RCT", "ES-CVTMLE", "A-TMLE", "A-TMLE (psi tilde)"),
#   est = (c(-0.013, -0.0123, -0.0107374, -0.0121455)*100),
#   lower = (c(-0.026, -0.0253, -0.0199295, -0.0236712)*100),
#   upper = (c(0, -0.000185, -0.0015454, -0.0006199)*100)
# )

res_df <- data.frame(
  Estimator = c("RCT", "ES-CVTMLE", "A-TMLE"),
  est = (c(-0.013, -0.0123, -0.0107374)*100),
  lower = (c(-0.026, -0.0253, -0.0199295)*100),
  upper = (c(0, -0.000185, -0.0015454)*100)
)

ggplot(res_df, aes(x = est, y = Estimator)) +
  geom_point(position = position_dodge(width = 0.5), size = 3.5) +
  geom_errorbar(aes(xmin = lower, xmax = upper),
                position = position_dodge(width = 0.5), width = 0.3, size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.5) +
  labs(title = "Estimated Difference in 1-Year Risk of MACE",
       x = "Risk difference (%)", y = "", color = "") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 14))
