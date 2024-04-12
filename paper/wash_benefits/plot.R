load("out/wash_results.RData")

res_df_unbiased$bias <- "unbiased external"
res_df_biased$bias <- "biased external"
res_df <- rbind(res_df_unbiased, res_df_biased)
res_df$estimator <- factor(res_df$estimator, levels = c("ES-CVTMLE", "A-TMLE"))

# plot (unbiased external)
p_unbiased <- ggplot(res_df_unbiased, aes(x = estimator, y = est)) +
  geom_point(position = position_dodge(width = 0.5), size = 3.5, color = "dodgerblue") +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                position = position_dodge(width = 0.5), width = 0.3, size = 1.5, color = "dodgerblue") +
  scale_y_continuous(limits = c(-0.3, 0.3), breaks = seq(-0.3, 0.3, 0.1),
                     labels = c(-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.5) +
  labs(title = "Unbiased external", x = "", y = "Estimate", color = "") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 14))

# plot (biased external)
p_biased <- ggplot(res_df_biased, aes(x = estimator, y = est)) +
  geom_point(position = position_dodge(width = 0.5), size = 3.5, color = "darkorange") +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                position = position_dodge(width = 0.5), width = 0.3, size = 1.5, color = "darkorange") +
  scale_y_continuous(limits = c(-0.3, 0.3), breaks = seq(-0.3, 0.3, 0.1),
                     labels = c(-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.5) +
  labs(title = "Biased external", x = "", y = "Estimate", color = "") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 14))

wash_plt <- ggarrange(p_unbiased, p_biased, nrow = 1, ncol = 2)
ggsave(filename = "wash_results.pdf", plot = wash_plt, device = "pdf",
       path = "figs", width = 10, height = 8, dpi = 300)

# get widths of the confidence intervals
res_df_unbiased$width <- res_df_unbiased$upper - res_df_unbiased$lower
res_df_biased$width <- res_df_biased$upper - res_df_biased$lower
atmle_unbiased_ci_width <- res_df_unbiased[res_df_unbiased$estimator == "A-TMLE", "width"]
escvtmle_unbiased_ci_width <- res_df_unbiased[res_df_unbiased$estimator == "ES-CVTMLE", "width"]
atmle_biased_ci_width <- res_df_biased[res_df_biased$estimator == "A-TMLE", "width"]
escvtmle_biased_ci_width <- res_df_biased[res_df_biased$estimator == "ES-CVTMLE", "width"]

print("Unbiased, A-TMLE relative to ES-CVTMLE: " %+%
        ((escvtmle_unbiased_ci_width-atmle_unbiased_ci_width)/escvtmle_unbiased_ci_width))
print("Biased, A-TMLE relative to ES-CVTMLE: " %+%
        ((escvtmle_biased_ci_width-atmle_biased_ci_width)/escvtmle_biased_ci_width))
