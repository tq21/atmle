ggplot(res_seq_atmle_df, aes(x = index, y = psi)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  facet_wrap(~eic_method) +
  geom_hline(yintercept = truth, linetype = "dashed", color = "red") +
  labs(x = "Index", y = "Psi Estimate",
       title = "Sequence of Psi Estimates with 95% Confidence Intervals") +
  theme_minimal()

abs(res_seq_atmle_df$PnEIC) <= res_seq_atmle_df$sn
