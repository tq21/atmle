library(ggplot2)
source("sim_data.R")

res_seq_atmle_df <- read.csv("out/atmle_seq_ate_0216_101151.csv")
truth <- get_truth()

ggplot(res_seq_atmle_df, aes(x = index, y = psi)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  facet_wrap(~selector) +
  geom_hline(yintercept = truth, linetype = "dashed", color = "red") +
  labs(x = "Index", y = "Psi Estimate",
       title = "Sequence of Psi Estimates with 95% Confidence Intervals") +
  theme_minimal()

abs(res_seq_atmle_df$PnEIC) <= res_seq_atmle_df$sn
