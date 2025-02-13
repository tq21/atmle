library(dplyr)
res_df <- read.csv("out/run_bias_c.csv")

res_df %>%
  summarize(mse = mean((psi-4.2)^2),
            ci_length = mean(upper - lower),
            coverage = mean((lower <= 4.2) & (4.2 <= upper)),
            .by = "est_name")
