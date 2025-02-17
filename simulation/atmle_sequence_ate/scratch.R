library(dplyr)
source("sim_data.R")
truth <- get_truth()
res_df <- read.csv("out/atmle_seq_ate_0216_101151.csv")
res_df %>% summarize(abs_bias = abs(mean(psi)-truth),
                     se = sd(psi),
                     mse = mean((psi-truth)^2),
                     cover = mean(lower <= truth & truth <= upper),
                     oracle_cover = mean(psi-1.96*sd(psi) <= truth & truth <= psi+1.96*sd(psi)),
                     ci_width = mean(upper - lower),
                     .by = c("n", "selector"))
