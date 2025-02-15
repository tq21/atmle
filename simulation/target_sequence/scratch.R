library(dplyr)
res_df <- read.csv("out/run_bias_c.csv")

res_df %>%
  summarize(mse = mean((psi-4.2)^2),
            ci_length = mean(upper - lower),
            coverage = mean((lower <= 4.2) & (4.2 <= upper)),
            .by = "est_name")

library(doMC)
registerDoMC(cores = 7)
x = matrix(rnorm(1e+05 * 100), 1e+05, 100)
y = rnorm(1e+05)
system.time(cv.glmnet(x, y))
system.time(cv.glmnet(x, y, parallel = TRUE))
