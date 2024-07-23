source("sim_data.R")
library(dplyr)
load_all()

set.seed(123)
B <- 500
est_random <- numeric(B)
est_oracle <- numeric(B)
cover_random <- numeric(B)
cover_oracle <- numeric(B)

for (b in 1:B) {
  print(b)
  data <- sim_data(N_rwd = 100000,
                   n_rct = 300,
                   n_rwd = 2000,
                   g_rct = 0.5,
                   ate = 1.5)
  res_random <- atmle(data = data$random,
                      S = "S",
                      W = c("W1", "W2", "W3"),
                      A = "A",
                      Y = "Y",
                      controls_only = FALSE,
                      family = "gaussian",
                      theta_method = "glm",
                      Pi_method = "glm",
                      g_method = "glm",
                      theta_tilde_method = "glm",
                      bias_working_model = "glmnet",
                      pooled_working_model = "glmnet",
                      verbose = FALSE)
  res_oracle <- atmle(data = data$oracle,
                      S = "S",
                      W = c("W1", "W2", "W3"),
                      A = "A",
                      Y = "Y",
                      controls_only = FALSE,
                      family = "gaussian",
                      theta_method = "glm",
                      Pi_method = "glm",
                      g_method = "glm",
                      theta_tilde_method = "glm",
                      bias_working_model = "glmnet",
                      pooled_working_model = "glmnet",
                      verbose = FALSE)

  if (res_random$lower <= 1.5 & res_random$upper >= 1.5) {
    cover_random[b] <- 1
  } else {
    cover_random[b] <- 0
  }

  if (res_oracle$lower <= 1.5 & res_oracle$upper >= 1.5) {
    cover_oracle[b] <- 1
  } else {
    cover_oracle[b] <- 0
  }

  est_random[b] <- res_random$est
  est_oracle[b] <- res_oracle$est
}

mean((est_random - 1.5)^2)
mean((est_oracle - 1.5)^2)
