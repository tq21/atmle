library(devtools)
library(dplyr)
library(EScvtmle)
library(doMC)
library(furrr)
library(ggplot2)
load_all()
source("sim_data.R")

registerDoMC(cores = availableCores()-1)
set.seed(123)
data_rct <- sim_data(ate = 1.5,
                     n = 400,
                     rct = TRUE,
                     g_rct = 0.5,
                     bias = "b",
                     controls_only = FALSE)
data_rwd <- sim_data(ate = 1.5,
                     n = 1200,
                     rct = FALSE,
                     g_rct = 0.5,
                     bias = "b",
                     controls_only = FALSE)
data <- rbind(data_rct, data_rwd)

B <- 50
psi <- lower <- upper <- numeric(B)

res_df <- map_dfr(seq(B), function(.b) {
  print(.b)
  res <- atmle(data = data,
               S = "S",
               W = c("W1", "W2", "W3"),
               A = "A",
               Y = "Y",
               controls_only = FALSE,
               atmle_pooled = TRUE,
               family = "gaussian",
               theta_method = "glm",
               g_method = "glm",
               theta_tilde_method = "glm",
               bias_working_model = "HAL",
               pooled_working_model = "HAL",
               max_degree = 1,
               verbose = FALSE,
               browse = FALSE)

  return(data.frame(psi = res$est,
                    lower = res$lower,
                    upper = res$upper))
})

# plot point estimates with error bars
res_df %>%
  ggplot(aes(x = 1:B, y = psi)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  geom_hline(yintercept = 1.5, linetype = "dashed") +
  labs(title = "",
       x = "",
       y = "ATE estimate") +
  theme_minimal()
