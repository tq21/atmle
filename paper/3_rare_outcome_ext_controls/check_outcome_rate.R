source("sim_data.R")

n <- 100000
data_S1 <- sim_data(beta = 1.5,
                    n = n,
                    rct = TRUE,
                    g_rct = 0.67,
                    bias = "a")
data_S0 <- sim_data(beta = 1.5,
                    n = n,
                    rct = FALSE,
                    g_rct = 0.67,
                    bias = "a")

mean(data_S1$Y)
mean(data_S0$Y)
