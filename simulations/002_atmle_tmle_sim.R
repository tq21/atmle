library(purrr)
source("utils.R")

# n increase -------------------------------------------------------------------
B <- 500
n_min <- 200
n_max <- 3000
n_step <- 200
truth <- 1.5
n_seq <- seq(n_min, n_max, n_step)

no_bias <- run_sim_n_increase(B = B,
                              n_min = n_min,
                              n_max = n_max,
                              n_step = n_step,
                              bA = truth,
                              bias = 0,
                              nuisance_method = "glm",
                              working_model = "lasso",
                              verbose = TRUE,
                              atmle_both = FALSE)

small_bias <- run_sim_n_increase(B = B,
                                 n_min = n_min,
                                 n_max = n_max,
                                 n_step = n_step,
                                 bA = truth,
                                 bias = 0.5,
                                 nuisance_method = "glm",
                                 working_model = "lasso",
                                 verbose = TRUE,
                                 atmle_both = FALSE)

large_bias <- run_sim_n_increase(B = B,
                                 n_min = n_min,
                                 n_max = n_max,
                                 n_step = n_step,
                                 bA = truth,
                                 bias = 1.8,
                                 nuisance_method = "glm",
                                 working_model = "lasso",
                                 verbose = TRUE,
                                 atmle_both = FALSE)

param_bias <- run_sim_n_increase(B = B,
                                 n_min = n_min,
                                 n_max = n_max,
                                 n_step = n_step,
                                 bA = truth,
                                 bias = "parametric",
                                 nuisance_method = "glm",
                                 working_model = "lasso",
                                 verbose = TRUE,
                                 atmle_both = FALSE)

save(list = c("no_bias", "small_bias", "large_bias", "param_bias",
              "B", "truth", "n_seq"),
     file = "out/sim_atmle_tmle_" %+% format(Sys.time(), "%Y%m%d_%H%M%S") %+% ".RData")
