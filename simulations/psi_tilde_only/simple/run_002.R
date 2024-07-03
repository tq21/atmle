source("utils.R")
source("sim_data_002.R")

set.seed(283947)
B <- 500
n_seq <- seq(500, 2000, 500)
t0 <- 3

data_list <- make_data(B = B, n_seq = n_seq)
truth <- get_truth(t0)

res_list <- run_sim(t0 = t0,
                    truth = truth,
                    data_list = data_list,
                    g_method = "glm",
                    lambda_method = "glm",
                    working_model = "glmnet",
                    cross_fit_nuisance = TRUE)
saveRDS(res_list, "out/res_list_002.rds")
