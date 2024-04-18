#' Simulations comparing the following estimators:
#' 1. A-TMLE, both psi_pound and psi_tilde
#' 2. ES-CVTMLE
#' 3. Nonparametric TMLE
#' 4. RCT only

source("utils.R")

# simulation parameters --------------------------------------------------------
set.seed(123)
B <- 500 # number of runs for each sample size n
n_rct_seq <- seq(1000, 1500, 100)
n_rwd_seq <- n_rct_seq * 5
beta <- 0.9 # beta
bias <- "a" # true bias
nuisance_method <- "glm"
working_model <- "glmnet"
g_rct <- 0.5
f_name <- bias %+% "_bias"
date_name <- "0418"

# make data
data_list <- make_data(B = B,
                       n_rct_seq = n_rct_seq,
                       n_rwd_seq = n_rwd_seq,
                       beta = beta,
                       bias = bias)

# check outcome rate
outcome_rate <- sapply(unlist(data_list, recursive = FALSE), function(x) mean(x$Y))

# 1. A-TMLE
atmle_both_res <- run_sim(data_list = data_list,
                          beta = beta,
                          nuisance_method = nuisance_method,
                          working_model = working_model,
                          g_rct = g_rct,
                          method = "atmle")
saveRDS(atmle_both_res,
        file = "out/atmle_both_" %+% f_name %+% "_" %+% date_name %+% ".RDS")

# 2. ESCVTMLE
# escvtmle_res <- run_sim(data_list = data_list,
#                         beta = beta,
#                         nuisance_method = nuisance_method,
#                         working_model = working_model,
#                         g_rct = g_rct,
#                         method = "escvtmle")
# saveRDS(escvtmle_res,
#         file = "out/escvtmle_" %+% f_name %+% "_" %+% date_name %+% ".RDS")
#
# # 3. Nonparametric TMLE
# tmle_res <- run_sim(data_list = data_list,
#                     beta = beta,
#                     nuisance_method = nuisance_method,
#                     working_model = working_model,
#                     g_rct = g_rct,
#                     method = "tmle")
# saveRDS(tmle_res,
#         file = "out/tmle_" %+% f_name %+% "_" %+% date_name %+% ".RDS")
#
# # 4. RCT only
# rct_only_res <- run_sim(data_list = data_list,
#                         beta = beta,
#                         nuisance_method = nuisance_method,
#                         working_model = working_model,
#                         g_rct = g_rct,
#                         method = "rct_only")
# saveRDS(rct_only_res,
#         file = "out/rct_only_" %+% f_name %+% "_" %+% date_name %+% ".RDS")
#
# # save all environment variables
#  #save.image(file = "out/env_" %+% f_name %+% "_" %+% date_name %+% ".RData")
