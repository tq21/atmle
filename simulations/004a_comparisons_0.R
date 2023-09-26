#' Simulations comparing the following estimators:
#' 1. A-TMLE both psi_pound and psi_tilde
#' 2. A-TMLE psi_pound, regular TMLE psi_tilde
#' 3. ESCVTMLE
#' 4. Nonparametric TMLE
#' Data generating distribution scenario: no bias

source("utils.R")

# simulation parameters --------------------------------------------------------
B <- 500 # number of runs for each sample size n
n_min <- 500 # smallest sample size
n_max <- 3000 # largest sample size
n_step <- 500 # sample size increment
bA <- 1.5 # true ATE
bias <- 0 # true bias
nuisance_method <- "glm"
working_model <- "lasso"
pRCT <- 0.67
f_name <- "0_bias_glm_lasso"
date_name <- "926"

# 1. A-TMLE both psi_pound and psi_tilde
atmle_both_res <- run_sim_n_increase(B = B,
                                     n_min = n_min,
                                     n_max = n_max,
                                     n_step = n_step,
                                     bA = bA,
                                     bias = bias,
                                     pRCT = pRCT,
                                     nuisance_method = nuisance_method,
                                     working_model = working_model,
                                     verbose = TRUE,
                                     method = "atmle")
saveRDS(atmle_both_res,
        file = "out/atmle_both_" %+% f_name %+% "_" %+% date_name %+% ".RDS")

# 2. A-TMLE psi_pound, regular TMLE psi_tilde
atmle_tmle_res <- run_sim_n_increase(B = B,
                                     n_min = n_min,
                                     n_max = n_max,
                                     n_step = n_step,
                                     bA = bA,
                                     bias = bias,
                                     pRCT = pRCT,
                                     nuisance_method = nuisance_method,
                                     working_model = working_model,
                                     verbose = TRUE,
                                     method = "atmle_tmle")
saveRDS(atmle_tmle_res,
        file = "out/atmle_tmle_" %+% f_name %+% "_" %+% date_name %+% ".RDS")

# 3. ESCVTMLE
escvtmle_res <- run_sim_n_increase(B = B,
                                   n_min = n_min,
                                   n_max = n_max,
                                   n_step = n_step,
                                   bA = bA,
                                   bias = bias,
                                   pRCT = pRCT,
                                   nuisance_method = nuisance_method,
                                   working_model = working_model,
                                   verbose = TRUE,
                                   method = "escvtmle")
saveRDS(escvtmle_res,
        file = "out/escvtmle_" %+% f_name %+% "_" %+% date_name %+% ".RDS")

# 4. RCT only
rct_only_res <- run_sim_n_increase(B = B,
                                   n_min = n_min,
                                   n_max = n_max,
                                   n_step = n_step,
                                   bA = bA,
                                   bias = bias,
                                   pRCT = pRCT,
                                   nuisance_method = nuisance_method,
                                   working_model = working_model,
                                   verbose = TRUE,
                                   method = "rct_only")
saveRDS(rct_only_res,
        file = "out/rct_only_" %+% f_name %+% "_" %+% date_name %+% ".RDS")

# # 4. Nonparametric TMLE
# tmle_res <- run_sim_n_increase(B = B,
#                                n_min = n_min,
#                                n_max = n_max,
#                                n_step = n_step,
#                                bA = bA,
#                                bias = bias,
#                                nuisance_method = nuisance_method,
#                                working_model = working_model,
#                                verbose = TRUE,
#                                method = "tmle")
# saveRDS(tmle_res,
#         file = "out/tmle_1.8_bias_glm_lasso_" %+% format(Sys.time(), "%Y%m%d_%H%M%S") %+% ".RDS")
