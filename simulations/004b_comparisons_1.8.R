#' Simulations comparing the following estimators:
<<<<<<< Updated upstream:simulations/004b_comparisons_1.8.R
#' 1. A-TMLE both psi_pound and psi_tilde
#' 2. A-TMLE psi_pound, regular TMLE psi_tilde
#' 3. ESCVTMLE
#' 4. Nonparametric TMLE
#' Data generating distribution scenario: 1.8 constant bias
=======
#' 1. A-TMLE, both psi_pound and psi_tilde
#' 2. ES-CVTMLE
#' 3. Nonparametric TMLE
#' 4. RCT only
>>>>>>> Stashed changes:paper/simulation_2/a_HAL_1.R

source("utils.R")

# simulation parameters --------------------------------------------------------
<<<<<<< Updated upstream:simulations/004b_comparisons_1.8.R
B <- 500 # number of runs for each sample size n
n_min <- 1000 # smallest sample size
n_max <- 3000 # largest sample size
n_step <- 500 # sample size increment
bA <- 1.5 # true ATE
bias <- 1.8 # true bias
nuisance_method <- "glm"
working_model <- "glmnet"
g_rct <- 0.67
f_name <- "1.8_bias_glm_glmnet"
date_name <- "1006_pos"
controls_only <- FALSE

# 1. A-TMLE both psi_pound and psi_tilde
atmle_both_res <- run_sim_n_increase(B = B,
                                     n_min = n_min,
                                     n_max = n_max,
                                     n_step = n_step,
                                     bA = bA,
                                     bias = bias,
                                     g_rct = g_rct,
                                     controls_only = controls_only,
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
                                     g_rct = g_rct,
                                     controls_only = controls_only,
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
                                   g_rct = g_rct,
                                   controls_only = controls_only,
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
                                   g_rct = g_rct,
                                   controls_only = controls_only,
                                   nuisance_method = nuisance_method,
                                   working_model = working_model,
                                   verbose = TRUE,
                                   method = "rct_only")
=======
set.seed(123)
B <- 200 # number of runs for each sample size n
n_rct_seq <- 400 #c(400, 500, 600, 700, 800)
n_rwd_seq <- n_rct_seq * 3
ate <- 1.5 # true ATE
bias <- "HAL_1" # true bias
controls_only <- FALSE
nuisance_method <- "glm"
working_model <- "glmnet"
g_rct <- 0.67
f_name <- bias %+% "_bias"
date_name <- "0218"

# make data
data_list <- make_data(B = B,
                       n_rct_seq = n_rct_seq,
                       n_rwd_seq = n_rwd_seq,
                       ate = ate,
                       bias = bias,
                       controls_only = controls_only)

# 1. A-TMLE
atmle_both_res <- run_sim(data_list = data_list,
                          ate = ate,
                          controls_only = controls_only,
                          nuisance_method = nuisance_method,
                          working_model = working_model,
                          g_rct = g_rct,
                          method = "atmle")
saveRDS(atmle_both_res,
        file = "out/atmle_both_" %+% f_name %+% "_" %+% date_name %+% ".RDS")

# 2. ESCVTMLE
escvtmle_res <- run_sim(data_list = data_list,
                        ate = ate,
                        controls_only = controls_only,
                        nuisance_method = nuisance_method,
                        working_model = working_model,
                        g_rct = g_rct,
                        method = "escvtmle")
saveRDS(escvtmle_res,
        file = "out/escvtmle_" %+% f_name %+% "_" %+% date_name %+% ".RDS")

# 3. Nonparametric TMLE
tmle_res <- run_sim(data_list = data_list,
                    ate = ate,
                    controls_only = controls_only,
                    nuisance_method = nuisance_method,
                    working_model = working_model,
                    g_rct = g_rct,
                    method = "tmle")
saveRDS(tmle_res,
        file = "out/tmle_" %+% f_name %+% "_" %+% date_name %+% ".RDS")

# 4. RCT only
rct_only_res <- run_sim(data_list = data_list,
                        ate = ate,
                        controls_only = controls_only,
                        nuisance_method = nuisance_method,
                        working_model = working_model,
                        g_rct = g_rct,
                        method = "rct_only")
>>>>>>> Stashed changes:paper/simulation_2/a_HAL_1.R
saveRDS(rct_only_res,
        file = "out/rct_only_" %+% f_name %+% "_" %+% date_name %+% ".RDS")

# 5. Nonparametric TMLE
tmle_res <- run_sim_n_increase(B = B,
                               n_min = n_min,
                               n_max = n_max,
                               n_step = n_step,
                               bA = bA,
                               bias = bias,
                               g_rct = g_rct,
                               controls_only = controls_only,
                               nuisance_method = nuisance_method,
                               working_model = working_model,
                               verbose = TRUE,
                               method = "tmle")
saveRDS(tmle_res,
        file = "out/tmle_" %+% f_name %+% "_" %+% date_name %+% ".RDS")
