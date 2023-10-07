#' Simulations comparing the following estimators:
#' 1. A-TMLE both psi_pound and psi_tilde
#' 2. A-TMLE psi_pound, regular TMLE psi_tilde
#' 3. ESCVTMLE
#' 4. Nonparametric TMLE
#' Data generating distribution scenario: simple parametric bias

source("utils.R")

# simulation parameters --------------------------------------------------------
B <- 500 # number of runs for each sample size n
n_min <- 1000 # smallest sample size
n_max <- 3000 # largest sample size
n_step <- 500 # sample size increment
bA <- 1.5 # true ATE
bias <- "param_simple" # true bias
nuisance_method <- "glm"
working_model <- "glmnet"
g_rct <- 0.67
f_name <- "param_simple_glm_glmnet"
date_name <- "1007"
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
