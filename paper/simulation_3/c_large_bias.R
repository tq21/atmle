#' Simulations comparing the following estimators:
#' 1. A-TMLE, both psi_pound and psi_tilde
#' 2. A-TMLE psi_pound, regular TMLE psi_tilde
#' 3. ESCVTMLE
#' 4. Nonparametric TMLE
#' 5. RCT only
#' Data generating distribution scenario: large bias

source("utils.R")
source("sim_data.R")

# simulation parameters --------------------------------------------------------
B <- 1000 # number of runs each time
rct_prop_min <- 0.05 # smallest proportion of RCT data
rct_prop_max <- 0.5 # largest proportion of RCT data
rct_prop_step <- 0.05 # increment of proportion of RCT data
ate <- 1.5 # true ATE
bias <- "large_bias" # true bias
controls_only <- FALSE
nuisance_method <- "glm"
working_model <- "glmnet"
g_rct <- 0.67
num_covs <- 4
f_name <- "large_bias"
date_name <- "1217"

# 1. A-TMLE
atmle_both_res <- run_sim(B = B,
                          rct_prop_min = rct_prop_min,
                          rct_prop_max = rct_prop_max,
                          rct_prop_step = rct_prop_step,
                          ate = ate,
                          bias = bias,
                          controls_only = controls_only,
                          nuisance_method = nuisance_method,
                          working_model = working_model,
                          g_rct = g_rct,
                          num_covs = num_covs,
                          method = "atmle")
saveRDS(atmle_both_res,
        file = "out/atmle_both_" %+% f_name %+% "_" %+% date_name %+% ".RDS")

# 2. A-TMLE psi_pound, regular TMLE psi_tilde
atmle_tmle_res <- run_sim(B = B,
                          rct_prop_min = rct_prop_min,
                          rct_prop_max = rct_prop_max,
                          rct_prop_step = rct_prop_step,
                          ate = ate,
                          bias = bias,
                          controls_only = controls_only,
                          nuisance_method = nuisance_method,
                          working_model = working_model,
                          g_rct = g_rct,
                          num_covs = num_covs,
                          method = "atmle_tmle")
saveRDS(atmle_tmle_res,
        file = "out/atmle_tmle_" %+% f_name %+% "_" %+% date_name %+% ".RDS")

# 3. ESCVTMLE
escvtmle_res <- run_sim(B = B,
                        rct_prop_min = rct_prop_min,
                        rct_prop_max = rct_prop_max,
                        rct_prop_step = rct_prop_step,
                        ate = ate,
                        bias = bias,
                        controls_only = controls_only,
                        nuisance_method = nuisance_method,
                        working_model = working_model,
                        g_rct = g_rct,
                        num_covs = num_covs,
                        method = "escvtmle")
saveRDS(escvtmle_res,
        file = "out/escvtmle_" %+% f_name %+% "_" %+% date_name %+% ".RDS")

# 4. Nonparametric TMLE
tmle_res <- run_sim(B = B,
                    rct_prop_min = rct_prop_min,
                    rct_prop_max = rct_prop_max,
                    rct_prop_step = rct_prop_step,
                    ate = ate,
                    bias = bias,
                    controls_only = controls_only,
                    nuisance_method = nuisance_method,
                    working_model = working_model,
                    g_rct = g_rct,
                    num_covs = num_covs,
                    method = "tmle")
saveRDS(tmle_res,
        file = "out/tmle_" %+% f_name %+% "_" %+% date_name %+% ".RDS")

# 5. RCT only
rct_only_res <- run_sim(B = B,
                        rct_prop_min = rct_prop_min,
                        rct_prop_max = rct_prop_max,
                        rct_prop_step = rct_prop_step,
                        ate = ate,
                        bias = bias,
                        controls_only = controls_only,
                        nuisance_method = nuisance_method,
                        working_model = working_model,
                        g_rct = g_rct,
                        num_covs = num_covs,
                        method = "rct_only")
saveRDS(rct_only_res,
        file = "out/rct_only_" %+% f_name %+% "_" %+% date_name %+% ".RDS")

# save all environment variables
save.image(file = "out/env_" %+% f_name %+% "_" %+% date_name %+% ".RData")
