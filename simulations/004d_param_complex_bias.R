#' Simulations comparing the following estimators:
#' 1. A-TMLE both psi_pound and psi_tilde
#' 2. A-TMLE psi_pound, regular TMLE psi_tilde
#' 3. ESCVTMLE
#' 4. Nonparametric TMLE
#' Data generating distribution scenario: complex parametric bias

source("utils.R")
source("sim_data.R")

# simulation parameters --------------------------------------------------------
B <- 1000 # number of runs for each sample size n
n_min <- 1000 # smallest sample size
n_max <- 3000 # largest sample size
n_step <- 500 # sample size increment
bA <- 1.5 # true ATE
bias <- "param_complex" # true bias
nuisance_method <- "glm"
working_model <- "glmnet"
g_rct <- 0.67
f_name <- "param_complex_glm_glmnet"
date_name <- "1028"
controls_only <- FALSE
num_covs <- 4
var_method <- "ic"
family <- "gaussian"

# 1. A-TMLE, both psi_pound and psi_tilde
atmle_both_res <- run_sim_n_increase(B = B,
                                     n_min = n_min,
                                     n_max = n_max,
                                     n_step = n_step,
                                     bA = bA,
                                     bias = bias,
                                     g_rct = g_rct,
                                     num_covs = num_covs,
                                     var_method = var_method,
                                     controls_only = controls_only,
                                     nuisance_method = nuisance_method,
                                     working_model = working_model,
                                     verbose = TRUE,
                                     family = family,
                                     method = "atmle")
saveRDS(atmle_both_res,
        file = "out/atmle_both_" %+% f_name %+% "_" %+% date_name %+% ".RDS")

# 2. Oracle A-TMLE, both psi_pound and psi_tilde
# oracle_atmle_both_res <- run_sim_n_increase(B = B,
#                                             n_min = n_min,
#                                             n_max = n_max,
#                                             n_step = n_step,
#                                             bA = bA,
#                                             bias = bias,
#                                             g_rct = g_rct,
#                                             num_covs = num_covs,
#                                             var_method = var_method,
#                                             controls_only = controls_only,
#                                             nuisance_method = nuisance_method,
#                                             working_model = working_model,
#                                             verbose = TRUE,
#                                             method = "oracle_atmle")
# saveRDS(oracle_atmle_both_res,
#         file = "out/oracle_atmle_both_" %+% f_name %+% "_" %+% date_name %+% ".RDS")

# 3. A-TMLE psi_pound, regular TMLE psi_tilde
# atmle_tmle_res <- run_sim_n_increase(B = B,
#                                      n_min = n_min,
#                                      n_max = n_max,
#                                      n_step = n_step,
#                                      bA = bA,
#                                      bias = bias,
#                                      g_rct = g_rct,
#                                      num_covs = num_covs,
#                                      var_method = var_method,
#                                      controls_only = controls_only,
#                                      nuisance_method = nuisance_method,
#                                      working_model = working_model,
#                                      family = family,
#                                      verbose = TRUE,
#                                      method = "atmle_tmle")
# saveRDS(atmle_tmle_res,
#         file = "out/atmle_tmle_" %+% f_name %+% "_" %+% date_name %+% ".RDS")

# 4. Oracle A-TMLE psi_pound, regular TMLE psi_tilde
# oracle_atmle_tmle_res <- run_sim_n_increase(B = B,
#                                             n_min = n_min,
#                                             n_max = n_max,
#                                             n_step = n_step,
#                                             bA = bA,
#                                             bias = bias,
#                                             g_rct = g_rct,
#                                             num_covs = num_covs,
#                                             var_method = var_method,
#                                             controls_only = controls_only,
#                                             nuisance_method = nuisance_method,
#                                             working_model = working_model,
#                                             verbose = TRUE,
#                                             method = "oracle_atmle_tmle")
# saveRDS(oracle_atmle_tmle_res,
#         file = "out/oracle_atmle_tmle_" %+% f_name %+% "_" %+% date_name %+% ".RDS")

# 5. ESCVTMLE
escvtmle_res <- run_sim_n_increase(B = B,
                                   n_min = n_min,
                                   n_max = n_max,
                                   n_step = n_step,
                                   bA = bA,
                                   bias = bias,
                                   g_rct = g_rct,
                                   num_covs = num_covs,
                                   var_method = var_method,
                                   controls_only = controls_only,
                                   nuisance_method = nuisance_method,
                                   working_model = working_model,
                                   family = family,
                                   verbose = TRUE,
                                   method = "escvtmle")
saveRDS(escvtmle_res,
        file = "out/escvtmle_" %+% f_name %+% "_" %+% date_name %+% ".RDS")

# 6. Nonparametric TMLE
tmle_res <- run_sim_n_increase(B = B,
                               n_min = n_min,
                               n_max = n_max,
                               n_step = n_step,
                               bA = bA,
                               bias = bias,
                               g_rct = g_rct,
                               num_covs = num_covs,
                               var_method = var_method,
                               controls_only = controls_only,
                               nuisance_method = nuisance_method,
                               working_model = working_model,
                               family = family,
                               verbose = TRUE,
                               method = "tmle")
saveRDS(tmle_res,
        file = "out/tmle_" %+% f_name %+% "_" %+% date_name %+% ".RDS")

# 7. RCT only
rct_only_res <- run_sim_n_increase(B = B,
                                   n_min = n_min,
                                   n_max = n_max,
                                   n_step = n_step,
                                   bA = bA,
                                   bias = bias,
                                   g_rct = g_rct,
                                   num_covs = num_covs,
                                   var_method = var_method,
                                   controls_only = controls_only,
                                   nuisance_method = nuisance_method,
                                   working_model = working_model,
                                   family = family,
                                   verbose = TRUE,
                                   method = "rct_only")
saveRDS(rct_only_res,
        file = "out/rct_only_" %+% f_name %+% "_" %+% date_name %+% ".RDS")
