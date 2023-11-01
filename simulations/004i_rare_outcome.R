#' Simulations comparing the following estimators:
#' 1. A-TMLE both psi_pound and psi_tilde
#' 2. A-TMLE psi_pound, regular TMLE psi_tilde
#' 3. ESCVTMLE
#' 4. Nonparametric TMLE
#' Data generating distribution scenario: complex parametric bias

source("utils.R")
source("sim_data.R")

# simulation parameters --------------------------------------------------------
B <- 200 # number of runs for each sample size n
n_min <- 1000 # smallest sample size
n_max <- 1500 # largest sample size
n_step <- 500 # sample size increment

# monte-carlo to evaluate the true ATE
Y1 <- sim_rare_outcomes(2.1, 5000000, 1, 1, 0, FALSE)
Y0 <- sim_rare_outcomes(2.1, 5000000, 1, 0, 0, FALSE)
bA <- mean(Y1$Y) - mean(Y0$Y) # true ATE
bias <- "param_complex"
nuisance_method <- "glm"
working_model <- "glmnet"
g_rct <- 0.67
f_name <- "param_complex_glm_glmnet_rare"
date_name <- "1031"
controls_only <- FALSE
num_covs <- 4
var_method <- "ic"
family <- "binomial"
type <- "rare"

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
                                     type = type,
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

# 2. ES-CVTMLE
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
                                   verbose = FALSE,
                                   type = type,
                                   method = "escvtmle")
saveRDS(escvtmle_res,
        file = "out/escvtmle_" %+% f_name %+% "_" %+% date_name %+% ".RDS")

# 3. Nonparametric TMLE
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
                               verbose = FALSE,
                               method = "tmle")
saveRDS(tmle_res,
        file = "out/tmle_" %+% f_name %+% "_" %+% date_name %+% ".RDS")

# 4. RCT only
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
