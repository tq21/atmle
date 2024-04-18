#' Simulations comparing the following estimators:
#' 1. A-TMLE for psi_tilde
#' 2. TMLE for psi_tilde
#' Data generating distribution scenario: severe positivity violation

source("utils.R")
source("utils_positivity.R")

# simulation parameters --------------------------------------------------------
B <- 500 # number of runs for each sample size n
n_min <- 500 # smallest sample size
n_max <- 5000 # largest sample size
n_step <- 500 # sample size increment
bA <- 1.5 # true ATE
bias <- 0 # true bias
nuisance_method <- "glm"
working_model <- "lasso"
pRCT <- 0.67
f_name <- "psi_tilde_positivity"
date_name <- "922"

# 1. A-TMLE for psi_tilde
atmle_res <- run_sim_n_increase(B = B,
                                n_min = n_min,
                                n_max = n_max,
                                n_step = n_step,
                                bA = bA,
                                bias = bias,
                                pRCT = pRCT,
                                nuisance_method = nuisance_method,
                                working_model = working_model,
                                verbose = TRUE,
                                method = "atmle psi_tilde")
saveRDS(atmle_res,
        file = "out/atmle_" %+% f_name %+% "_" %+% date_name %+% ".RDS")

# 2. TMLE for psi_tilde
tmle_res <- run_sim_n_increase(B = B,
                               n_min = n_min,
                               n_max = n_max,
                               n_step = n_step,
                               bA = bA,
                               bias = bias,
                               pRCT = pRCT,
                               nuisance_method = nuisance_method,
                               working_model = working_model,
                               verbose = TRUE,
                               method = "tmle psi_tilde")
saveRDS(tmle_res,
        file = "out/tmle_" %+% f_name %+% "_" %+% date_name %+% ".RDS")
