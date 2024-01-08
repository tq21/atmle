options(sl3.verbose = TRUE)
source("utils.R")
source("sim_data.R")
set.seed(29857)

B <- 200
n <- 2000
ate <- 1.5
bias <- "small_bias"
nuisance_method <- "glm"
working_model <- "glmnet"
g_rct <- 0.67
verbose <- TRUE
controls_only <- FALSE
num_covs <- 4

tmp_1 <- run_sim(B = B,
                 n = n,
                 ate = ate,
                 bias = bias,
                 nuisance_method = nuisance_method,
                 working_model = working_model,
                 Q_method = "sl3",
                 g_rct = g_rct,
                 controls_only = controls_only,
                 var_method = "ic",
                 num_covs = num_covs,
                 verbose = verbose,
                 family = "gaussian",
                 method = "atmle_tmle",
                 type = "NULL")
mean(tmp_1$psi_coverage)
var(tmp_1$psi_est)+(mean(tmp_1$psi_est)-1.5)^2

tmp_2 <- run_sim(B = B,
                 n = n,
                 ate = ate,
                 bias = bias,
                 nuisance_method = nuisance_method,
                 working_model = working_model,
                 g_rct = g_rct,
                 controls_only = controls_only,
                 var_method = "ic",
                 num_covs = num_covs,
                 verbose = verbose,
                 family = "gaussian",
                 method = "escvtmle",
                 type = "NULL")
mean(tmp_2$psi_coverage)
var(tmp$psi_est)
hist(tmp$psi_est)
mean(tmp$psi_est)-1.5
var(tmp_2$psi_est)+(mean(tmp_2$psi_est)-ate)^2
tmp$psi_ci_upper-tmp$psi_ci_lower
