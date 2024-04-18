source("utils.R")
source("sim_data.R")
set.seed(29857)

B <- 200
n <- 1000
ate <- 1.5
bias <- 1.8
nuisance_method <- "glm"
working_model <- "glmnet"
g_rct <- 0.67
verbose <- TRUE
controls_only <- TRUE
num_covs <- 4

tmp_1 <- run_sim(B = B,
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
                 method = "atmle",
                 type = "NULL")

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
mean(tmp_1$psi_coverage)
var(tmp_1$psi_est)+(mean(tmp_1$psi_est)-1.5)^2

mean(tmp_2$escvtmle_prop_selected)
mean(tmp_2$psi_coverage)
var(tmp_2$psi_est)+(mean(tmp_2$psi_est)-ate)^2
