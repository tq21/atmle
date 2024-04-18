# source("utils.R")
# source("sim_data.R")
options(sl3.verbose = TRUE)

g_rct <- 0.67
bias <- "HAL_3"
data <- sim_two_covs(1.5, 1000, 0.2, g_rct, bias, FALSE)

S_node <- 1
W_node <- 2:3
A_node <- 4
Y_node <- 5
controls_only <- FALSE
var_method <- "ic"
theta_method <- "glm"
Pi_method <- "glm"
g_method <- "glm"
theta_tilde_method <- "glm"
Q_method <- "glmnet"
bias_working_model <- "HAL"
pooled_working_model <- "glmnet"
verbose <- TRUE
family <- "gaussian"

atmle_res <- atmle(data = data,
                   S_node = S_node,
                   W_node = W_node,
                   A_node = A_node,
                   Y_node = Y_node,
                   atmle_pooled = TRUE,
                   controls_only = controls_only,
                   theta_method = theta_method,
                   Pi_method = Pi_method,
                   g_method = g_method,
                   theta_tilde_method = theta_tilde_method,
                   Q_method = Q_method,
                   bias_working_model = bias_working_model,
                   pooled_working_model = pooled_working_model,
                   g_rct = g_rct,
                   family = family,
                   verbose = FALSE,
                   min_working_model = TRUE,
                   max_iter = 1)

escvtmle_res <- ES.cvtmle(txinrwd = !controls_only,
                          data = data,
                          study = "S",
                          covariates = c("W1", "W2", "W3"),
                          treatment_var = "A",
                          treatment = 1,
                          outcome = "Y",
                          pRCT = g_rct,
                          family = family,
                          Q.SL.library = c("SL.glm"),
                          g.SL.library = c("SL.glm"),
                          Q.discreteSL = TRUE,
                          g.discreteSL = TRUE,
                          V = 5)

# B <- 500
# n <- 2000
# ate <- 1.5
# bias <- "param_simple"
# nuisance_method <- "glmnet"
# working_model <- "glmnet"
# g_rct <- 0.67
# verbose <- TRUE
# controls_only <- FALSE
# num_covs <- 4
#
# tmp_1 <- run_sim(B = B,
#                  n = n,
#                  ate = ate,
#                  bias = bias,
#                  nuisance_method = nuisance_method,
#                  working_model = working_model,
#                  g_rct = g_rct,
#                  controls_only = controls_only,
#                  var_method = "ic",
#                  num_covs = num_covs,
#                  verbose = verbose,
#                  family = "gaussian",
#                  method = "atmle",
#                  type = "NULL")
# mean(tmp_1$psi_coverage)
# hist(tmp_1$psi_est)
# var(tmp_1$psi_est)+(mean(tmp_1$psi_est)-ate)^2
#
# tmp_2 <- run_sim(B = B,
#                  n = n,
#                  ate = ate,
#                  bias = bias,
#                  nuisance_method = nuisance_method,
#                  working_model = working_model,
#                  g_rct = g_rct,
#                  controls_only = controls_only,
#                  var_method = "ic",
#                  num_covs = num_covs,
#                  verbose = verbose,
#                  family = "gaussian",
#                  method = "escvtmle",
#                  type = "NULL")
# mean(tmp_2$psi_coverage)
# hist(tmp_2$psi_est)
# var(tmp_2$psi_est)+(mean(tmp_2$psi_est)-ate)^2

atmle_res$upper-atmle_res$lower
#atmle_res_bounded$upper-atmle_res_bounded$lower
as.numeric(escvtmle_res$CI$b2v[2]-escvtmle_res$CI$b2v[1])
# tmle_res$upper-tmle_res$lower
# rct_only_res$upper-rct_only_res$lower
escvtmle_res$proportionselected
#rct_only_res$upper-rct_only_res$lower
