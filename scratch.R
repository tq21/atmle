library(devtools)
library(EScvtmle)
load_all()
source("sim_data.R")

g_rct <- 0.67
bias <- "b"
ate <- 1.5
n_rct <- 500
n_rwd <- 2000
controls_only <- FALSE
data_rct <- sim_data(ate = ate,
                     n = n_rct,
                     rct = TRUE,
                     g_rct = g_rct,
                     bias = bias,
                     controls_only = controls_only)
data_rwd <- sim_data(ate = ate,
                     n = n_rwd,
                     rct = FALSE,
                     g_rct = g_rct,
                     bias = bias,
                     controls_only = controls_only)
data <- rbind(data_rct, data_rwd)

S_node <- 1
W_node <- 2:4
A_node <- 5
Y_node <- 6
theta_method <- "glm"
Pi_method <- "glm"
g_method <- "glm"
theta_tilde_method <- "glm"
Q_method <- "glm"
bias_working_model <- "glmnet"
pooled_working_model <- "glmnet"
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
                   verbose = FALSE)

procova_est <- procova(data = data,
                       S_node = S_node,
                       W_node = W_node,
                       A_node = A_node,
                       Y_node = Y_node,
                       controls_only = controls_only,
                       family = family,
                       g_rct = g_rct,
                       Q_method = Q_method,
                       g_method = g_method,
                       v_folds = 5,
                       verbose = FALSE)

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

tmle_res <- nonparametric(data = data,
                          S_node = S_node,
                          W_node = W_node,
                          A_node = A_node,
                          Y_node = Y_node,
                          controls_only = controls_only,
                          family = "gaussian",
                          atmle_pooled = TRUE,
                          theta_method = theta_method,
                          Pi_method = Pi_method,
                          g_method = g_method,
                          theta_tilde_method = theta_tilde_method,
                          Q_method = Q_method,
                          bias_working_model = bias_working_model,
                          pooled_working_model = pooled_working_model,
                          g_rct = g_rct,
                          verbose = FALSE)
# (mean(atmle_both_res$all_psi_est[[1]])-ate)^2+var(atmle_both_res$all_psi_est[[1]])
# (mean(escvtmle_res$all_psi_est[[1]])-ate)^2+var(escvtmle_res$all_psi_est[[1]])
# (mean(tmle_res$all_psi_est[[1]])-ate)^2+var(tmle_res$all_psi_est[[1]])

# mean(atmle_both_res$all_psi_ci_upper[[1]]-atmle_both_res$all_psi_ci_lower[[1]])
# mean(escvtmle_res$all_psi_ci_upper[[1]]-escvtmle_res$all_psi_ci_lower[[1]])
# mean(tmle_res$all_psi_ci_upper[[1]]-tmle_res$all_psi_ci_lower[[1]])
# rct_only_res$upper-rct_only_res$lower
#rct_only_res$upper-rct_only_res$lower

atmle_res$upper-atmle_res$lower
procova_est$upper-procova_est$lower
as.numeric(escvtmle_res$CI$b2v[2]-escvtmle_res$CI$b2v[1])
escvtmle_res$proportionselected
tmle_res$upper-tmle_res$lower
