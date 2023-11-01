source("utils.R")
source("sim_data.R")
options(sl3.verbose = TRUE)

g_rct <- 0.67
data <- sim_rare_outcomes(2.1, 2000, 0.2, g_rct, "param_complex", TRUE)
# mean(data$Y)
# mean(data$Y[data$S == 1])
# mean(data$Y[data$S == 0])

Y1 <- sim_rare_outcomes(2.1, 5000000, 1, 1, 0, FALSE)
Y0 <- sim_rare_outcomes(2.1, 5000000, 1, 0, 0, FALSE)
truth <- mean(Y1$Y) - mean(Y0$Y)

S_node <- 1
W_node <- c(2, 3, 4, 5)
A_node <- 6
Y_node <- 7
controls_only <- TRUE
var_method <- "ic"
theta_method <- "glm"
Pi_method <- "glm"
g_method <- "glm"
theta_tilde_method <- "glm"
Q_method <- "glm"
bias_working_model <- "glmnet"
pooled_working_model <- "glmnet"
verbose <- TRUE
family <- "binomial"

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

atmle_res_bounded <- atmle(data = data,
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
                           g_bound = c(0.025, 0.975),
                           Pi_bound = c(0.025, 0.975),
                           theta_bound = c(0, 1),
                           verbose = FALSE)

escvtmle_res <- ES.cvtmle(txinrwd = !controls_only,
                          data = data,
                          study = "S",
                          covariates = c("W1", "W2", "W3", "W4"),
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



# tmle_res <- nonparametric(data = data,
#                           S_node = S_node,
#                           W_node = W_node,
#                           A_node = A_node,
#                           Y_node = Y_node,
#                           g_rct = g_rct,
#                           nuisance_method = nuisance_method,
#                           working_model = working_model,
#                           verbose = verbose)
#
# rct_only_res <- rct_only(data = data,
#                          S_node = S_node,
#                          W_node = W_node,
#                          A_node = A_node,
#                          Y_node = Y_node,
#                          g_rct = g_rct,
#                          nuisance_method = nuisance_method,
#                          verbose = verbose)

#atmle_oracle_res$upper-atmle_oracle_res$lower

atmle_res$upper-atmle_res$lower
atmle_res_bounded$upper-atmle_res_bounded$lower
as.numeric(escvtmle_res$CI$b2v[2]-escvtmle_res$CI$b2v[1])
# tmle_res$upper-tmle_res$lower
# rct_only_res$upper-rct_only_res$lower
escvtmle_res$proportionselected
#rct_only_res$upper-rct_only_res$lower
