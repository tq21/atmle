source("utils.R")
source("sim_data.R")
options(sl3.verbose = TRUE)

g_rct <- 0.67
data <- sim_four_covs(1.5, 2000, 0.2, g_rct, 1.8, FALSE)

S_node <- 1
W_node <- c(2, 3, 4, 5)
A_node <- 6
Y_node <- 7
controls_only <- FALSE
var_method <- "ic"
theta_method <- "glm"
Pi_method <- "glm"
g_method <- "glm"
theta_tilde_method <- "glm"
Q_method <- "glm"
bias_working_model <- "glmnet"
pooled_working_model <- "glmnet"
g_rct <- 0.67
verbose <- TRUE
family <- "gaussian"

atmle_res <- atmle(data = data,
                   S_node = S_node,
                   W_node = W_node,
                   A_node = A_node,
                   Y_node = Y_node,
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

escvtmle_res <- ES.cvtmle(txinrwd = !controls_only,
                          data = data,
                          study = "S",
                          covariates = c("W1", "W2", "W3", "W4"),
                          treatment_var = "A",
                          treatment = 1,
                          outcome = "Y",
                          pRCT = g_rct,
                          family = "gaussian",
                          Q.SL.library = c("SL.glm"),
                          g.SL.library = c("SL.glm"),
                          Q.discreteSL = TRUE,
                          g.discreteSL = TRUE,
                          V = 5)

atmle_res$upper-atmle_res$lower
as.numeric(escvtmle_res$CI$b2v[2]-escvtmle_res$CI$b2v[1])
escvtmle_res$proportionselected
