source("utils.R")

set.seed(23223)

data <- generate_two_covs(1.5, 800, 0.2, g_rct = 0.67, bias = "HAL", FALSE)

S_node = 1
W_node = c(2, 3)
A_node = 4
Y_node = 5
nuisance_method = "glm"
working_model = "HAL"
g_rct = 0.67
verbose = TRUE
controls_only = FALSE

atmle_res <- atmle(data,
                   S_node = S_node,
                   W_node = W_node,
                   A_node = A_node,
                   Y_node = Y_node,
                   controls_only = controls_only,
                   atmle_pooled = TRUE,
                   var_method = "ic",
                   nuisance_method = nuisance_method,
                   working_model = working_model,
                   g_rct = g_rct,
                   verbose = verbose)

escvtmle_res <- ES.cvtmle(txinrwd = !controls_only,
                          data = data,
                          study = "S",
                          covariates = c("W1", "W2"),
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

tmle_res <- nonparametric(data = data,
                          S_node = S_node,
                          W_node = W_node,
                          A_node = A_node,
                          Y_node = Y_node,
                          g_rct = g_rct,
                          nuisance_method = nuisance_method,
                          working_model = working_model,
                          verbose = verbose)

atmle_res$upper-atmle_res$lower
escvtmle_res$psi_ci_upper-escvtmle_res$psi_ci_lower
tmp_3$psi_ci_upper-tmp_3$psi_ci_lower
tmp_4$psi_ci_upper-tmp_4$psi_ci_lower
tmp_2$escvtmle_prop_selected
