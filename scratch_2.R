source("utils.R")

large_n <- 10000000
W1 <- rnorm(large_n, 0, 1)
W2 <- rnorm(large_n, 0, 1)
W3 <- rnorm(large_n, 0, 1)
W4 <- rnorm(large_n, 0, 1)
ate <- mean(plogis(1.5+0.8*W1-1.1*W2+0.9*W3-1.3*W4)-plogis(-2+0.8*W1-1.1*W2+0.9*W3-1.3*W4))

data <- sim_binary_outcome(1.5, 5000, 0.2, 0.67, "param_complex", FALSE)

S_node = 1
W_node = c(2, 3, 4, 5)
A_node = 6
Y_node = 7
nuisance_method = "glm"
working_model = "glmnet"
g_rct = 0.67
verbose = TRUE
controls_only = FALSE

# atmle_oracle_res <- atmle_oracle(data,
#                                  S_node = S_node,
#                                  W_node = W_node,
#                                  A_node = A_node,
#                                  Y_node = Y_node,
#                                  controls_only = controls_only,
#                                  bias = bias,
#                                  atmle_pooled = FALSE,
#                                  var_method = "ic",
#                                  nuisance_method = nuisance_method,
#                                  working_model = working_model,
#                                  g_rct = g_rct,
#                                  verbose = verbose)

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

tmle_res <- nonparametric(data = data,
                          S_node = S_node,
                          W_node = W_node,
                          A_node = A_node,
                          Y_node = Y_node,
                          g_rct = g_rct,
                          nuisance_method = nuisance_method,
                          working_model = working_model,
                          verbose = verbose)

rct_only_res <- rct_only(data = data,
                         S_node = S_node,
                         W_node = W_node,
                         A_node = A_node,
                         Y_node = Y_node,
                         g_rct = g_rct,
                         nuisance_method = nuisance_method,
                         verbose = verbose)

#atmle_oracle_res$upper-atmle_oracle_res$lower

atmle_res$upper-atmle_res$lower
as.numeric(escvtmle_res$CI$b2v[2]-escvtmle_res$CI$b2v[1])
tmle_res$upper-tmle_res$lower
rct_only_res$upper-rct_only_res$lower
escvtmle_res$proportionselected
ate


