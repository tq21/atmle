#source("utils.R")
source("sim_data.R")

B <- 1000
cover <- numeric(length = B)
est <- numeric(length = B)

for (i in 1:B) {
  print(i)
  data <- sim_four_covs(1.5, 2000, 0.2, 0.67, "large_bias", TRUE)

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
  g_rct <- 0.67
  verbose <- FALSE
  atmle_pooled <- TRUE
  v_folds <- 5

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
                     family = "gaussian",
                     verbose = verbose)

  if (atmle_res$upper >= 1.5 & atmle_res$lower <= 1.5) {
    cover[i] <- 1
    print("covered")
  } else {
    cover[i] <- 0
    print("not covered")
  }

  est[i] <- atmle_res$est
}
