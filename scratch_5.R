#source("utils.R")
source("sim_data.R")

B <- 1000
cover <- numeric(length = B)
est <- numeric(length = B)

for (i in 1:B) {
  print(i)
  data <- sim_four_covs(1.5, 2000, 0.2, 0.67, 0, FALSE)

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
  verbose <- FALSE
  atmle_pooled <- TRUE
  v_folds <- 5

  escvtmle_res <- ES.cvtmle(txinrwd = TRUE,
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

  if (escvtmle_res$CI$b2v[2] >= 1.5 & escvtmle_res$CI$b2v[1] <= 1.5) {
    cover[i] <- 1
    print("covered")
  } else {
    cover[i] <- 0
    print("not covered")
  }

  est[i] <- escvtmle_res$ATE$b2v
}
