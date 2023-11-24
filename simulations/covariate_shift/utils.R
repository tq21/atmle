library(EScvtmle)
devtools::load_all()
source("sim_data.R")
`%+%` <- function(a, b) paste0(a, b)

run_sim <- function(ate,
                    g_rct,
                    bias,
                    gamma,
                    nuisance_method,
                    working_model,
                    B) {

  # results
  res <- data.frame(estimator = c("A-TMLE", "ES-CVTMLE", "TMLE"),
                    gamma = rep(gamma, 3),
                    bias = numeric(3),
                    variance = numeric(3),
                    mse = numeric(3),
                    coverage = numeric(3))

  atmle_ests <- numeric(B)
  atmle_cover <- numeric(B)
  escvtmle_ests <- numeric(B)
  escvtmle_cover <- numeric(B)
  tmle_ests <- numeric(B)
  tmle_cover <- numeric(B)

  for (i in 1:B) {
    print("iteration: " %+% i)

    # RCT
    data <- sim_data(ate, 1000, 0.2, g_rct, bias, FALSE, gamma)

    # atmle
    atmle_res <- atmle(data = data,
                       S_node = 1,
                       W_node = 2:5,
                       A_node = 6,
                       Y_node = 7,
                       controls_only = FALSE,
                       family = "gaussian",
                       atmle_pooled = TRUE,
                       var_method = "ic",
                       theta_method = nuisance_method,
                       Pi_method = nuisance_method,
                       g_method = nuisance_method,
                       theta_tilde_method = nuisance_method,
                       Q_method = nuisance_method,
                       bias_working_model = working_model,
                       pooled_working_model = working_model,
                       g_rct = g_rct,
                       verbose = FALSE)

    atmle_ests[i] <- atmle_res$est
    atmle_cover[i] <- (atmle_res$lower < ate) & (ate < atmle_res$upper)

    # escvtmle
    tmp <- ES.cvtmle(txinrwd = TRUE,
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
    escvtmle_res <- list(est = tmp$ATE$b2v,
                         lower = as.numeric(tmp$CI$b2v[1]),
                         upper = as.numeric(tmp$CI$b2v[2]))

    escvtmle_ests[i] <- escvtmle_res$est
    escvtmle_cover[i] <- (escvtmle_res$lower < ate) & (ate < escvtmle_res$upper)

    # TODO: add nonparametric TMLE
    tmle_res <- nonparametric(data = data,
                              S_node = 1,
                              W_node = 2:5,
                              A_node = 6,
                              Y_node = 7,
                              controls_only = FALSE,
                              family = "gaussian",
                              atmle_pooled = FALSE,
                              var_method = "ic",
                              theta_method = nuisance_method,
                              Pi_method = nuisance_method,
                              g_method = nuisance_method,
                              theta_tilde_method = nuisance_method,
                              Q_method = nuisance_method,
                              bias_working_model = working_model,
                              pooled_working_model = working_model,
                              g_rct = g_rct,
                              verbose = FALSE)

    tmle_ests[i] <- tmle_res$est
    tmle_cover[i] <- (tmle_res$lower < ate) & (ate < tmle_res$upper)
  }

  res$bias[1] <- abs(mean(atmle_ests) - ate)
  res$bias[2] <- abs(mean(escvtmle_ests) - ate)
  res$bias[3] <- abs(mean(tmle_ests) - ate)
  res$variance[1] <- var(atmle_ests)
  res$variance[2] <- var(escvtmle_ests)
  res$variance[3] <- var(tmle_ests)
  res$mse[1] <- mean((atmle_ests - ate)^2)
  res$mse[2] <- mean((escvtmle_ests - ate)^2)
  res$mse[3] <- mean((tmle_ests - ate)^2)
  res$coverage[1] <- mean(atmle_cover)
  res$coverage[2] <- mean(escvtmle_cover)
  res$coverage[3] <- mean(tmle_cover)

  return(res)
}
