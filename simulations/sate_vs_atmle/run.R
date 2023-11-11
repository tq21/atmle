source("sim_data.R")
devtools::load_all()
`%+%` <- function(a, b) paste0(a, b)

# parameters
ate <- 1.5
B <- 100
n <- 1000
g_rct_seq <- 0.67
bias <- "param_simple"

res <- data.frame(g_rct = g_rct_seq,
                  rct_prop = g_rct_seq,
                  sate_mse = numeric(length = length(g_rct_seq)),
                  sate_cover = numeric(length = length(g_rct_seq)),
                  atmle_mse = numeric(length = length(g_rct_seq)),
                  atmle_cover = numeric(length = length(g_rct_seq)))

for (i in 1:length(g_rct_seq)) {
  sate_est <- numeric(length = B)
  atmle_est <- numeric(length = B)
  sate_cover <- numeric(length = B)
  atmle_cover <- numeric(length = B)
  escvtmle_est <- numeric(length = B)
  escvtmle_cover <- numeric(length = B)
  print("g_rct = " %+% g_rct_seq[i])

  for (j in 1:B) {
    print("B = " %+% j %+% " of " %+% B %+% " for g_rct = " %+% g_rct_seq[i])
    data <- sim_four_covs(ate, n, 0.2, g_rct_seq[i], bias, FALSE)

    # sample average treatment effect
    sate_res <- sate(data = data,
                     S_node = 1,
                     W_node = c(2, 3, 4, 5),
                     A_node = 6,
                     Y_node = 7,
                     g_rct = g_rct_seq[i],
                     family = "gaussian")

    # atmle
    atmle_res <- atmle(data = data,
                       S_node = 1,
                       W_node = c(2, 3, 4, 5),
                       A_node = 6,
                       Y_node = 7,
                       controls_only = FALSE,
                       family = "gaussian",
                       atmle_pooled = TRUE,
                       r_loss = TRUE,
                       theta_method = "glm",
                       Pi_method = "glm",
                       g_method = "glm",
                       theta_tilde_method = "glm",
                       Q_method = "glm",
                       bias_working_model = "glmnet",
                       pooled_working_model = "glmnet",
                       g_rct = g_rct_seq[i],
                       verbose = FALSE)

    # es-cvtmle
    escvtmle_res <- ES.cvtmle(txinrwd = TRUE,
                              data = data,
                              study = "S",
                              covariates = c("W1", "W2", "W3", "W4"),
                              treatment_var = "A",
                              treatment = 1,
                              outcome = "Y",
                              pRCT = g_rct_seq[i],
                              family = "gaussian",
                              Q.SL.library = c("SL.glm"),
                              g.SL.library = c("SL.glm"),
                              Q.discreteSL = TRUE,
                              g.discreteSL = TRUE,
                              target.gwt = TRUE,
                              V = 5)

    # store results
    sate_est[j] <- sate_res$est
    atmle_est[j] <- atmle_res$est
    escvtmle_est[j] <- escvtmle_res$ATE$b2v
    sate_cover[j] <- ifelse(ate >= sate_res$lower & ate <= sate_res$upper, 1, 0)
    atmle_cover[j] <- ifelse(ate >= atmle_res$lower & ate <= atmle_res$upper, 1, 0)
    escvtmle_cover[j] <- ifelse(ate >= escvtmle_res$CI$b2v[1] & ate <= escvtmle_res$CI$b2v[2], 1, 0)
  }

  res$sate_mse[i] <- mean((sate_est - ate)^2)
  res$sate_cover[i] <- mean(sate_cover)
  res$atmle_mse[i] <- mean((atmle_est - ate)^2)
  res$atmle_cover[i] <- mean(atmle_cover)
  res$escvtmle_mse[i] <- mean((escvtmle_est - ate)^2)
  res$escvtmle_cover[i] <- mean(escvtmle_cover)
}

