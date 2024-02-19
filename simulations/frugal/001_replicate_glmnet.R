library(EScvtmle)
library(devtools)
load_all()

`%+%` <- function(a, b) paste0(a, b)

# DGP: SCM VERSION
sim_data <- function(n, gamma, rct) {
  W1 <- rnorm(n, 0, 1)
  W2 <- rbinom(n, 1, 0.5)
  W3 <- rnorm(n, 0.1*W2, 1)
  U1 <- rnorm(n, 0, 1)
  U2 <- rnorm(n, 0, 0.35)
  A <- numeric(length = n)

  if (rct) {
    S <- rep(1, n)
    A <- rbinom(n, 1, 0.5)
  } else {
    S <- rep(0, n)
    A <- rbinom(n, 1, plogis(W1+W2+gamma*U1+gamma*U2))
  }

  Y <- rnorm(n, (0.2+0.1*W3)*A+U1+U2, 1)

  data <- data.frame(S = S,
                     W1 = W1,
                     W2 = W2,
                     W3 = W3,
                     A = A,
                     Y = Y)

  return(data)
}

# SCM SIMULATIONS
set.seed(124)
B <- 200
gamma_seq <- seq(0, 5, 0.2)
ate <- 0.205
g_rct <- 0.5
n_rct <- 300
n_rwd <- 1200
controls_only <- FALSE
S_node <- 1
W_node <- 2:4
A_node <- 5
Y_node <- 6

est_pooled <- matrix(nrow = B, ncol = length(gamma_seq))
est_atmle <- matrix(nrow = B, ncol = length(gamma_seq))
est_escvtmle <- matrix(nrow = B, ncol = length(gamma_seq))
est_rct_only <- matrix(nrow = B, ncol = length(gamma_seq))
est_tmle <- matrix(nrow = B, ncol = length(gamma_seq))
cover_pooled <- matrix(nrow = B, ncol = length(gamma_seq))
cover_atmle <- matrix(nrow = B, ncol = length(gamma_seq))
cover_escvtmle <- matrix(nrow = B, ncol = length(gamma_seq))
cover_rct_only <- matrix(nrow = B, ncol = length(gamma_seq))
cover_tmle <- matrix(nrow = B, ncol = length(gamma_seq))

for (b in 1:B) {
  print("run: " %+% b)

  # simulate data
  rct_data_list <- lapply(gamma_seq, function(bias) {
    return(sim_data(n_rct, bias, TRUE))
  })
  rwd_data_list <- lapply(gamma_seq, function(bias) {
    return(sim_data(n_rwd, bias, FALSE))
  })

  for (i in 1:length(gamma_seq)) {
    data <- rbind(rct_data_list[[i]], rwd_data_list[[i]])

    # tmle on pooled
    tmle_pooled_res <- tmle(Y = data[, Y_node],
                            A = data[, A_node],
                            W = data[, W_node],
                            Qform = "Y ~ A + A*W3",
                            g.SL.library = c("SL.glm"),
                            Q.SL.library = c("SL.glm"),
                            family = "gaussian")

    # rct only
    rct_only_res <- rct_only(data = data,
                             S_node = S_node,
                             W_node = W_node,
                             A_node = A_node,
                             Y_node = Y_node,
                             g_rct = g_rct,
                             nuisance_method = "glm",
                             family = "gaussian",
                             verbose = FALSE)

    # tmle on RWD
    nonparam_tmle <- nonparametric(data = data,
                                   S_node = S_node,
                                   W_node = W_node,
                                   A_node = A_node,
                                   Y_node = Y_node,
                                   controls_only = controls_only,
                                   family = "gaussian",
                                   atmle_pooled = TRUE,
                                   theta_method = "glm",
                                   Pi_method = "glm",
                                   g_method = "glm",
                                   theta_tilde_method = "glm",
                                   Q_method = "glm",
                                   bias_working_model = "glmnet",
                                   pooled_working_model = "glmnet",
                                   g_rct = g_rct,
                                   verbose = FALSE)

    # atmle with HAL
    atmle_res <- atmle(data = data,
                       S_node = S_node,
                       W_node = W_node,
                       A_node = A_node,
                       Y_node = Y_node,
                       atmle_pooled = TRUE,
                       controls_only = controls_only,
                       theta_method = "glm",
                       Pi_method = "glm",
                       g_method = "glm",
                       theta_tilde_method = "glm",
                       Q_method = "glm",
                       bias_working_model = "glmnet",
                       pooled_working_model = "glmnet",
                       g_rct = g_rct,
                       family = "gaussian",
                       min_working_model = TRUE,
                       verbose = FALSE)

    # ES-CVTMLE
    escvtmle_res <- ES.cvtmle(txinrwd = !controls_only,
                              data = data,
                              study = "S",
                              covariates = c("W1", "W2", "W3"),
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

    if (atmle_res$lower <= ate & atmle_res$upper >= ate) {
      print("covered")
    } else {
      print("not covered")
    }

    # collect results
    est_pooled[b,i] <- tmle_pooled_res$estimates$ATE$psi
    est_atmle[b,i] <- atmle_res$est
    est_rct_only[b,i] <- rct_only_res$est
    est_escvtmle[b,i] <- escvtmle_res$ATE$b2v
    est_tmle[b,i] <- nonparam_tmle$est
    cover_pooled[b,i] <- tmle_pooled_res$estimates$ATE$CI[1] <= ate &  tmle_pooled_res$estimates$ATE$CI[2]  >= ate
    cover_atmle[b,i] <- atmle_res$lower <= ate & atmle_res$upper >= ate
    cover_rct_only[b,i] <- rct_only_res$lower <= ate & rct_only_res$upper >= ate
    cover_escvtmle[b,i] <- escvtmle_res$CI$b2v[1] <= ate & escvtmle_res$CI$b2v[2] >= ate
    cover_tmle[b,i] <- nonparam_tmle$lower <= ate & nonparam_tmle$upper >= ate
  }
}

save(list = c("est_pooled",
              "est_atmle",
              "est_escvtmle",
              "est_rct_only",
              "est_tmle",
              "cover_pooled",
              "cover_atmle",
              "cover_escvtmle",
              "cover_rct_only",
              "cover_tmle"), file = "out/001_replicate_glmnet.RData")
