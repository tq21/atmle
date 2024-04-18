library(ggplot2)
library(devtools)
library(EScvtmle)
library(sl3)
library(purrr)

load_all()
`%+%` <- function(a, b) paste0(a, b)
source("sim_data.R")

get_truth <- function() {
  data_A1 <- sim_data(beta = beta,
                      n = 100000,
                      rct = TRUE,
                      g_rct = 0,
                      bias = bias,
                      A_counter = 1)
  data_A0 <- sim_data(beta = beta,
                      n = 100000,
                      rct = TRUE,
                      g_rct = 0,
                      bias = bias,
                      A_counter = 0)

  return(mean(data_A1$Y - data_A0$Y))
}

make_data <- function(B,
                      n_rct_seq,
                      n_rwd_seq,
                      beta,
                      bias) {

  data_list <- map2(n_rct_seq, n_rwd_seq, function(n_rct, n_rwd) {
    map(1:B, function(i) {
      rct_data <- sim_data(beta = beta,
                           n = n_rct,
                           rct = TRUE,
                           g_rct = g_rct,
                           bias = bias)
      rwd_data <- sim_data(beta = beta,
                           n = n_rwd,
                           rct = FALSE,
                           g_rct = g_rct,
                           bias = bias)
      rbind(rct_data, rwd_data)
    })
  })

  return(data_list)
}

run_sim <- function(data_list,
                    beta = beta,
                    bias = bias,
                    nuisance_method,
                    working_model,
                    g_rct,
                    method) {

  # results
  all_psi_est <- vector(mode = "list", length = length(data_list))
  all_psi_coverage <- vector(mode = "list", length = length(data_list))
  all_psi_ci_lower <- vector(mode = "list", length = length(data_list))
  all_psi_ci_upper <- vector(mode = "list", length = length(data_list))
  all_escvtmle_prop_selected <- vector(mode = "list", length = length(data_list))

  # evaluate truth
  ate <- get_truth()

  for (i in 1:length(data_list)) {
    cur_data_list <- data_list[[i]]

    print("rct: " %+% sum(cur_data_list[[1]]$S) %+%
            ", rwd: " %+% sum(cur_data_list[[1]]$S == 0))

    psi_est <- vector(length = length(cur_data_list))
    psi_coverage <- vector(length = length(cur_data_list))
    psi_ci_lower <- vector(length = length(cur_data_list))
    psi_ci_upper <- vector(length = length(cur_data_list))
    escvtmle_prop_selected <- vector(length = length(cur_data_list))

    for (j in 1:length(cur_data_list)) {
      data <- cur_data_list[[j]]

      S_node <- 1
      W_node <- 2:4
      A_node <- 5
      Y_node <- 6

      # fit
      res <- NULL
      if (method == "atmle") {
        res <- atmle(data = data,
                     S_node = S_node,
                     W_node = W_node,
                     A_node = A_node,
                     Y_node = Y_node,
                     controls_only = TRUE,
                     family = "binomial",
                     atmle_pooled = TRUE,
                     theta_method = nuisance_method,
                     Pi_method = nuisance_method,
                     g_method = nuisance_method,
                     theta_tilde_method = nuisance_method,
                     Q_method = nuisance_method,
                     bias_working_model = working_model,
                     pooled_working_model = working_model,
                     g_rct = g_rct,
                     verbose = FALSE,
                     var_method = "ic",
                     enumerate_basis_args = list(
                       num_knots = hal9001:::num_knots_generator(
                         max_degree = 3,
                         smoothness_orders = 1,
                         base_num_knots_0 = 200,
                         base_num_knots_1 = 20)))
      } else if (method == "atmle_tmle") {
        res <- atmle(data = data,
                     S_node = S_node,
                     W_node = W_node,
                     A_node = A_node,
                     Y_node = Y_node,
                     controls_only = TRUE,
                     family = "binomial",
                     atmle_pooled = FALSE,
                     theta_method = nuisance_method,
                     Pi_method = nuisance_method,
                     g_method = nuisance_method,
                     theta_tilde_method = nuisance_method,
                     Q_method = nuisance_method,
                     bias_working_model = working_model,
                     pooled_working_model = working_model,
                     g_rct = g_rct,
                     verbose = FALSE,
                     enumerate_basis_args = list(
                       num_knots = hal9001:::num_knots_generator(
                         max_degree = 3,
                         smoothness_orders = 1,
                         base_num_knots_0 = 200,
                         base_num_knots_1 = 20)))
      } else if (method == "escvtmle") {
        covariates <- c("W1", "W2", "W3")
        tmp <- ES.cvtmle(txinrwd = FALSE,
                         data = data,
                         study = "S",
                         covariates = covariates,
                         treatment_var = "A",
                         treatment = 1,
                         outcome = "Y",
                         pRCT = g_rct,
                         family = "binomial",
                         Q.SL.library = c("SL.glm"),
                         g.SL.library = c("SL.glm"),
                         Q.discreteSL = TRUE,
                         g.discreteSL = TRUE,
                         V = 5)
        res <- list(est = tmp$ATE$b2v,
                    lower = as.numeric(tmp$CI$b2v[1]),
                    upper = as.numeric(tmp$CI$b2v[2]))
        escvtmle_prop_selected[j] <- tmp$proportionselected$b2v
      } else if (method == "tmle") {
        res <- nonparametric(data = data,
                             S_node = S_node,
                             W_node = W_node,
                             A_node = A_node,
                             Y_node = Y_node,
                             controls_only = TRUE,
                             family = "binomial",
                             atmle_pooled = TRUE,
                             theta_method = nuisance_method,
                             Pi_method = nuisance_method,
                             g_method = nuisance_method,
                             theta_tilde_method = nuisance_method,
                             Q_method = nuisance_method,
                             bias_working_model = working_model,
                             pooled_working_model = working_model,
                             g_rct = g_rct,
                             verbose = FALSE)
      } else if (method == "rct_only") {
        res <- rct_only(data,
                        S_node = S_node,
                        W_node = W_node,
                        A_node = A_node,
                        Y_node = Y_node,
                        nuisance_method = nuisance_method,
                        family = "binomial",
                        g_rct = g_rct,
                        verbose = FALSE)
      }

      if (res$lower <= ate & res$upper >= ate) {
        print("psi covered")
        psi_coverage[j] <- 1
      } else {
        print("psi not covered")
        psi_coverage[j] <- 0
      }

      print("current coverage: " %+% round(mean(psi_coverage[1:j]), 2))

      psi_est[j] <- res$est
      psi_ci_lower[j] <- res$lower
      psi_ci_upper[j] <- res$upper
    }

    all_psi_est[[i]] <- psi_est
    all_psi_coverage[[i]] <- psi_coverage
    all_psi_ci_lower[[i]] <- psi_ci_lower
    all_psi_ci_upper[[i]] <- psi_ci_upper

    if (method == "escvtmle") {
      all_escvtmle_prop_selected[[i]] <- escvtmle_prop_selected
    }
  }

  return(list(all_psi_est = all_psi_est,
              all_psi_coverage = all_psi_coverage,
              all_psi_ci_lower = all_psi_ci_lower,
              all_psi_ci_upper = all_psi_ci_upper,
              all_escvtmle_prop_selected = all_escvtmle_prop_selected))
}
