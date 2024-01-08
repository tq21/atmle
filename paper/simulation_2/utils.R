library(ggplot2)
library(devtools)
library(EScvtmle)
library(sl3)

load_all()
`%+%` <- function(a, b) paste0(a, b)

run_sim_n_increase <- function(B,
                               n_min,
                               n_max,
                               n_step,
                               rct_prop,
                               ate,
                               bias,
                               controls_only,
                               nuisance_method,
                               working_model,
                               g_rct,
                               num_covs,
                               method,
                               verbose = TRUE) {

  n_seq <- seq(n_min, n_max, n_step)

  # results
  all_psi_est <- vector(mode = "list", length = length(n_seq))
  all_psi_coverage <- vector(mode = "list", length = length(n_seq))
  all_psi_ci_lower <- vector(mode = "list", length = length(n_seq))
  all_psi_ci_upper <- vector(mode = "list", length = length(n_seq))
  all_escvtmle_prop_selected <- vector(mode = "list", length = length(n_seq))

  for (i in 1:length(n_seq)) {
    n <- n_seq[i]
    print(n)

    psi_est <- vector(length = B)
    psi_coverage <- vector(length = B)
    psi_ci_lower <- vector(length = B)
    psi_ci_upper <- vector(length = B)
    escvtmle_prop_selected <- vector(length = B)

    for (j in 1:B) {
      # simulate data
      data <- NULL
      S_node <- 1
      W_nodes <- NULL
      if (num_covs == 2) {
        data <- sim_two_covs(ate, n, rct_prop, g_rct, bias, controls_only)
        W_node <- 2:3
        A_node <- 4
        Y_node <- 5
      } else if (num_covs == 4) {
        data <- sim_four_covs(ate, n, rct_prop, g_rct, bias, controls_only)
        W_node <- 2:5
        A_node <- 6
        Y_node <- 7
      }

      # fit
      res <- NULL
      if (method == "atmle") {
        res <- atmle(data = data,
                     S_node = S_node,
                     W_node = W_node,
                     A_node = A_node,
                     Y_node = Y_node,
                     controls_only = controls_only,
                     family = "gaussian",
                     atmle_pooled = TRUE,
                     theta_method = "glmnet",
                     Pi_method = nuisance_method,
                     g_method = nuisance_method,
                     theta_tilde_method = nuisance_method,
                     Q_method = nuisance_method,
                     bias_working_model = working_model,
                     pooled_working_model = "glmnet",
                     g_rct = g_rct,
                     verbose = FALSE)
      } else if (method == "atmle_tmle") {
        res <- atmle(data = data,
                     S_node = S_node,
                     W_node = W_node,
                     A_node = A_node,
                     Y_node = Y_node,
                     controls_only = controls_only,
                     family = "gaussian",
                     atmle_pooled = FALSE,
                     theta_method = "glmnet",
                     Pi_method = nuisance_method,
                     g_method = nuisance_method,
                     theta_tilde_method = nuisance_method,
                     Q_method = nuisance_method,
                     bias_working_model = working_model,
                     pooled_working_model = "glmnet",
                     g_rct = g_rct,
                     verbose = FALSE)
      } else if (method == "escvtmle") {
        covariates <- NULL
        if (num_covs == 2) {
          covariates <- c("W1", "W2")
        } else if (num_covs == 4) {
          covariates <- c("W1", "W2", "W3", "W4")
        }
        tmp <- ES.cvtmle(txinrwd = !controls_only,
                         data = data,
                         study = "S",
                         covariates = covariates,
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
                             controls_only = controls_only,
                             family = "gaussian",
                             atmle_pooled = TRUE,
                             theta_method = "glmnet",
                             Pi_method = nuisance_method,
                             g_method = nuisance_method,
                             theta_tilde_method = nuisance_method,
                             Q_method = nuisance_method,
                             bias_working_model = working_model,
                             pooled_working_model = "glmnet",
                             g_rct = g_rct,
                             verbose = FALSE)
      } else if (method == "rct_only") {
        res <- rct_only(data,
                        S_node = S_node,
                        W_node = W_node,
                        A_node = A_node,
                        Y_node = Y_node,
                        nuisance_method = nuisance_method,
                        family = "gaussian",
                        g_rct = g_rct,
                        verbose = FALSE)
      }

      if (res$lower <= ate & res$upper >= ate) {
        if (verbose) print("psi covered")
        psi_coverage[j] <- 1
      } else {
        if (verbose) print("psi not covered")
        psi_coverage[j] <- 0
      }

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
              escvtmle_prop_selected = escvtmle_prop_selected))
}
