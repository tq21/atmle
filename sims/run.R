library(purrr)
library(data.table)
library(sl3)
library(DataFusionDGPs)
library(doMC)
library(devtools)
load_all()
#registerDoMC(cores = 5)
`%+%` <- function(a, b) paste0(a, b)
set.seed(123)

run <- function(n,
                sim_data,
                sim_data_args,
                family,
                seed,
                truth,
                truth_2,
                browse = FALSE) {

  B <- 100 # number of Monte-Carlo runs
  SL_lib <- list(Lrnr_glm$new(),
                 Lrnr_xgboost$new(),
                 Lrnr_ranger$new(),
                 Lrnr_dbarts$new())
  lrnr_glm <- list(learners = list(Lrnr_cv$new(Lrnr_glm$new())))
  SL_con <- list(learners = SL_lib,
                 metalearner = Lrnr_cv_selector$new(loss_squared_error))
  SL_bin <- list(learners = SL_lib,
                 metalearner = Lrnr_cv_selector$new(loss_loglik_binomial))

  res_df <- map_dfr(n, function(.n) {

    atmle_psi <- atmle_lower <- atmle_upper <- rep(NA, B)
    atmle_psi_2 <- atmle_lower_2 <- atmle_upper_2 <- rep(NA, B)
    tmle_psi <- tmle_lower <- tmle_upper <- rep(NA, B)
    tmle_psi_2 <- tmle_lower_2 <- tmle_upper_2 <- rep(NA, B)

    walk(seq(B), function(.b) {

      if (browse) browser()

      cat("n: " %+% .n %+% ", b: " %+% .b %+% "...\n")
      set.seed(seed+.b)
      sim_data_args <- c(sim_data_args, list(n = .n))
      data <- do.call(sim_data, sim_data_args)

      # common arguments
      args <- list(data = data,
                   S_node = "S",
                   W_nodes = names(data)[grep("^W", names(data))],
                   A_node = "A",
                   Y_node = "Y",
                   family = family,
                   n_folds = 10)

      # A-TMLE -----------------------------------------------------------------
      tryCatch({
        atmle_obj <- do.call(atmle_ate_fusion$new, args)
        atmle_obj$run(Q_method = SL_con,
                      Pi_bar_method = lrnr_glm,
                      g_method = SL_bin,
                      target_method = "tmle",
                      max_iter = 50,
                      n_lambda = 1,
                      verbose = FALSE)
        atmle_psi[.b] <<- atmle_obj$results$psi[1]
        atmle_lower[.b] <<- atmle_obj$results$lower[1]
        atmle_upper[.b] <<- atmle_obj$results$upper[1]
        atmle_psi_2[.b] <<- atmle_obj$results$psi[2]
        atmle_lower_2[.b] <<- atmle_obj$results$lower[2]
        atmle_upper_2[.b] <<- atmle_obj$results$upper[2]
      }, error = function(e) {
        atmle_psi[.b] <<- atmle_lower[.b] <<- atmle_upper[.b] <<- NA
        atmle_psi_2[.b] <<- atmle_lower_2[.b] <<- atmle_upper_2[.b] <<- NA
      })
      cur_bias_atmle <- abs(mean(atmle_psi-truth, na.rm = TRUE))
      cur_se_atmle <- sd(atmle_psi, na.rm = TRUE)
      cur_mse_atmle <- mean((atmle_psi-truth)^2, na.rm = TRUE)
      cur_cover_atmle <- mean(atmle_lower <= truth & atmle_upper >= truth, na.rm = TRUE)
      cur_oracle_cover_atmle <- mean(atmle_psi+qnorm(0.025)*cur_se_atmle <= truth & atmle_psi+qnorm(0.975)*cur_se_atmle >= truth, na.rm = TRUE)
      cur_bias_atmle_2 <- abs(mean(atmle_psi_2-truth_2, na.rm = TRUE))
      cur_se_atmle_2 <- sd(atmle_psi_2, na.rm = TRUE)
      cur_mse_atmle_2 <- mean((atmle_psi_2-truth_2)^2, na.rm = TRUE)
      cur_cover_atmle_2 <- mean(atmle_lower_2 <= truth_2 & atmle_upper_2 >= truth_2, na.rm = TRUE)
      cur_oracle_cover_atmle_2 <- mean(atmle_psi_2+qnorm(0.025)*cur_se_atmle_2 <= truth_2 & atmle_psi_2+qnorm(0.975)*cur_se_atmle_2 >= truth_2, na.rm = TRUE)

      # nonparametric TMLE -----------------------------------------------------
      tryCatch({
        tmle_obj <- do.call(np_tmle_R6$new, args)
        tmle_obj$Q <- list(Q1WA = atmle_obj$Q$S1,
                           Q1W1 = atmle_obj$Q$S1A1,
                           Q1W0 = atmle_obj$Q$S1A0)
        tmle_obj$Pi_bar <- atmle_obj$Pi_bar
        tmle_obj$g11W <- atmle_obj$g$S1
        tmle_obj$run(Q_method = NULL,
                     Pi_bar_method = NULL,
                     g_method = NULL)
        tmle_psi[.b] <<- tmle_obj$results$psi[1]
        tmle_lower[.b] <<- tmle_obj$results$lower[1]
        tmle_upper[.b] <<- tmle_obj$results$upper[1]
        tmle_psi_2[.b] <<- tmle_obj$results$psi[2]
        tmle_lower_2[.b] <<- tmle_obj$results$lower[2]
        tmle_upper_2[.b] <<- tmle_obj$results$upper[2]
      }, error = function(e) {
        tmle_psi[.b] <<- tmle_lower[.b] <<- tmle_upper[.b] <<- NA
        tmle_psi_2[.b] <<- tmle_lower_2[.b] <<- tmle_upper_2[.b] <<- NA
      })
      cur_bias_tmle <- abs(mean(tmle_psi-truth, na.rm = TRUE))
      cur_se_tmle <- sd(tmle_psi, na.rm = TRUE)
      cur_mse_tmle <- mean((tmle_psi-truth)^2, na.rm = TRUE)
      cur_cover_tmle <- mean(tmle_lower <= truth & tmle_upper >= truth, na.rm = TRUE)
      cur_oracle_cover_tmle <- mean(tmle_psi+qnorm(0.025)*cur_se_tmle <= truth & tmle_psi+qnorm(0.975)*cur_se_tmle >= truth, na.rm = TRUE)
      cur_bias_tmle_2 <- abs(mean(tmle_psi_2-truth_2, na.rm = TRUE))
      cur_se_tmle_2 <- sd(tmle_psi_2, na.rm = TRUE)
      cur_mse_tmle_2 <- mean((tmle_psi_2-truth_2)^2, na.rm = TRUE)
      cur_cover_tmle_2 <- mean(tmle_lower_2 <= truth_2 & tmle_upper_2 >= truth_2, na.rm = TRUE)
      cur_oracle_cover_tmle_2 <- mean(tmle_psi_2+qnorm(0.025)*cur_se_tmle_2 <= truth_2 & tmle_psi_2+qnorm(0.975)*cur_se_tmle_2 >= truth_2, na.rm = TRUE)

      # print running results
      cat("Param: Avg. over pooled W\n")
      cat("A-TMLE: bias: " %+% format(round(cur_bias_atmle, 5), nsmall=5) %+%
            ", se: " %+% format(round(cur_se_atmle, 5), nsmall=5) %+%
            ", mse: " %+% format(round(cur_mse_atmle, 5), nsmall=5) %+%
            ", coverage: " %+% format(round(cur_cover_atmle, 2), nsmall=2) %+%
            ", oracle coverage: " %+% format(round(cur_oracle_cover_atmle, 2), nsmall=2) %+% "\n")
      cat("TMLE:   bias: " %+% format(round(cur_bias_tmle, 5), nsmall=5) %+%
            ", se: " %+% format(round(cur_se_tmle, 5), nsmall=5) %+%
            ", mse: " %+% format(round(cur_mse_tmle, 5), nsmall=5) %+%
            ", coverage: " %+% format(round(cur_cover_tmle, 2), nsmall=2) %+%
            ", oracle coverage: " %+% format(round(cur_oracle_cover_tmle, 2), nsmall=2) %+% "\n")
      cat("Param: Avg. over RCT W\n")
      cat("A-TMLE: bias: " %+% format(round(cur_bias_atmle_2, 5), nsmall=5) %+%
            ", se: " %+% format(round(cur_se_atmle_2, 5), nsmall=5) %+%
            ", mse: " %+% format(round(cur_mse_atmle_2, 5), nsmall=5) %+%
            ", coverage: " %+% format(round(cur_cover_atmle_2, 2), nsmall=2) %+%
            ", oracle coverage: " %+% format(round(cur_oracle_cover_atmle_2, 2), nsmall=2) %+% "\n")
      cat("TMLE:   bias: " %+% format(round(cur_bias_tmle_2, 5), nsmall=5) %+%
            ", se: " %+% format(round(cur_se_tmle_2, 5), nsmall=5) %+%
            ", mse: " %+% format(round(cur_mse_tmle_2, 5), nsmall=5) %+%
            ", coverage: " %+% format(round(cur_cover_tmle_2, 2), nsmall=2) %+%
            ", oracle coverage: " %+% format(round(cur_oracle_cover_tmle_2, 2), nsmall=2) %+% "\n\n")
    })

    return(rbind(data.frame(n = .n, b = seq_len(B), est_name = "atmle",
                            param = c("Avg. over pooled", "Avg. over RCT"),
                            psi = c(atmle_psi, atmle_psi_2),
                            lower = c(atmle_lower, atmle_lower_2),
                            upper = c(atmle_upper, atmle_upper_2)),
                 data.frame(n = .n, b = seq_len(B), est_name = "tmle",
                            param = c("Avg. over pooled", "Avg. over RCT"),
                            psi = c(tmle_psi, tmle_psi_2),
                            lower = c(tmle_lower, tmle_lower_2),
                            upper = c(tmle_upper, tmle_upper_2))))
  })

  return(res_df)
}
