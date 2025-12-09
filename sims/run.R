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

      # A-TMLE
      tryCatch({
        atmle_obj <- do.call(atmle_ate_fusion$new, args)
        atmle_obj$run(Q_method = SL_con,
                      Pi_bar_method = lrnr_glm,
                      g_method = SL_bin,
                      target_method = "tmle",
                      max_iter = 50,
                      n_lambda = 1,
                      verbose = FALSE)
        atmle_psi[.b] <<- atmle_obj$results$psi
        atmle_lower[.b] <<- atmle_obj$results$lower
        atmle_upper[.b] <<- atmle_obj$results$upper
      }, error = function(e) {
        atmle_psi[.b] <<- atmle_lower[.b] <<- atmle_upper[.b] <<- NA
      })
      cur_bias_atmle <- abs(mean(atmle_psi-truth, na.rm = TRUE))
      cur_se_atmle <- sd(atmle_psi, na.rm = TRUE)
      cur_mse_atmle <- mean((atmle_psi-truth)^2, na.rm = TRUE)
      cur_cover_atmle <- mean(atmle_lower <= truth & atmle_upper >= truth, na.rm = TRUE)
      cur_oracle_cover_atmle <- mean(atmle_psi+qnorm(0.025)*cur_se_atmle <= truth & atmle_psi+qnorm(0.975)*cur_se_atmle >= truth, na.rm = TRUE)

      # print running results
      cat("A-TMLE: bias: " %+% format(round(cur_bias_atmle, 5), nsmall=5) %+%
            ", se: " %+% format(round(cur_se_atmle, 5), nsmall=5) %+%
            ", mse: " %+% format(round(cur_mse_atmle, 5), nsmall=5) %+%
            ", coverage: " %+% format(round(cur_cover_atmle, 2), nsmall=2) %+%
            ", oracle coverage: " %+% format(round(cur_oracle_cover_atmle, 2), nsmall=2) %+% "\n")
    })

    return(rbind(data.frame(n = .n, b = seq_len(B),
                            est_name = "atmle",
                            psi = atmle_psi, lower = atmle_lower, upper = atmle_upper)))
  })

  return(res_df)
}
