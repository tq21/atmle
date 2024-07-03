library(devtools)
library(purrr)
load_all()

`%+%` <- function(a, b) paste0(a, b)

make_data <- function(B, n_seq) {
  data_list <- map(n_seq, function(n) {
    map(1:B, function(i) {
      data <- sim_data(n = n)
      return(data)
    })
  })

  return(data_list)
}

get_truth <- function(t0) {
  data_A1 <- sim_data(100000, A_counter = 1)
  data_A0 <- sim_data(100000, A_counter = 0)
  truth <- mean(data_A1$T_tilde > t0)-mean(data_A0$T_tilde > t0)

  return(truth)
}

run_sim <- function(t0,
                    truth,
                    data_list,
                    g_method,
                    lambda_method,
                    working_model,
                    cross_fit_nuisance) {

  # results
  # R-learner
  all_r_learner <- vector("list", length(data_list))
  all_r_learner_cover <- vector("list", length(data_list))
  all_r_learner_lower <- vector("list", length(data_list))
  all_r_learner_upper <- vector("list", length(data_list))

  # TMLE for survival outcomes
  all_tmle <- vector("list", length(data_list))
  all_tmle_cover <- vector("list", length(data_list))
  all_tmle_lower <- vector("list", length(data_list))
  all_tmle_upper <- vector("list", length(data_list))

  # R-learner + TMLE lambda
  all_r_learner_tmle_lambda <- vector("list", length(data_list))
  all_r_learner_tmle_lambda_cover <- vector("list", length(data_list))
  all_r_learner_tmle_lambda_lower <- vector("list", length(data_list))
  all_r_learner_tmle_lambda_upper <- vector("list", length(data_list))

  # R-learner + TMLE projection
  all_r_learner_tmle_proj <- vector("list", length(data_list))
  all_r_learner_tmle_proj_cover <- vector("list", length(data_list))
  all_r_learner_tmle_proj_lower <- vector("list", length(data_list))
  all_r_learner_tmle_proj_upper <- vector("list", length(data_list))

  for (i in 1:length(data_list)) {

    cur_data_list <- data_list[[i]]

    print("sample size: " %+% nrow(cur_data_list[[1]]))

    # R-learner
    r_learner <- vector(length = length(cur_data_list))
    r_learner_cover <- vector(length = length(cur_data_list))
    r_learner_lower <- vector(length = length(cur_data_list))
    r_learner_upper <- vector(length = length(cur_data_list))

    # TMLE for survival outcomes
    tmle <- vector(length = length(cur_data_list))
    tmle_cover <- vector(length = length(cur_data_list))
    tmle_lower <- vector(length = length(cur_data_list))
    tmle_upper <- vector(length = length(cur_data_list))

    # R-learner + TMLE lambda
    r_learner_tmle_lambda <- vector(length = length(cur_data_list))
    r_learner_tmle_lambda_cover <- vector(length = length(cur_data_list))
    r_learner_tmle_lambda_lower <- vector(length = length(cur_data_list))
    r_learner_tmle_lambda_upper <- vector(length = length(cur_data_list))

    # R-learner + TMLE projection
    r_learner_tmle_proj <- vector(length = length(cur_data_list))
    r_learner_tmle_proj_cover <- vector(length = length(cur_data_list))
    r_learner_tmle_proj_lower <- vector(length = length(cur_data_list))
    r_learner_tmle_proj_upper <- vector(length = length(cur_data_list))

    for (j in 1:length(cur_data_list)) {

      data <- cur_data_list[[j]]

      res <- atmle_surv(data = data,
                        W = c("W1", "W2"),
                        A = "A",
                        T_tilde = "T_tilde",
                        Delta = "Delta",
                        tau = 5,
                        t0 = t0,
                        g_method = g_method,
                        lambda_method = lambda_method,
                        cate_working_model = working_model,
                        cross_fit_nuisance = cross_fit_nuisance)

      # R-learner
      r_learner[j] <- res$r_learner
      r_learner_cover[j] <- as.numeric(res$r_learner_lower <= truth & truth <= res$r_learner_upper)
      r_learner_lower[j] <- res$r_learner_lower
      r_learner_upper[j] <- res$r_learner_upper

      # TMLE for survival outcomes
      tmle[j] <- res$tmle
      tmle_cover[j] <- as.numeric(res$tmle_lower <= truth & truth <= res$tmle_upper)
      tmle_lower[j] <- res$tmle_lower
      tmle_upper[j] <- res$tmle_upper

      # R-learner + TMLE lambda
      r_learner_tmle_lambda[j] <- res$r_learner_tmle_lambda
      r_learner_tmle_lambda_cover[j] <- as.numeric(res$r_learner_tmle_lambda_lower <= truth & truth <= res$r_learner_tmle_lambda_upper)
      r_learner_tmle_lambda_lower[j] <- res$r_learner_tmle_lambda_lower
      r_learner_tmle_lambda_upper[j] <- res$r_learner_tmle_lambda_upper

      # R-learner + TMLE projection
      r_learner_tmle_proj[j] <- res$r_learner_tmle_proj
      r_learner_tmle_proj_cover[j] <- as.numeric(res$r_learner_tmle_proj_lower <= truth & truth <= res$r_learner_tmle_proj_upper)
      r_learner_tmle_proj_lower[j] <- res$r_learner_tmle_proj_lower
      r_learner_tmle_proj_upper[j] <- res$r_learner_tmle_proj_upper

      print("(" %+% j %+% "/" %+% length(cur_data_list) %+% ") " %+%
              "R-learner: " %+% (round(mean(r_learner_cover[1:j]), 2)) %+%
              ", TMLE: " %+% (round(mean(tmle_cover[1:j]), 2)) %+%
              ", R-learner + TMLE lambda: " %+% (round(mean(r_learner_tmle_lambda_cover[1:j]), 2)) %+%
              ", R-learner + TMLE projection: " %+% (round(mean(r_learner_tmle_proj_cover[1:j]), 2)))
    }

    # R-learner
    all_r_learner[[i]] <- r_learner
    all_r_learner_cover[[i]] <- r_learner_cover
    all_r_learner_lower[[i]] <- r_learner_lower
    all_r_learner_upper[[i]] <- r_learner_upper

    # TMLE for survival outcomes
    all_tmle[[i]] <- tmle
    all_tmle_cover[[i]] <- tmle_cover
    all_tmle_lower[[i]] <- tmle_lower
    all_tmle_upper[[i]] <- tmle_upper

    # R-learner + TMLE lambda
    all_r_learner_tmle_lambda[[i]] <- r_learner_tmle_lambda
    all_r_learner_tmle_lambda_cover[[i]] <- r_learner_tmle_lambda_cover
    all_r_learner_tmle_lambda_lower[[i]] <- r_learner_tmle_lambda_lower
    all_r_learner_tmle_lambda_upper[[i]] <- r_learner_tmle_lambda_upper

    # R-learner + TMLE projection
    all_r_learner_tmle_proj[[i]] <- r_learner_tmle_proj
    all_r_learner_tmle_proj_cover[[i]] <- r_learner_tmle_proj_cover
    all_r_learner_tmle_proj_lower[[i]] <- r_learner_tmle_proj_lower
    all_r_learner_tmle_proj_upper[[i]] <- r_learner_tmle_proj_upper
  }

  return(list(all_r_learner = all_r_learner,
              all_r_learner_cover = all_r_learner_cover,
              all_r_learner_lower = all_r_learner_lower,
              all_r_learner_upper = all_r_learner_upper,
              all_tmle = all_tmle,
              all_tmle_cover = all_tmle_cover,
              all_tmle_lower = all_tmle_lower,
              all_tmle_upper = all_tmle_upper,
              all_r_learner_tmle_lambda = all_r_learner_tmle_lambda,
              all_r_learner_tmle_lambda_cover = all_r_learner_tmle_lambda_cover,
              all_r_learner_tmle_lambda_lower = all_r_learner_tmle_lambda_lower,
              all_r_learner_tmle_lambda_upper = all_r_learner_tmle_lambda_upper,
              all_r_learner_tmle_proj = all_r_learner_tmle_proj,
              all_r_learner_tmle_proj_cover = all_r_learner_tmle_proj_cover,
              all_r_learner_tmle_proj_lower = all_r_learner_tmle_proj_lower,
              all_r_learner_tmle_proj_upper = all_r_learner_tmle_proj_upper))
}
