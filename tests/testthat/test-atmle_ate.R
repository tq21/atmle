library(dgps)
library(R6)
library(sl3)
library(devtools)
load_all()

lrnr_bin <- list(learners = list(Lrnr_dbarts$new(),
                                 Lrnr_xgboost$new()),
                 metalearner = Lrnr_cv_selector$new(loss_loglik_binomial))
lrnr_con <- list(learners = list(Lrnr_dbarts$new(),
                                 Lrnr_xgboost$new()),
                 metalearner = Lrnr_cv_selector$new(loss_squared_error))
lrnr_glm <- list(learners = list(Lrnr_cv$new(Lrnr_glm$new())))

data <- sim_data("dgp_gaussian_1", n = 500)

obj <- atmle_ate$new(data = data,
                     W_nodes = c("W1", "W2", "W3"),
                     A_node = "A",
                     Y_node = "Y",
                     n_folds = 5,
                     family = "gaussian",
                     seed = 123)
obj$run(theta_method = lrnr_con,
        g_method = lrnr_bin,
        family = "gaussian",
        learn_theta_via_Q = TRUE,
        n_lambda = 1,
        cate_hal_args = list(max_degree = 3L,
                             smoothness_orders = 1L,
                             num_knots = 20L),
        target_method = "tmle",
        parallel = FALSE,
        browse = TRUE)

obj$results
