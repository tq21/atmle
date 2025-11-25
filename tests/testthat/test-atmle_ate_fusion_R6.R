library(dgps)
library(R6)
library(sl3)
library(devtools)
load_all()
source("tests/testthat/utils.R")
set.seed(123)
lrnr_bin <- list(learners = list(Lrnr_dbarts$new(),
                                 Lrnr_xgboost$new()),
                 metalearner = Lrnr_cv_selector$new(loss_loglik_binomial))
lrnr_con <- list(learners = list(Lrnr_dbarts$new(),
                                 Lrnr_xgboost$new()),
                 metalearner = Lrnr_cv_selector$new(loss_squared_error))
lrnr_glm <- list(learners = list(Lrnr_cv$new(Lrnr_glm$new())))

data <- sim_data(n=1000, controls_only=FALSE, family="gaussian", prop_miss=0)

obj <- atmle_ate_fusion$new(data = data,
                            S_node = "S",
                            W_nodes = c("W1", "W2"),
                            A_node = "A",
                            Y_node = "Y",
                            family = "gaussian",
                            n_folds = 5)
obj$run(Q_method = lrnr_con,
        Pi_bar_method = lrnr_bin,
        g_method = lrnr_glm,
        family = "gaussian",
        target_method = "tmle",
        max_iter = 50,
        n_lambda = 50,
        browse = FALSE)

