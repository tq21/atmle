#' @title Adaptive-TMLE for Data Integration
#'
#' @import R6
atmle_ate_fusion <- R6Class(
  classname = "A-TMLE for RCT + RWD",
  inherit = atmle_ate,
  public = list(

    S_node = NULL,

    initialize = function(data,
                          S_node,
                          W_nodes,
                          A_node,
                          Y_node,
                          family,
                          n_folds = NULL,
                          seed = 123) {

      self$S_node <- S_node
      super$initialize(data = data,
                       W_nodes = W_nodes,
                       A_node = A_node,
                       Y_node = Y_node,
                       family = family,
                       n_folds = n_folds,
                       seed = seed)

    },

    target_Pi = function(target_gwt) {

      self$Pi_star <- self$Pi

      if (self$controls_only) {
        # only controls in external data
        if (target_gwt) {
          wt <- (1-self$A)/self$g0W
          HAW <- -self$tau_S$cate_W0
        } else {
          wt <- rep(1, length(self$A))
          HAW <- -(1-self$A)/self$g0W*self$tau_S$cate_W0
        }

        # logistic submodel, controls only
        epsilon <- as.numeric(coef(glm(S ~ -1+offset(qlogis(self$Pi$A0))+HAW,
                                       family = "quasibinomial", weights = wt))
        epsilon[is.na(epsilon)] <- 0

        # TMLE update
        if (target_gwt) {
          self$Pi_star$A0 <- .bound(plogis(qlogis(self$Pi$A0)+epsilon[1]*HAW), self$Pi_bounds)
          self$Pi_star$pred[self$A == 0] <- self$Pi_star$A0[self$A == 0]
        } else {
          self$Pi_star$A0 <- .bound(plogis(qlogis(self$Pi$A0)+epsilon[1]*(-1/self$g0W*self$tau_S$cate_W0)), self$Pi_bounds)
          self$Pi_star$pred[A == 0] <- self$Pi_star$A0[self$A == 0]
        }
      } else {
        # both treated and controls in external data
        if (target_gwt) {
          wt <- self$A/self$g1W+(1-self$A)/self$g0W
          H1_n <- self$tau_S$cate_W1*self$A
          H0_n <- self$tau_S$cate_W0*(1-self$A)
        } else {
          wt <- rep(1, length(self$A))
          H1_n <- self$A/self$g1W*self$tau_S$cate_W1
          H0_n <- (1-self$A)/self$g0W*self$tau_S$cate_W0
        }

        # logistic submodel, both treated and controls
        epsilon <- as.numeric(coef(glm(S ~ -1+offset(qlogis(self$Pi$pred))+H0_n+H1_n,
                                       family = "quasibinomial", weights = wt)))
        epsilon[is.na(epsilon)] <- 0

        # TMLE updates
        if (target_gwt) {
          self$Pi_star$pred <- .bound(plogis(qlogis(self$Pi$pred)+epsilon[1]*H0_n+epsilon[2]*H1_n), self$Pi_bounds)
          self$Pi_star$A0 <- .bound(plogis(qlogis(self$Pi$A0)+epsilon[1]*self$tau_S$cate_W0), self$Pi_bounds)
          self$Pi_star$A1 <- .bound(plogis(qlogis(self$Pi$A1)+epsilon[2]*self$tau_S$cate_W1), self$Pi_bounds)
        } else {
          self$Pi_star$pred <- .bound(plogis(qlogis(self$Pi$pred)+epsilon[1]*H0_n+epsilon[2]*H1_n), self$Pi_bounds)
          self$Pi_star$A0 <- .bound(plogis(qlogis(self$Pi$A0)+epsilon[1]*self$tau_S$cate_W0/self$g0W), self$Pi_bounds)
          self$Pi_star$A1 <- .bound(plogis(qlogis(self$Pi$A1)+epsilon[2]*self$tau_S$cate_W1/self$g1W), self$Pi_bounds)
        }
      }

      # update relevant parts of tau_S
      tau_S$pseudo_outcome <- ifelse(abs(self$S[self$Delta == 1]-self$Pi_star$pred[self$Delta == 1]) < 1e-10, 0,
                                     (self$Y[self$Delta == 1]-self$theta_WA[self$Delta == 1])/(self$S[self$Delta == 1]-self$Pi_star$pred[self$Delta == 1]))
      tau_S$pseudo_weights <- (self$S[self$Delta == 1]-self$Pi_star$pred[self$Delta == 1])^2*weights[self$Delta == 1]

      return(list(Pi = Pi_star,
                  tau_S = tau_S))
    },

    eval_Pi = function() {
      numer <- self$g$S1*self$Pi_bar
      denom <- numer+self$g$S0*(1-self$Pi_bar)
      self$Pi <- numer/denom

      return(invisible(self$Pi))
    },

    eval_g_bar = function() {
      self$g_bar <- self$g$S1*self$Pi_bar+self$g$S0*(1-self$Pi_bar)

      return(invisible(self$g_bar))
    },

    eval_Q_bar = function() {
      self$Q_bar$A <- self$Q$S1*self$Pi+self$Q$S0*(1-self$Pi)
      self$Q_bar$A1 <- self$Q$S1A1*self$Pi+self$Q$S0A1*(1-self$Pi)
      self$Q_bar$A0 <- self$Q$S1A0*self$Pi+self$Q$S0A0*(1-self$Pi)

      return(invisible(self$Q_bar))
    },

    eval_theta = function() {
      self$theta <- self$Q_bar$A1*self$g_bar+self$Q_bar$A0*(1-self$g_bar)

      return(invisible(self$theta))
    },

    run_init_est = function(Q_bound) {
      # estimate Q(S,W,A)=E(Y|S,W,A)
      self$Q_fit <- self$fit_regression(method = Q_method,
                                        folds = self$folds_obs,
                                        covariate_nodes = c(self$S_node,
                                                            self$W_nodes,
                                                            self$A_node),
                                        outcome_node = self$Y_node,
                                        bound = Q_bound)
      data_1WA <- self$data; data_1WA[[self$S_node]] <- 1
      data_1W1 <- self$data; data_1W1[[self$S_node]] <- 1; data_1W1[[self$A_node]] <- 1
      data_1W0 <- self$data; data_1W0[[self$S_node]] <- 1; data_1W0[[self$A_node]] <- 0
      data_0W1 <- self$data; data_0W1[[self$S_node]] <- 0; data_0W1[[self$A_node]] <- 1
      data_0W0 <- self$data; data_0W0[[self$S_node]] <- 0; data_0W0[[self$A_node]] <- 0
      Q1WA_task <- sl3_Task$new(data = data_1WA,
                                covariates = c(self$S_node, self$W_nodes, self$A_node),
                                folds = self$folds,
                                outcome = self$Y_node)
      Q1W1_task <- sl3_Task$new(data = data_1W1,
                                covariates = c(self$S_node, self$W_nodes, self$A_node),
                                folds = self$folds,
                                outcome = self$Y_node)
      Q1W0_task <- sl3_Task$new(data = data_1W0,
                                covariates = c(self$S_node, self$W_nodes, self$A_node),
                                folds = self$folds,
                                outcome = self$Y_node)
      Q0W1_task <- sl3_Task$new(data = data_0W1,
                                covariates = c(self$S_node, self$W_nodes, self$A_node),
                                folds = self$folds,
                                outcome = self$Y_node)
      Q0W0_task <- sl3_Task$new(data = data_0W0,
                                covariates = c(self$S_node, self$W_nodes, self$A_node),
                                folds = self$folds,
                                outcome = self$Y_node)
      self$Q$S1 <- self$Q_fit$predict(Q1WA_task)
      self$Q$S1A1 <- self$Q_fit$predict(Q1W1_task)
      self$Q$S1A0 <- self$Q_fit$predict(Q1W0_task)
      self$Q$S0A1 <- self$Q_fit$predict(Q0W1_task)
      self$Q$S0A0 <- self$Q_fit$predict(Q0W0_task)

      # estimate \bar{Pi}=P(S=1|W)
      self$fit_regression(method = Pi_bar_method,
                          folds = self$folds,
                          covariate_nodes = self$W_nodes,
                          outcome_node = self$S_node)

      # estimate g(1|S,W)=P(A=1|S,W)
      self$fit_regression(method = g_method,
                          folds = self$folds,
                          covariate_nodes = c(self$S_node, self$W_nodes),
                          outcome_node = self$A_node)
      data_A1 <- self$data; data_A1[[self$A_node]] <- 1
      data_A0 <- self$data; data_A0[[self$A_node]] <- 0
      g1W_task <- sl3_Task$new(data = data_A1,
                               covariates = c(self$S_node, self$W_nodes),
                               folds = self$folds,
                               outcome = self$A_node)
      g0W_task <- sl3_Task$new(data = data_A0,
                               covariates = c(self$S_node, self$W_nodes),
                               folds = self$folds,
                               outcome = self$A_node)
      self$g$S1 <- self$g_fit$predict(g1W_task)
      self$g$S0 <- self$g_fit$predict(g0W_task)

      # evaluate Pi, g_bar, Q_bar, theta
      self$eval_Pi(); self$eval_g_bar(); self$eval_Q_bar(); self$eval_theta()

    },

    fit_cate = function(W,
                        A,
                        Y,
                        g1W,
                        theta,
                        hal_args,
                        parallel) {
      # R-learner
      cate_fit <- list()
      cate_fit$pseudo_outcome <- ifelse(abs(A-g1W) < 1e-10, 0, (Y-theta)/(A-g1W))
      cate_fit$pseudo_weights <- (A-g1W)^2

      # make design matrix
      blist <- enumerate_basis(x = as.matrix(W),
                               max_degree = hal_args$max_degree,
                               smoothness_orders = hal_args$smoothness_orders,
                               num_knots = hal_args$num_knots)
      cate_fit$phi <- make_design_matrix(X = as.matrix(W), blist = blist)

      # fit HAL
      cate_fit$fit <- cv.glmnet(x = cate_fit$phi,
                                y = cate_fit$pseudo_outcome,
                                weights = cate_fit$pseudo_weights,
                                family = "gaussian",
                                alpha = 1,
                                foldid = self$foldsid,
                                parallel = parallel)

      return(invisible(cate_fit))

    },

    run = function(Q_method,
                   theta_method,
                   g_method,
                   Pi_method,
                   family = c("gaussian", "binomial"),
                   metalearner = Lrnr_cv_selector$new(),
                   A_cate_hal_args = list(max_degree = 3L,
                                          smoothness_orders = 1L,
                                          num_knots = 20L),
                   S_cate_hal_args = A_cate_hal_args,
                   target_method = "tmle",
                   browse = FALSE) {

      if (browse) browser()
      family <- match.arg(family)

      # initial estimation -----------------------------------------------------
      self$run_init_est(Q_bound = Q_bound)

      # obtain CATE working model ----------------------------------------------
      self$tau_A <- self$fit_cate(W = self$W,
                                  A = self$A,
                                  Y = self$Y,
                                  g1W = self$g_bar,
                                  theta = self$theta,
                                  hal_args = tau_A_hal_args,
                                  parallel = parallel)

      # obtain CARE working model ----------------------------------------------
      self$tau_S <- self$fit_cate(W = as.matrix(cbind(self$W, A=self$A)),
                                  A = self$S,
                                  Y = self$Y,
                                  g1W = self$Pi$A,
                                  theta = self$Q_bar$A,
                                  hal_args = S_cate_hal_args,
                                  parallel = parallel)

      # TODO: targeting

    }
  )
)
