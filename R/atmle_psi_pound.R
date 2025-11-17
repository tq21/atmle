atmle_psi_pound <- R6Class(
  classname = "A-TMLE for Psi Pound",
  inherit = atmle_ate,
  public = list(
    S_node = NULL,
    controls_only = NULL,
    target_wt = NULL,

    initialize = function(data,
                          S_node,
                          W_nodes,
                          A_node,
                          Y_node,
                          family,
                          n_folds = NULL,
                          seed = 123) {

      super$initialize(data = data,
                       W_nodes = c(W_nodes, A_node),
                       A_node = S_node,
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

    run = function(theta_method,
                   Pi_method,
                   family = c("gaussian", "binomial"),
                   S_cate_hal_args = list(max_degree = 3L,
                                          smoothness_orders = 1L,
                                          num_knots = 20L),
                   target_method = "tmle",
                   learn_theta_via_Q = TRUE,
                   g_bound = NULL,
                   verbose = FALSE,
                   browse = FALSE) {

      if (browse) browser()
      family <- match.arg(family)

      # initial estimation -----------------------------------------------------
      self$run_init_est(theta_method = theta_method,
                        g_method = Pi_method,
                        learn_theta_via_Q = learn_theta_via_Q,
                        g_bound = g_bound,
                        verbose = verbose,
                        browse = browse)

      # obtain CATE working model ----------------------------------------------
      self$tau_S$fit <- self$fit_cate(hal_args = cate_hal_args,
                                      parallel = parallel)

      # obtain tau_S(W,1) and tau_S(w,0)
      self$tau_S$phi_W1 <- make_design_matrix(X = as.matrix(cbind(self$W, A=1)),
                                              blist = self$tau_S$blist)
      self$tau_S$phi_W0 <- make_design_matrix(X = as.matrix(cbind(self$W, A=0)),
                                              blist = self$tau_S$blist)
      self$tau_S$pred_W1 <- as.numeric(cbind(1, self$tau_S$phi_W1) %*% self$tau_S$beta)
      self$tau_S$pred_W0 <- as.numeric(cbind(1, self$tau_S$phi_W0) %*% self$tau_S$beta)


      self$psi_pound_obj$fit_Q(Q_method = Q_pound_method,
                               family = family)
      self$psi_pound_obj$fit_g(g_method = g_method)
      self$psi_pound_obj$theta_pred <-
        self$psi_pound_obj$Q$Q1W*self$psi_pound_obj$g1W+
        self$psi_pound_obj$Q$Q0W*(1-self$psi_pound_obj$g1W)
      self$psi_pound_obj$fit_cate(hal_args = S_cate_hal_args)
      self$target_Pi(self$psi_pound_obj)

    }
  )
)
