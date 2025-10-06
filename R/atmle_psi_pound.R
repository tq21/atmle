atmle_psi_pound <- R6Class::R6Class(
  classname = "A-TMLE Psi Pound",
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
                          n_folds = NULL,
                          g_bound = NULL,
                          Pi_bound = NULL,
                          seed = 123) {
      super$initialize(data = data,
                       W_nodes = W_nodes,
                       A_node = A_node,
                       Y_node = Y_node,
                       n_folds = NULL,
                       g_bound = g_bound,
                       seed = seed)
    },

    target_Pi = function() {

      if (self$controls_only) {
        if (self$target_wt) {
          wt <- (1-self$A)/(1-self$g1W)
          HAW <- -tau_S$cate_W0
        } else {
          wt <- rep(1, length(A))
          HAW <- -(1-A)/(1-g1W)*tau_S$cate_W0
        }

        # logistic submodel, controls only
        epsilon <- coef(glm(S ~ -1 + offset(qlogis(Pi$A0)) + HAW,
                            family = "quasibinomial", weights = wt
        ))
        epsilon[is.na(epsilon)] <- 0

        # TMLE update
        if (target_gwt) {
          Pi_star$A0 <- .bound(plogis(qlogis(Pi$A0) + epsilon[1] * HAW), Pi_bounds)
          Pi_star$pred[A == 0] <- Pi_star$A0[A == 0]
        } else {
          Pi_star$A0 <- .bound(plogis(qlogis(Pi$A0) + epsilon[1] * (-1/(1-g1W)*tau_S$cate_W0)), Pi_bounds)
          Pi_star$pred[A == 0] <- Pi_star$A0[A == 0]
        }
      } else {
        if (target_gwt) {
          wt <- A/g1W+(1-A)/(1-g1W)
          H1_n <- tau_S$cate_W1*A
          H0_n <- tau_S$cate_W0*(1-A)
        } else {
          wt <- rep(1, length(A))
          H1_n <- A/g1W*tau_S$cate_W1
          H0_n <- (1-A)/(1-g1W)*tau_S$cate_W0
        }

        # logistic submodel, both treated and controls
        epsilon <- coef(glm(S ~ -1 + offset(qlogis(Pi$pred)) + H0_n + H1_n,
                            family = "quasibinomial", weights = wt
        ))
        epsilon[is.na(epsilon)] <- 0

        # TMLE updates
        if (target_gwt) {
          Pi_star$pred <- .bound(plogis(qlogis(Pi$pred) + epsilon[1] * H0_n + epsilon[2] * H1_n), Pi_bounds)
          Pi_star$A0 <- .bound(plogis(qlogis(Pi$A0) + epsilon[1] * tau_S$cate_W0), Pi_bounds)
          Pi_star$A1 <- .bound(plogis(qlogis(Pi$A1) + epsilon[2] * tau_S$cate_W1), Pi_bounds)
        } else {
          Pi_star$pred <- .bound(plogis(qlogis(Pi$pred) + epsilon[1] * H0_n + epsilon[2] * H1_n), Pi_bounds)
          Pi_star$A0 <- .bound(plogis(qlogis(Pi$A0) + epsilon[1] * tau_S$cate_W0/(1-g1W)), Pi_bounds)
          Pi_star$A1 <- .bound(plogis(qlogis(Pi$A1) + epsilon[2] * tau_S$cate_W1/g1W), Pi_bounds)
        }
      }

      # update relevant parts of tau_S
      tau_S$pseudo_outcome <- ifelse(abs(S[delta == 1]-Pi_star$pred[delta == 1]) < 1e-10, 0, (Y[delta == 1]-theta_WA[delta == 1])/(S[delta == 1]-Pi_star$pred[delta == 1]))
      tau_S$pseudo_weights <- (S[delta == 1]-Pi_star$pred[delta == 1])^2*weights[delta == 1]

      return(list(Pi = Pi_star,
                  tau_S = tau_S))
    },

    estimate = function(Q_method = list(learner = list(Lrnr_glmnet$new(),
                                                       Lrnr_xgboost$new()),
                                        discrete_SL = TRUE),
                        theta_method = Q_method,
                        g_method = Q_method,
                        Pi_method = Q_method,
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

      # estimate psi tilde parameter
      self$psi_tilde_obj$estimate(Q_method = Q_tilde_method,
                                  g_method = g_method,
                                  family = family,
                                  cate_hal_args = A_cate_hal_args,
                                  target_method = target_method,
                                  browse = browse)

      # estimate psi pound parameter
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
