#' @title Adaptive-TMLE for Data Integration
#'
#' @import R6
atmle_ate_fusion <- R6Class(
  classname = "A-TMLE for RCT + RWD",
  inherit = tmle_R6,
  public = list(

    S_node = NULL,
    Pi_star = NULL,
    tau_S = NULL,
    tau_A = NULL,
    controls_only = NULL,
    target_gwt = NULL,

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

    target_Pi = function() {

      self$Pi_star <- self$Pi

      if (self$controls_only) {
        # only controls in external data
        if (self$target_gwt) {
          wt <- (1-self$A)/self$g_bar0
          H0W <- -self$tau_S$cate_W0
        } else {
          wt <- rep(1, length(self$A))
          H0W <- -(1-self$A)/self$g_bar0*self$tau_S$cate_W0
        }

        # logistic submodel
        epsilon <- as.numeric(coef(glm(S ~ -1+offset(qlogis(self$Pi$A))+H0W,
                                       family = "quasibinomial", weights = wt)))
        epsilon[is.na(epsilon)] <- 0

        # update
        if (self$target_gwt) {
          self$Pi_star$A0 <- plogis(qlogis(self$Pi$A0)+epsilon[1]*H0W)
          self$Pi_star$A[self$A == 0] <- self$Pi_star$A0[self$A == 0]
        } else {
          self$Pi_star$A0 <- plogis(qlogis(self$Pi$A0)+epsilon[1]*(-1/self$g_bar0*self$tau_S$cate_W0))
          self$Pi_star$A[A == 0] <- self$Pi_star$A0[self$A == 0]
        }
      } else {
        # both treated and controls in external data
        if (self$target_gwt) {
          wt <- self$A/self$g_bar+(1-self$A)/self$g_bar0
          H1W <- self$tau_S$cate_W1*self$A
          H0W <- self$tau_S$cate_W0*(1-self$A)
        } else {
          wt <- rep(1, length(self$A))
          H1W <- self$A/self$g_bar*self$tau_S$cate_W1
          H0W <- -(1-self$A)/self$g_bar0*self$tau_S$cate_W0
        }

        # logistic submodel
        epsilon <- as.numeric(coef(glm(S ~ -1+offset(qlogis(self$Pi$A))+H0W+H1W,
                                       family = "quasibinomial", weights = wt)))
        epsilon[is.na(epsilon)] <- 0

        # updates
        if (self$target_gwt) {
          self$Pi_star$A <- plogis(qlogis(self$Pi$A)+epsilon[1]*H0W+epsilon[2]*H1W)
          self$Pi_star$A0 <- plogis(qlogis(self$Pi$A0)+epsilon[1]*self$tau_S$cate_W0)
          self$Pi_star$A1 <- plogis(qlogis(self$Pi$A1)+epsilon[2]*self$tau_S$cate_W1)
        } else {
          self$Pi_star$A <- plogis(qlogis(self$Pi$A)+epsilon[1]*H0W+epsilon[2]*H1W)
          self$Pi_star$A0 <- plogis(qlogis(self$Pi$A0)+epsilon[1]*self$tau_S$cate_W0/self$g_bar0)
          self$Pi_star$A1 <- plogis(qlogis(self$Pi$A1)+epsilon[2]*self$tau_S$cate_W1/self$g_bar)
        }
      }

      # update relevant parts of tau_S
      tau_S$pseudo_outcome <- ifelse(abs(self$S[self$Delta == 1]-self$Pi_star$A[self$Delta == 1]) < 1e-10, 0,
                                     (self$Y[self$Delta == 1]-self$Q_bar$A[self$Delta == 1])/(self$S[self$Delta == 1]-self$Pi_star$A[self$Delta == 1]))
      tau_S$pseudo_weights <- (self$S[self$Delta == 1]-self$Pi_star$A[self$Delta == 1])^2*self$weights[self$Delta == 1]

    },

    target_tau = function(n_lambda,
                          cate_fit,
                          tau_A) {

      # target in a sequence of working models (or cv selected WM if n_lambda = 1)
      cv_lambda <- cate_fit$fit$lambda.min
      lambda_seq <- cate_fit$fit$lambda
      lambda_seq <- lambda_seq[lambda_seq <= cv_lambda]
      lambda_seq <- lambda_seq[1:min(n_lambda, length(lambda_seq))]

      # extract a sequence of working models indexed by lambda
      wm_seq <- map(lambda_seq, function(.lambda) {
        list(non_zero = which(as.numeric(coef(cate_fit$fit, s = .lambda))[-1] != 0),
             lambda = .lambda)
      })

      # use the CV selected fit as initial fit
      intercept <- as.numeric(coef(cate_fit$fit, s = cv_lambda))[1]
      beta <- as.numeric(coef(cate_fit$fit, s = cv_lambda))[-1]

      # perform targeting in each working model
      res_list <- map(seq_along(wm_seq), function(.j) {
        cur_wm <- wm_seq[[.j]]
        phi_select_j <- cate_fit$phi[, cur_wm$non_zero, drop = FALSE]
        phi_select_j <- cbind(1, phi_select_j)
        beta_j <- c(intercept, beta[cur_wm$non_zero])

        if (length(cur_wm$non_zero) > 0) {
          if (self$beta_target_method == "relaxed") {
            beta_star <- self$target_relaxed(pseudo_outcome = self$cate_fit$pseudo_outcome,
                                             pseudo_weights = self$cate_fit$pseudo_weights,
                                             phi_W = phi_select_j)
          } else if (self$beta_target_method == "tmle") {
            if (tau_A) {
              beta_star <- self$target_tau_A_tmle(phi_W = phi_select_j,
                                                  beta = beta_j)
            } else {
              phi_W1 <- cate_fit$phi_W1[, cur_wm$non_zero, drop = FALSE]
              phi_W1 <- cbind(1, phi_W1)
              phi_W0 <- cate_fit$phi_W0[, cur_wm$non_zero, drop = FALSE]
              phi_W0 <- cbind(1, phi_W0)
              beta_star <- self$target_tau_S_tmle(phi_WA = phi_select_j,
                                                  phi_W1 = phi_W1,
                                                  phi_W0 = phi_W0,
                                                  beta = beta_j)
            }
          }
        } else {
          beta_star <- mean(cate_fit$pseudo_outcome)
        }
        beta_star[is.na(beta_star)] <- 0
        cate_pred <- as.numeric(phi_select_j %*% beta_star)

        if (tau_A) {
          return(list(beta_star = beta_star,
                      phi_W = phi_select_j,
                      cate_pred = cate_pred))
        } else {
          cate_pred_A1 <- as.numeric(phi_W1 %*% beta_star)
          cate_pred_A0 <- as.numeric(phi_W0 %*% beta_star)
          return(list(beta_star = beta_star,
                      phi_WA = phi_select_j,
                      phi_W1 = phi_W1,
                      phi_W0 = phi_W0,
                      cate_pred = list(A = cate_pred,
                                       A1 = cate_pred_A1,
                                       A0 = cate_pred_A0)))
        }
      })

      return(res_list)
    },

    #' TMLE targeting of beta for tau_A
    target_tau_A_tmle = function(phi_W,
                                 beta) {

      phi_W <- as.matrix(phi_W)
      IM <- t(phi_W) %*% diag(self$g_bar*self$g_bar0) %*% phi_W / nrow(phi_W)
      IM_inv <- mat_inverse(IM)
      clever_cov <- as.vector(IM_inv %*% colMeans(phi_W))
      H <- (self$A-self$g_bar)*as.vector(phi_W %*% clever_cov)
      tau <- as.numeric(phi_W %*% beta)
      R <- self$Y-self$theta-(self$A-self$g1W)*tau
      epsilon <- sum(H*R)/sum(H*H)
      beta <- beta+epsilon*clever_cov

      return(beta)

    },

    #' TMLE targeting of beta for tau_S
    target_tau_S_tmle = function(phi_WA,
                                 phi_W1,
                                 phi_W0,
                                 beta) {

      phi_WA <- as.matrix(phi_WA)
      phi_W1 <- as.matrix(phi_W1)
      phi_W0 <- as.matrix(phi_W0)
      IM <- t(phi_WA) %*% diag(Pi$A*(1-Pi$A)) %*% phi_WA / nrow(phi_WA)
      IM_inv <- mat_inverse(IM)
      if (self$controls_only) {
        clever_cov <- as.vector(IM_inv %*% colMeans((1-Pi$A0)*phi_W0))
      } else {
        clever_cov <- as.vector(IM_inv %*% colMeans((1-Pi$A0)*phi_W0-(1-Pi$A1)*phi_W1))
      }
      H <- (S-Pi$A)*as.vector(phi_WA %*% clever_cov)
      tau <- as.numeric(phi_WA %*% beta)
      R <- self$Y-self$Q_bar$A-(S-Pi$A)*tau
      epsilon <- sum(H*R)/sum(H*H)
      beta <- beta+epsilon*clever_cov

      return(beta)

    },

    target_tau_relaxed = function(pseudo_outcome,
                                  pseudo_weights,
                                  phi) {

      relax_fit <- glm(pseudo_outcome ~ .,
                       family = "gaussian",
                       data = data.frame(as.matrix(phi[, 2:ncol(phi), drop=FALSE])),
                       weights = pseudo_weights)
      beta <- as.numeric(coef(relax_fit))

      return(beta)
    },

    #' Iterative targeting of Pi and beta_S
    target = function(max_iter,
                      verbose) {

      cur_iter <- 1
      PnEIC <- Inf
      sn <- 0
      while (cur_iter <= max_iter & abs(PnEIC) > sn) {
        # target Pi
        self$target_Pi()

        # TODO: do we need this?
        # evaluate EIC of psi pound
        # self$eic_psi_pound <- eic_psi_pound_wm(S = self$S,
        #                                        Y = self$Y,
        #                                        A = self$A,
        #                                        g1W = g_bar,
        #                                        theta_WA = self$Q_bar$A,
        #                                        Pi = self$Pi_star,
        #                                        tau_S = self$tau_S,
        #                                        weights = self$weights,
        #                                        controls_only = self$controls_only)
        # PnEIC <- mean(self$eic_psi_pound)
        # sn <- 0.001*sqrt(var(self$eic_psi_pound))/(sqrt(length(self$Y))*log(length(self$Y)))
        # TODO: should we allow break here?
        # if (abs(PnEIC) <= sn) {
        #   break
        # }

        # target beta_S
        self$target_tau(n_lambda = self$n_lambda,
                        cate_fit = self$tau_S,
                        tau_A = FALSE)
        # TODO: need to do this for every tau in the wm seq

        # re-evaluate EIC of psi pound
        self$eic_psi_pound <- eic_psi_pound_wm(S = self$S,
                                               Y = self$Y,
                                               A = self$A,
                                               g1W = g_bar,
                                               theta_WA = self$Q_bar$A,
                                               Pi = self$Pi_star,
                                               tau_S = self$tau_S,
                                               weights = self$weights,
                                               controls_only = self$controls_only)
        PnEIC <- mean(self$eic_psi_pound)
        sn <- sqrt(var(self$eic_psi_pound))/(sqrt(length(self$Y))*log(length(self$Y)))
        cur_iter <- cur_iter + 1
        if (verbose) print(round(PnEIC, 5))
      }

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
