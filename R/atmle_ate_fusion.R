#' @title Adaptive-TMLE for Data Integration
#'
#' @import R6
#'
#' @export
atmle_ate_fusion <- R6Class(
  classname = "A-TMLE for RCT + RWD",
  inherit = tmle_R6,
  public = list(

    S_node = NULL,
    S = NULL,
    Pi_star = NULL,
    Pi_star_avg_over_S1 = NULL,
    tau_S = NULL,
    tau_A = NULL,
    controls_only = NULL,
    target_gwt = NULL,
    beta_target_method = NULL,
    tau_A_star = NULL,
    tau_S_star = NULL,
    tau_A_star_avg_over_S1 = NULL,
    tau_S_star_avg_over_S1 = NULL,
    Delta = NULL,
    weights = NULL,
    Pi_bar_fit = NULL,
    g = list(S1 = NULL, S0 = NULL),
    Q = list(S1 = NULL, S1A1 = NULL, S1A0 = NULL, S0A1 = NULL, S0A0 = NULL),
    Pi = list(A = NULL, A1 = NULL, A0 = NULL),
    Pi_bar = NULL,
    g_bar = NULL,
    g_bar0 = NULL,
    Q_bar = list(A = NULL, A1 = NULL, A0 = NULL),
    theta = NULL,
    results = NULL,

    initialize = function(data,
                          S_node,
                          W_nodes,
                          A_node,
                          Y_node,
                          family,
                          n_folds = NULL) {

      super$initialize(data = data,
                       W_nodes = W_nodes,
                       A_node = A_node,
                       Y_node = Y_node,
                       family = family,
                       n_folds = n_folds)
      self$S_node <- S_node
      self$S <- data[[S_node]]
      self$Delta <- rep(1, nrow(data))
      self$weights <- rep(1, nrow(data))

    },

    #' Pi(1|A,W)=P(S=1|A,W)
    eval_Pi = function() {
      numer_A1 <- as.numeric(self$g$S1*self$Pi_bar)
      denom_A1 <- as.numeric(numer_A1+self$g$S0*(1-self$Pi_bar))
      numer_A0 <- as.numeric((1-self$g$S1)*self$Pi_bar)
      denom_A0 <- as.numeric(numer_A0+(1-self$g$S0)*(1-self$Pi_bar))
      self$Pi$A1 <- numer_A1/denom_A1
      self$Pi$A0 <- numer_A0/denom_A0
      self$Pi$A <- self$A*self$Pi$A1+(1-self$A)*self$Pi$A0

      return(invisible(self$Pi))
    },

    eval_g_bar = function() {
      self$g_bar <- as.numeric(self$g$S1*self$Pi_bar+self$g$S0*(1-self$Pi_bar))
      self$g_bar0 <- 1-self$g_bar

      return(invisible(self$g_bar))
    },

    #' Q_bar(A,W)=E(Y|A,W)
    eval_Q_bar = function() {
      self$Q_bar$A <- as.numeric(self$Q$S1*self$Pi$A+self$Q$S0*(1-self$Pi$A))
      self$Q_bar$A1 <- as.numeric(self$Q$S1A1*self$Pi$A1+self$Q$S0A1*(1-self$Pi$A0))
      self$Q_bar$A0 <- as.numeric(self$Q$S1A0*self$Pi$A1+self$Q$S0A0*(1-self$Pi$A0))

      return(invisible(self$Q_bar))
    },

    eval_theta = function() {
      self$theta <- as.numeric(self$Q_bar$A1*self$g_bar+self$Q_bar$A0*(1-self$g_bar))

      return(invisible(self$theta))
    },

    run_init_est = function(Q_method,
                            Pi_bar_method,
                            g_method,
                            Q_bound) {
      # estimate Q(S,W,A)=E(Y|S,W,A)
      self$Q_fit <- self$fit_regression(method = Q_method,
                                        folds = self$folds_obs,
                                        covariate_nodes = c(self$S_node,
                                                            self$W_nodes,
                                                            self$A_node),
                                        outcome_node = self$Y_node,
                                        bound = Q_bound)
      data_1WA <- self$data; data_1WA[[self$S_node]] <- 1
      data_0WA <- self$data; data_0WA[[self$S_node]] <- 0
      data_1W1 <- self$data; data_1W1[[self$S_node]] <- 1; data_1W1[[self$A_node]] <- 1
      data_1W0 <- self$data; data_1W0[[self$S_node]] <- 1; data_1W0[[self$A_node]] <- 0
      data_0W1 <- self$data; data_0W1[[self$S_node]] <- 0; data_0W1[[self$A_node]] <- 1
      data_0W0 <- self$data; data_0W0[[self$S_node]] <- 0; data_0W0[[self$A_node]] <- 0
      Q1WA_task <- sl3_Task$new(data = data_1WA,
                                covariates = c(self$S_node, self$W_nodes, self$A_node),
                                folds = self$folds,
                                outcome = self$Y_node)
      Q0WA_task <- sl3_Task$new(data = data_0WA,
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
      self$Q$S0 <- self$Q_fit$predict(Q0WA_task)
      self$Q$S1A1 <- self$Q_fit$predict(Q1W1_task)
      self$Q$S1A0 <- self$Q_fit$predict(Q1W0_task)
      self$Q$S0A1 <- self$Q_fit$predict(Q0W1_task)
      self$Q$S0A0 <- self$Q_fit$predict(Q0W0_task)

      # estimate \bar{Pi}=P(S=1|W)
      self$Pi_bar_fit <- self$fit_regression(method = Pi_bar_method,
                                             folds = self$folds,
                                             covariate_nodes = self$W_nodes,
                                             outcome_node = self$S_node)
      self$Pi_bar <- self$Pi_bar_fit$predict()

      # estimate g(1|S,W)=P(A=1|S,W)
      self$g_fit <- self$fit_regression(method = g_method,
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
      cate_fit$blist <- enumerate_basis(x = as.matrix(W),
                                        max_degree = hal_args$max_degree,
                                        smoothness_orders = hal_args$smoothness_orders,
                                        num_knots = hal_args$num_knots)
      cate_fit$phi <- make_design_matrix(X = as.matrix(W), blist = cate_fit$blist)

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

    target_Pi = function(cate_fit,
                         avg_over_S1) {

      if (is.null(self$Pi_star)) {
        self$Pi_star <- self$Pi
      }

      if (self$controls_only) {
        # only controls in external data
        if (self$target_gwt) {
          wt <- (1-self$A)/self$g_bar0
          H0W <- -(1-self$A)*cate_fit$cate_W0
          if (avg_over_S1) {
            wt <- self$S/mean(self$S)*wt
          }
        } else {
          wt <- rep(1, length(self$A))
          H0W <- -(1-self$A)/self$g_bar0*cate_fit$cate_W0
          if (avg_over_S1) {
            H0W <- self$S/mean(self$S)*H0W
          }
        }

        # logistic submodel
        epsilon <- as.numeric(coef(glm(S ~ -1+offset(qlogis(self$Pi_star$A))+H0W,
                                       family = "quasibinomial", weights = wt)))
        epsilon[is.na(epsilon)] <- 0

        # update
        if (self$target_gwt) {
          self$Pi_star$A0 <- plogis(qlogis(self$Pi_star$A0)+epsilon[1]*H0W)
          self$Pi_star$A[self$A == 0] <- self$Pi_star$A0[self$A == 0]
        } else {
          if (avg_over_S1) {
            self$Pi_star$A0 <- plogis(qlogis(self$Pi_star$A0)+epsilon[1]*(-self$S/mean(self$S)/self$g_bar0*cate_fit$cate_W0))
          } else {
            self$Pi_star$A0 <- plogis(qlogis(self$Pi_star$A0)+epsilon[1]*(-1/self$g_bar0*cate_fit$cate_W0))
          }
          self$Pi_star$A[A == 0] <- self$Pi_star$A0[self$A == 0]
        }
      } else {
        # both treated and controls in external data
        if (self$target_gwt) {
          wt <- self$A/self$g_bar+(1-self$A)/self$g_bar0
          H1W <- cate_fit$cate_W1*self$A
          H0W <- cate_fit$cate_W0*(1-self$A)
          if (avg_over_S1) {
            wt <- self$S/mean(self$S)*wt
          }
        } else {
          wt <- rep(1, length(self$A))
          H1W <- self$A/self$g_bar*cate_fit$cate_W1
          H0W <- (1-self$A)/self$g_bar0*cate_fit$cate_W0
          if (avg_over_S1) {
            H1W <- self$S/mean(self$S)*H1W
            H0W <- self$S/mean(self$S)*H0W
          }
        }

        # logistic submodel
        epsilon <- as.numeric(coef(glm(self$S ~ -1+offset(qlogis(self$Pi_star$A))+H0W+H1W,
                                       family = "quasibinomial", weights = wt)))
        epsilon[is.na(epsilon)] <- 0

        # updates
        if (self$target_gwt) {
          self$Pi_star$A <- plogis(qlogis(self$Pi_star$A)+epsilon[1]*H0W+epsilon[2]*H1W)
          self$Pi_star$A0 <- plogis(qlogis(self$Pi_star$A0)+epsilon[1]*cate_fit$cate_W0)
          self$Pi_star$A1 <- plogis(qlogis(self$Pi_star$A1)+epsilon[2]*cate_fit$cate_W1)
        } else {
          self$Pi_star$A <- plogis(qlogis(self$Pi_star$A)+epsilon[1]*H0W+epsilon[2]*H1W)
          if (avg_over_S1) {
            self$Pi_star$A0 <- plogis(qlogis(self$Pi_star$A0)+epsilon[1]*self$S/mean(self$S)/self$g_bar0*cate_fit$cate_W0)
            self$Pi_star$A1 <- plogis(qlogis(self$Pi_star$A1)+epsilon[2]*self$S/mean(self$S)/self$g_bar*cate_fit$cate_W1)
          } else {
            self$Pi_star$A0 <- plogis(qlogis(self$Pi_star$A0)+epsilon[1]*cate_fit$cate_W0/self$g_bar0)
            self$Pi_star$A1 <- plogis(qlogis(self$Pi_star$A1)+epsilon[2]*cate_fit$cate_W1/self$g_bar)
          }
        }
      }

      # update relevant parts of tau_S
      pseudo_outcome <- ifelse(abs(self$S[self$Delta == 1]-self$Pi_star$A[self$Delta == 1]) < 1e-10, 0,
                               (self$Y[self$Delta == 1]-self$Q_bar$A[self$Delta == 1])/(self$S[self$Delta == 1]-self$Pi_star$A[self$Delta == 1]))
      pseudo_weights <- (self$S[self$Delta == 1]-self$Pi_star$A[self$Delta == 1])^2*self$weights[self$Delta == 1]

      return(list(pseudo_outcome = pseudo_outcome,
                  pseudo_weights = pseudo_weights))
    },

    target_tau = function(cate_fit,
                          tau_A,
                          avg_over_S1) {

      if (tau_A) {
        phi <- cate_fit$phi_W
      } else {
        phi <- cate_fit$phi_WA
      }

      # perform targeting in the provided working model
      if (self$beta_target_method == "relaxed") {
        obj <- self$target_tau_relaxed(pseudo_outcome = cate_fit$pseudo_outcome,
                                       pseudo_weights = cate_fit$pseudo_weights,
                                       phi = phi)
      } else if (self$beta_target_method == "tmle") {
        if (tau_A) {
          obj <- self$target_tau_A_tmle(phi_W = phi,
                                        beta = cate_fit$beta,
                                        avg_over_S1 = avg_over_S1)
        } else {
          obj <- self$target_tau_S_tmle(phi_WA = phi,
                                        phi_W1 = cate_fit$phi_W1,
                                        phi_W0 = cate_fit$phi_W0,
                                        beta = cate_fit$beta,
                                        avg_over_S1 = avg_over_S1)
        }
      }
      beta_star <- obj$beta
      eic <- obj$eic
      beta_star[is.na(beta_star)] <- 0

      return(list(idx = cate_fit$idx,
                  lambda = cate_fit$lambda,
                  beta = beta_star,
                  eic = eic,
                  cate_W = obj$cate_W))

    },

    #' TMLE targeting of beta for tau_A
    target_tau_A_tmle = function(phi_W,
                                 beta,
                                 avg_over_S1) {

      phi_W <- as.matrix(phi_W)
      IM <- t(phi_W) %*% diag(self$g_bar*self$g_bar0) %*% phi_W / nrow(phi_W)
      IM_inv <- mat_inverse(IM)
      if (avg_over_S1) {
        # parameter that avgs over S=1 covariate distribution
        clever_cov <- as.vector(IM_inv %*% colMeans((self$S/mean(self$S))*phi_W))
      } else {
        # parameter that avgs over pooled covariate distribution
        clever_cov <- as.vector(IM_inv %*% colMeans(phi_W))
      }
      H <- (self$A-self$g_bar)*as.vector(phi_W %*% clever_cov)
      tau <- as.numeric(phi_W %*% beta)
      R <- self$Y-self$theta-(self$A-self$g_bar)*tau
      epsilon <- sum(H*R)/sum(H*H)
      beta <- beta+epsilon*clever_cov

      # compute EIC
      tau_star <- as.numeric(phi_W %*% beta)
      if (avg_over_S1) {
        D_beta <- as.vector(phi_W%*%IM_inv%*%colMeans(self$S/mean(self$S)*phi_W)*(self$A-self$g_bar)*
                              (self$Y-self$theta-(self$A-self$g_bar)*tau_star))
        W_comp <- self$S/mean(self$S)*(tau_star-mean(self$S/mean(self$S)*tau_star))
      } else {
        D_beta <- as.vector(phi_W%*%IM_inv%*%colMeans(phi_W)*(self$A-self$g_bar)*
                              (self$Y-self$theta-(self$A-self$g_bar)*tau_star))
        W_comp <- tau_star-mean(tau_star)
      }
      eic <- D_beta+W_comp

      return(list(beta = beta,
                  eic = eic,
                  cate_W = tau_star))

    },

    #' TMLE targeting of beta for tau_S
    target_tau_S_tmle = function(phi_WA,
                                 phi_W1,
                                 phi_W0,
                                 beta,
                                 avg_over_S1) {

      phi_WA <- as.matrix(phi_WA)
      phi_W1 <- as.matrix(phi_W1)
      phi_W0 <- as.matrix(phi_W0)
      IM <- t(phi_WA) %*% diag(self$Pi_star$A*(1-self$Pi_star$A)) %*% phi_WA / nrow(phi_WA)
      IM_inv <- mat_inverse(IM)
      if (self$controls_only) {
        if (avg_over_S1) {
          clever_cov <- as.vector(IM_inv %*% colMeans(self$S/mean(self$S)*(1-self$Pi_star$A0)*phi_W0))
        } else {
          clever_cov <- as.vector(IM_inv %*% colMeans((1-self$Pi_star$A0)*phi_W0))
        }
      } else {
        if (avg_over_S1) {
          clever_cov <- as.vector(IM_inv %*% colMeans(self$S/mean(self$S)*((1-self$Pi_star$A0)*phi_W0-(1-self$Pi_star$A1)*phi_W1)))
        } else {
          clever_cov <- as.vector(IM_inv %*% colMeans((1-self$Pi_star$A0)*phi_W0-(1-self$Pi_star$A1)*phi_W1))
        }
      }
      H <- (self$S-self$Pi_star$A)*as.vector(phi_WA %*% clever_cov)
      tau <- as.numeric(phi_WA %*% beta)
      R <- self$Y-self$Q_bar$A-(self$S-self$Pi_star$A)*tau
      epsilon <- sum(H*R)/sum(H*H)
      beta <- beta+epsilon*clever_cov

      # compute EIC
      tau_S <- list(phi_WA = phi_WA,
                    phi_W1 = phi_W1,
                    phi_W0 = phi_W0,
                    cate_WA = as.numeric(phi_WA %*% beta),
                    cate_W1 = as.numeric(phi_W1 %*% beta),
                    cate_W0 = as.numeric(phi_W0 %*% beta))
      eic <- eic_psi_pound_wm(S = self$S,
                              Y = self$Y,
                              A = self$A,
                              g1W = self$g_bar,
                              theta_WA = self$Q_bar$A,
                              Pi = self$Pi_star,
                              tau_S = tau_S,
                              weights = self$weights,
                              controls_only = self$controls_only,
                              IM_inv = IM_inv,
                              avg_over_S1 = avg_over_S1)

      return(list(beta = beta,
                  eic = eic,
                  cate_W = tau_S$cate_WA))

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

    #' Target beta_A over a sequence of working models
    target_beta_A_seq = function(n_lambda,
                                 avg_over_S1) {
      lambda_seq <- self$tau_A$fit$lambda
      lambda_cv <- self$tau_A$fit$lambda.min
      lambda_seq <- lambda_seq[lambda_seq <= lambda_cv]
      lambda_seq <- lambda_seq[1:min(n_lambda, length(lambda_seq))]
      beta_cv <- as.numeric(coef(self$tau_A$fit, s = lambda_cv))
      non_zero_cv <- which(beta_cv != 0)
      cate_W_cv <- as.numeric(cbind(1, self$tau_A$phi)[, non_zero_cv, drop = FALSE] %*% beta_cv[non_zero_cv])

      res_list <- map(seq_along(lambda_seq), function(.j) {

        # extract info on current working model
        lambda <- lambda_seq[.j]
        non_zero <- which(as.numeric(coef(self$tau_A$fit, s = lambda)) != 0)
        cate_fit <- list(phi_W = cbind(1, self$tau_A$phi)[, non_zero, drop = FALSE],
                         beta = beta_cv[non_zero],
                         pseudo_outcome = self$tau_A$pseudo_outcome,
                         pseudo_weights = self$tau_A$pseudo_weights)
        cur_res <- self$target_tau(cate_fit = cate_fit,
                                   tau_A = TRUE,
                                   avg_over_S1 = avg_over_S1)
        cur_res$idx <- .j
        cur_res$lambda <- lambda

        return(cur_res)
      })

      return(res_list)
    },

    #' Iterative targeting of Pi and beta_S
    target_Pi_beta_S = function(n_lambda,
                                avg_over_S1,
                                max_iter,
                                verbose) {

      # target in a sequence of working models (or cv selected WM if n_lambda = 1)
      cv_lambda <- self$tau_S$fit$lambda.min
      lambda_seq <- self$tau_S$fit$lambda
      lambda_seq <- lambda_seq[lambda_seq <= cv_lambda]
      lambda_seq <- lambda_seq[1:min(n_lambda, length(lambda_seq))]
      beta_cv <- as.numeric(coef(self$tau_S$fit, s = cv_lambda))
      non_zero_cv <- which(beta_cv != 0)
      phi_W1_cv <- cbind(1, self$tau_S$phi_W1)[, non_zero_cv, drop = FALSE]
      phi_W0_cv <- cbind(1, self$tau_S$phi_W0)[, non_zero_cv, drop = FALSE]
      phi_WA_cv <- cbind(1, self$tau_S$phi)[, non_zero_cv, drop = FALSE]
      cate_W1_cv <- as.numeric(phi_W1_cv %*% beta_cv[non_zero_cv])
      cate_W0_cv <- as.numeric(phi_W0_cv %*% beta_cv[non_zero_cv])
      cate_WA_cv <- as.numeric(phi_WA_cv %*% beta_cv[non_zero_cv])

      # for each working model, target Pi and beta_S iteratively
      res_list <- map(seq_along(lambda_seq), function(.j) {

        # extract info on current working model
        lambda <- lambda_seq[.j]
        non_zero <- which(as.numeric(coef(self$tau_S$fit, s = lambda)) != 0)
        cate_fit <- list(idx = .j,
                         lambda = lambda,
                         phi_W1 = cbind(1, self$tau_S$phi_W1)[, non_zero, drop = FALSE],
                         phi_W0 = cbind(1, self$tau_S$phi_W0)[, non_zero, drop = FALSE],
                         phi_WA = cbind(1, self$tau_S$phi)[, non_zero, drop = FALSE],
                         beta = beta_cv[non_zero],
                         cate_W1 = cate_W1_cv,
                         cate_W0 = cate_W0_cv,
                         cate_WA = cate_WA_cv,
                         pseudo_outcome = self$tau_S$pseudo_outcome,
                         pseudo_weights = self$tau_S$pseudo_weights)

        # target Pi and beta_S iteratively
        cur_iter <- 1
        PnEIC <- Inf
        sn <- 0
        while (cur_iter <= max_iter & abs(PnEIC) > sn) {
          # target Pi
          pseudo_list <- self$target_Pi(cate_fit, avg_over_S1)
          cate_fit$pseudo_outcome <- pseudo_list$pseudo_outcome
          cate_fit$pseudo_weights <- pseudo_list$pseudo_weights

          # target beta_S
          obj <- self$target_tau(cate_fit = cate_fit,
                                 tau_A = FALSE,
                                 avg_over_S1 = avg_over_S1)
          cate_fit$beta <- obj$beta
          cate_fit$eic <- obj$eic
          cate_fit$cate_WA <- obj$cate_W
          cate_fit$cate_W1 <- as.numeric(cate_fit$phi_W1 %*% cate_fit$beta)
          cate_fit$cate_W0 <- as.numeric(cate_fit$phi_W0 %*% cate_fit$beta)

          PnEIC <- mean(obj$eic)
          sn <- 1e-4*sqrt(var(obj$eic))/(sqrt(length(self$Y))*log(length(self$Y)))
          cur_iter <- cur_iter + 1
          if (verbose) print(round(PnEIC, 10))
        }

        cate_fit$Pi_star <- self$Pi_star
        self$Pi_star <- NULL

        return(cate_fit)
      })

      return(res_list)

    },

    inference = function(alpha = 0.05) {

      self$results <- map_dfr(self$tau_A_star, function(.tau_A) {
        df_psi <- map_dfr(self$tau_S_star, function(.tau_S) {
          # point estimate
          psi_tilde <- mean(.tau_A$cate_W)
          if (self$controls_only) {
            psi_pound <- mean((1-.tau_S$Pi_star$A0)*.tau_S$cate_W0)
          } else {
            psi_pound <- mean((1-.tau_S$Pi_star$A0)*.tau_S$cate_W0-(1-.tau_S$Pi_star$A1)*.tau_S$cate_W1)
          }
          psi <- psi_tilde-psi_pound

          # inference
          eic <- .tau_A$eic-.tau_S$eic
          se <- sqrt(var(eic, na.rm = TRUE)/nrow(self$data))
          lower <- psi+qnorm(alpha/2)*se
          upper <- psi+qnorm(1-alpha/2)*se

          return(data.frame(param = "Avg. over pooled",
                            tau_A_idx = .tau_A$idx,
                            tau_S_idx = .tau_S$idx,
                            tau_A_lambda = .tau_A$lambda,
                            tau_S_lambda = .tau_S$lambda,
                            psi_tilde = psi_tilde,
                            psi_pound = psi_pound,
                            psi = psi,
                            se = se,
                            lower = lower,
                            upper = upper,
                            alpha = alpha))
        })

        df_psi_avg_over_S1 <- map_dfr(self$tau_S_star_avg_over_S1, function(.tau_S) {
          # point estimate
          psi_tilde <- mean(self$S/mean(self$S)*.tau_A$cate_W)
          if (self$controls_only) {
            psi_pound <- mean(self$S/mean(self$S)*(1-.tau_S$Pi_star$A0)*.tau_S$cate_W0)
          } else {
            psi_pound <- mean(self$S/mean(self$S)*((1-.tau_S$Pi_star$A0)*.tau_S$cate_W0-(1-.tau_S$Pi_star$A1)*.tau_S$cate_W1))
          }
          psi <- psi_tilde-psi_pound

          # inference
          eic <- .tau_A$eic-.tau_S$eic
          se <- sqrt(var(eic, na.rm = TRUE)/nrow(self$data))
          lower <- psi+qnorm(alpha/2)*se
          upper <- psi+qnorm(1-alpha/2)*se

          return(data.frame(param = "Avg. over S=1",
                            tau_A_idx = .tau_A$idx,
                            tau_S_idx = .tau_S$idx,
                            tau_A_lambda = .tau_A$lambda,
                            tau_S_lambda = .tau_S$lambda,
                            psi_tilde = psi_tilde,
                            psi_pound = psi_pound,
                            psi = psi,
                            se = se,
                            lower = lower,
                            upper = upper,
                            alpha = alpha))
        })

        df_psi <- rbind(df_psi, df_psi_avg_over_S1)

        return(df_psi)
      })

      return(invisible(self$results))

    },

    run = function(Q_method,
                   Pi_bar_method,
                   g_method,
                   A_cate_hal_args = list(max_degree = 3L,
                                          smoothness_orders = 1L,
                                          num_knots = 20L),
                   S_cate_hal_args = list(max_degree = 3L,
                                          smoothness_orders = 1L,
                                          num_knots = 20L),
                   target_method = "tmle",
                   target_gwt = TRUE,
                   Q_bound = NULL,
                   max_iter = 50,
                   n_lambda = 10,
                   parallel = FALSE,
                   verbose = TRUE,
                   browse = FALSE) {

      if (browse) browser()

      self$controls_only <- all(self$A[self$S == 0] == 0)

      # initial estimation -----------------------------------------------------
      self$run_init_est(Q_method = Q_method,
                        Pi_bar_method = Pi_bar_method,
                        g_method = g_method,
                        Q_bound = Q_bound)

      # obtain CATE working model ----------------------------------------------
      self$tau_A <- self$fit_cate(W = self$W,
                                  A = self$A,
                                  Y = self$Y,
                                  g1W = self$g_bar,
                                  theta = self$theta,
                                  hal_args = A_cate_hal_args,
                                  parallel = parallel)

      # obtain CARE working model ----------------------------------------------
      self$tau_S <- self$fit_cate(W = as.matrix(cbind(self$W, A=self$A)),
                                  A = self$S,
                                  Y = self$Y,
                                  g1W = self$Pi$A,
                                  theta = self$Q_bar$A,
                                  hal_args = S_cate_hal_args,
                                  parallel = parallel)
      self$tau_S$phi_W1 <- make_design_matrix(X = as.matrix(cbind(self$W, A=1)),
                                              blist = self$tau_S$blist)
      self$tau_S$phi_W0 <- make_design_matrix(X = as.matrix(cbind(self$W, A=0)),
                                              blist = self$tau_S$blist)

      # 1. parameter that averages of pooled-covariate distribution ------------
      # target beta_A for each working model
      self$beta_target_method <- target_method
      self$tau_A_star <- self$target_beta_A_seq(n_lambda = n_lambda,
                                                avg_over_S1 = FALSE)

      # iterative targeting of Pi and beta_S
      self$target_gwt <- target_gwt
      self$tau_S_star <- self$target_Pi_beta_S(n_lambda = n_lambda,
                                               avg_over_S1 = FALSE,
                                               max_iter = max_iter,
                                               verbose = verbose)
      Pi_star_tmp <- self$Pi_star
      self$Pi_star <- NULL
      browser()

      # 2. parameter that averages over S=1 covariate distribution -------------
      # target beta_A for each working model
      self$tau_A_star_avg_over_S1 <- self$target_beta_A_seq(n_lambda = n_lambda,
                                                            avg_over_S1 = TRUE)

      # iterative targeting of Pi and beta_S
      self$target_gwt <- target_gwt
      self$tau_S_star_avg_over_S1 <- self$target_Pi_beta_S(n_lambda = n_lambda,
                                                           avg_over_S1 = TRUE,
                                                           max_iter = max_iter,
                                                           verbose = verbose)
      self$Pi_star_avg_over_S1 <- self$Pi_star
      self$Pi_star <- Pi_star_tmp

      # point estimate and inference -------------------------------------------
      self$inference()

    }
  )
)
