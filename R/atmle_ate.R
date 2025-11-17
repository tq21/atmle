#' @title Adaptive-TMLE for ATE
#'
#' @import R6
atmle_ate <- R6Class(
  classname = "A-TMLE ATE",
  public = list(
    data = NULL,
    W_nodes = NULL,
    A_node = NULL,
    Y_node = NULL,
    Delta_node = NULL,
    family = NULL,
    W = NULL,
    A = NULL,
    Y = NULL,
    Delta = NULL,
    n_folds = NULL,
    folds = NULL,
    foldid = NULL,
    folds_obs = NULL,
    foldid_obs = NULL,
    theta_fit = NULL,
    Q_fit = NULL,
    g_fit = NULL,
    H = list(HAW = NA, H1W = NA, H0W = NA),
    cate_fit = list(pseudo_outcome = NA,
                    pseudo_weights = NA,
                    fit = NA,
                    phi_W = NA),
    g1W = NULL,
    g0W = NULL,
    g_bound = NULL,
    parallel = NULL,
    browse = NULL,
    Q = list(QAW = NA, Q1W = NA, Q0W = NA),
    theta = NULL,
    wm_seq = NULL,
    eic = vector(mode = "list"),
    results = NULL,
    phi_WA_full = NULL,
    phi_WA_select = NULL,
    aug_fit = NULL,
    phi_W_aug = NULL,

    initialize = function(data,
                          W_nodes,
                          A_node,
                          Y_node,
                          family,
                          n_folds = NULL,
                          seed = 123) {
      self$data <- data
      self$W_nodes <- W_nodes
      self$A_node <- A_node
      self$Y_node <- Y_node
      self$W <- data[, W_nodes, drop = FALSE]
      self$A <- data[[A_node]]
      self$Y <- data[[Y_node]]
      self$family <- family
      self$n_folds <- n_folds

    },

    create_folds = function(n_eff,
                            n_folds,
                            strata_ids,
                            seed) {

      set.seed(seed)



      self$folds <- make_folds(n = nrow(self$data), V = n_folds,
                               strata_ids = strata_ids)
      self$foldid <- folds2foldvec(self$folds)

    },

    fit_regression = function(method,
                              folds,
                              covariate_nodes,
                              outcome_node,
                              weight_node = NULL,
                              subset = seq(nrow(self$data)),
                              bound = NULL) {

      if (length(method$learners) == 1) {
        lrnr <- method$learners[[1]]
        if (!is.null(bound)) {
          lrnr_bound <- Lrnr_bound$new(bound)
          lrnr <- Pipeline$new(lrnr, lrnr_bound)
        }
      } else {
        lrnr_stack <- Stack$new(method$learners)
        if (!is.null(bound)) {
          lrnr_bound <- Lrnr_bound$new(bound)
          lrnr_stack <- Pipeline$new(lrnr_stack, lrnr_bound)
        }
        lrnr <- make_learner(Pipeline, Lrnr_cv$new(lrnr_stack),
                             method$metalearner)
      }
      suppressWarnings({
        task <- sl3_Task$new(data = self$data[subset, , drop = FALSE],
                             covariates = covariate_nodes,
                             outcome = outcome_node,
                             weights = weight_node,
                             folds = folds)
      })
      suppressMessages(fit_obj <- lrnr$train(task))

      return(invisible(fit_obj))

    },

    fit_cate = function(hal_args,
                        parallel) {
      # R-learner
      self$cate_fit$pseudo_outcome <- ifelse(abs(self$A-self$g1W) < 1e-10, 0,
                                             (self$Y-self$theta)/(self$A-self$g1W))
      self$cate_fit$pseudo_weights <- (self$A-self$g1W)^2

      # make design matrix
      basis_list <- enumerate_basis(x = as.matrix(self$W),
                                    max_degree = hal_args$max_degree,
                                    smoothness_orders = hal_args$smoothness_orders,
                                    num_knots = hal_args$num_knots)
      self$cate_fit$phi_W_full <- make_design_matrix(X = as.matrix(self$W), blist = basis_list)

      # fit HAL
      self$cate_fit$fit <- cv.glmnet(x = self$cate_fit$phi_W_full,
                                     y = self$cate_fit$pseudo_outcome,
                                     weights = self$cate_fit$pseudo_weights,
                                     family = "gaussian",
                                     alpha = 1,
                                     foldid = folds2foldvec(self$folds),
                                     parallel = parallel)

    },

    augment_wm = function(target_method,
                          parallel) {

      self$g0W <- 1-self$g1W
      self$g1W <- .bound(self$g1W, self$g_bounds)
      self$g0W <- .bound(self$g0W, self$g_bounds)
      clever_cov <- (self$A/self$g1W-(1-self$A)/self$g0W)

      # regress clever covariate on bases of W to extract influential bases
      self$aug_fit <- cv.glmnet(x = self$cate_fit$phi_W_full,
                                y = clever_cov,
                                alpha = 1,
                                family = "gaussian",
                                foldid = folds2foldvec(self$folds),
                                parallel = parallel)

      # drop columns that are already in the working model
      non_zero_idx <- which(as.numeric(coef(self$aug_fit, s = "lambda.min"))[-1] != 0)
      union_idx <- union(non_zero_idx, self$wm_seq[[1]]$non_zero)
      self$phi_W_aug <- self$cate_fit$phi_W_full[, union_idx, drop = FALSE]

      # targeting in the augmented working model
      cv_lambda <- self$cate_fit$fit$lambda.min
      intercept <- as.numeric(coef(self$cate_fit$fit, s = cv_lambda))[1]
      beta <- as.numeric(coef(self$cate_fit$fit, s = cv_lambda))[-1]
      beta <- beta[union_idx]
      if (ncol(self$phi_W_aug) > 0) {
        if (target_method == "relaxed") {
          beta_star <- self$target_relaxed(pseudo_outcome = self$cate_fit$pseudo_outcome,
                                           pseudo_weights = self$cate_fit$pseudo_weights,
                                           phi_W = as.matrix(cbind(1, self$phi_W_aug)))
        } else if (target_method == "tmle") {
          beta_star <- self$target_tmle(beta = c(intercept, beta),
                                        phi_W = cbind(1, self$phi_W_aug))
        }
      } else {
        beta_star <- mean(self$cate_fit$pseudo_outcome)
      }
      beta_star[is.na(beta_star)] <- 0
      cate_pred <- as.numeric(as.matrix(cbind(1, self$phi_W_aug)) %*% beta_star)
      self$wm_seq[[length(self$wm_seq)+1]] <- list(non_zero = union_idx,
                                                   lambda = 0,
                                                   cate_pred = cate_pred)
    },

    target = function(n_lambda,
                      method) {

      # target in a sequence of working models
      cv_lambda <- self$cate_fit$fit$lambda.min
      lambda_seq <- self$cate_fit$fit$lambda
      lambda_seq <- lambda_seq[lambda_seq <= cv_lambda]
      lambda_seq <- lambda_seq[1:min(n_lambda, length(lambda_seq))]

      # extract a sequence of working models indexed by lambda
      self$wm_seq <- map(lambda_seq, function(.lambda) {
        list(non_zero = which(as.numeric(coef(self$cate_fit$fit, s = .lambda))[-1] != 0),
             lambda = .lambda)
      })

      # use the CV selected fit as initial fit
      intercept <- as.numeric(coef(self$cate_fit$fit, s = cv_lambda))[1]
      beta <- as.numeric(coef(self$cate_fit$fit, s = cv_lambda))[-1]

      # perform targeting in each working model
      res_list <- map(seq_along(self$wm_seq), function(.j) {
        cur_wm <- self$wm_seq[[.j]]
        phi_W_select_j <- self$cate_fit$phi_W_full[, cur_wm$non_zero, drop = FALSE]
        phi_W_select_j <- cbind(1, phi_W_select_j)
        beta_j <- c(intercept, beta[cur_wm$non_zero])

        if (length(cur_wm$non_zero) > 0) {
          if (method == "relaxed") {
            beta_star <- self$target_relaxed(pseudo_outcome = self$cate_fit$pseudo_outcome,
                                             pseudo_weights = self$cate_fit$pseudo_weights,
                                             phi_W = phi_W_select_j)
          } else if (method == "tmle") {
            beta_star <- self$target_tmle(beta = beta_j,
                                          phi_W = phi_W_select_j)
          }
        } else {
          beta_star <- mean(self$cate_fit$pseudo_outcome)
        }
        beta_star[is.na(beta_star)] <- 0
        self$wm_seq[[.j]]$cate_pred <- as.numeric(phi_W_select_j %*% beta_star)
      })
    },

    target_tmle = function(beta,
                           phi_W) {

      phi_W <- as.matrix(phi_W)
      IM <- t(phi_W) %*% diag(self$g1W*(1-self$g1W)) %*% phi_W / nrow(phi_W)
      IM_inv <- mat_inverse(IM)
      clever_cov <- as.vector(IM_inv %*% colMeans(phi_W))
      H <- (self$A-self$g1W)*as.vector(phi_W %*% clever_cov)
      tau <- as.numeric(phi_W %*% beta)
      R <- self$Y-self$theta-(self$A-self$g1W)*tau
      epsilon <- sum(H*R)/sum(H*H)
      beta <- beta+epsilon*clever_cov

      return(beta)
    },

    target_relaxed = function(pseudo_outcome,
                              pseudo_weights,
                              phi_W) {

      relax_fit <- glm(pseudo_outcome ~ .,
                       family = "gaussian",
                       data = data.frame(as.matrix(phi_W[, 2:ncol(phi_W), drop=FALSE])),
                       weights = pseudo_weights)
      beta <- as.numeric(coef(relax_fit))

      return(beta)
    },

    inference = function(alpha = 0.05,
                         small_diag = 1e-3,
                         fall_back_method = "svd_pseudo_inv") {

      self$results <- map_dfr(seq_along(self$wm_seq), function(.j) {
        cur_wm <- self$wm_seq[[.j]]
        phi_W_select_j <- self$cate_fit$phi_W_full[, cur_wm$non_zero, drop = FALSE]
        phi_W_select_j <- as.matrix(cbind(1, phi_W_select_j))
        self$eic[[.j]] <- eic_ate_atmle(Y = self$Y,
                                        A = self$A,
                                        g1W = self$g1W,
                                        theta = self$theta,
                                        phi_W = phi_W_select_j,
                                        cate_pred = cur_wm$cate_pred,
                                        small_diag = small_diag,
                                        fall_back_method = fall_back_method)
        psi <- mean(cur_wm$cate_pred)
        se <- sqrt(var(self$eic[[.j]], na.rm = TRUE)/nrow(self$data))
        lower <- psi+qnorm(alpha/2)*se
        upper <- psi+qnorm(1-alpha/2)*se
        return(data.frame(lambda = cur_wm$lambda,
                          psi = psi, lower = lower, upper = upper, se = se))
      })
    },

    run_init_est = function(theta_method,
                            g_method,
                            learn_theta_via_Q = TRUE,
                            g_bound = NULL,
                            verbose = FALSE,
                            browse = FALSE) {

      if (browse) browser()

      # compute Delta ----------------------------------------------------------
      self$Delta_node <- "__Delta"
      i <- 1
      while (self$Delta_node %in% colnames(self$data)) {
        self$Delta_node <- paste0("__Delta", i)
        i <- i + 1
      }
      self$Delta <- as.numeric(!is.na(self$Y))
      n_eff <- sum(self$Delta)

      # make folds -------------------------------------------------------------
      if (is.null(self$folds_obs) | is.null(self$folds)) {
        # rule of thumb from tmle R package
        if (is.null(self$n_folds)) {
          if (n_eff <= 30){
            self$n_folds <- n_eff
          } else if (n_eff <= 500) {
            self$n_folds <- 20
          } else if (n_eff <= 1000) {
            self$n_folds <- 10
          } else if (n_eff <= 10000){
            self$n_folds <- 5
          } else {
            self$n_folds <- 3 # at least 3 for cv.glmnet to work
          }
        }
        if (self$family == "binomial") {
          strata_ids <- paste0(self$Delta, "-", self$A, "-", self$Y)
          strata_ids_obs <- paste0(self$A[self$Delta == 1], "-", self$Y[self$Delta == 1])
        } else {
          strata_ids <- paste0(self$Delta, "-", self$A)
          strata_ids_obs <- paste0(self$A[self$Delta == 1])
        }
        self$folds <- make_folds(n = nrow(self$data), V = self$n_folds,
                                 strata_ids = strata_ids)
        self$foldid <- folds2foldvec(self$folds)
        self$folds_obs <- make_folds(n = sum(self$Delta), V = self$n_folds,
                                     strata_ids = strata_ids_obs)
        self$foldid_obs <- folds2foldvec(self$folds_obs)
      }

      # fit g(1|W)=P(A=1|W) ----------------------------------------------------
      if (is.null(g_bound)) {
        self$g_bound <- 5/sqrt(n_eff)/log(n_eff)
      } else {
        self$g_bound <- g_bound
      }
      if (is.null(self$g1W)) {
        if (verbose) cat("Fitting g(1|W)=P(A=1|W)...\n")
        self$g_fit <- self$fit_regression(method = g_method,
                                          folds = self$folds,
                                          covariate_nodes = self$W_nodes,
                                          outcome_node = self$A_node)
        self$g1W <- self$g_fit$predict()
      }
      self$g0W <- 1-self$g1W
      self$g1W <- .bound(self$g1W, c(self$g_bound, 1))
      self$g0W <- .bound(self$g0W, c(self$g_bound, 1))

      # compute clever covariate
      self$H$HAW <- self$A/self$g1W-(1-self$A)/self$g0W
      self$H$H1W <- 1/self$g1W
      self$H$H0W <- -1/self$g0W

      # fit theta(W)=E(Y|W) ----------------------------------------------------
      if (is.null(self$theta)) {
        if (verbose) cat("Fitting theta(W)=E(Y|W)...\n")
        if (learn_theta_via_Q) {
          # theta(W)=E(Y|A=1,W)P(A=1|W)+E(Y|A=0,W)P(A=0|W)
          self$Q_fit <- self$fit_regression(method = theta_method,
                                            folds = self$folds_obs,
                                            covariate_nodes = c(self$W_nodes,
                                                                self$A_node),
                                            outcome_node = self$Y_node)

          # make counterfactual data and predictions
          data_A0 <- self$data; data_A0[[self$A_node]] <- 0
          data_A1 <- self$data; data_A1[[self$A_node]] <- 1
          Q1W_task <- sl3_Task$new(data = data_A1,
                                   covariates = c(self$W_nodes, self$A_node),
                                   outcome = self$Y_node,
                                   folds = self$folds)
          Q0W_task <- sl3_Task$new(data = data_A0,
                                   covariates = c(self$W_nodes, self$A_node),
                                   outcome = self$Y_node,
                                   folds = self$folds)
          self$Q$QAW <- self$Q_fit$predict()
          self$Q$Q0W <- self$Q_fit$predict(Q0W_task)
          self$Q$Q1W <- self$Q_fit$predict(Q1W_task)
          self$theta <- self$Q$Q1W*self$g1W+self$Q$Q0W*self$g0W
        } else {
          self$theta_fit <- self$fit_regression(method = theta_method,
                                                folds = self$folds,
                                                covariate_nodes = self$W_nodes,
                                                outcome_node = self$Y_node)
          self$theta <- self$theta_fit$predict()
        }
      }

    },

    run = function(theta_method,
                   g_method,
                   family,
                   learn_theta_via_Q = TRUE,
                   n_lambda = 1,
                   cate_hal_args = list(max_degree = 3L,
                                        smoothness_orders = 1L,
                                        num_knots = 20L),
                   target_method = "tmle",
                   parallel = FALSE,
                   browse = FALSE) {

      # obtain initial estimators ----------------------------------------------
      self$run_init_est(theta_method = theta_method,
                        g_method = g_method,
                        learn_theta_via_Q = learn_theta_via_Q,
                        browse = browse)

      # obtain CATE working model ----------------------------------------------
      self$fit_cate(hal_args = cate_hal_args,
                    parallel = parallel)

      # perform targeting ------------------------------------------------------
      self$target(n_lambda = n_lambda,
                  method = target_method)

      # point estimate and inference -------------------------------------------
      self$inference()
    }
  )
)
