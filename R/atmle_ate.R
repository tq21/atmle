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
    W = NULL,
    A = NULL,
    Y = NULL,
    n_folds = NULL,
    folds = NULL,
    Q_method = NULL,
    cross_fit_Q = NULL,
    Q_fit = NULL,
    g_fit = NULL,
    cate_fit = list(pseudo_outcome = NA,
                    pseudo_weights = NA,
                    fit = NA,
                    phi_W = NA),
    g1W = NULL,
    g0W = NULL,
    g_bounds = NULL,
    parallel = NULL,
    browse = NULL,
    Q = list(QAW = NA, Q1W = NA, Q0W = NA),
    theta_pred = NULL,
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
                          n_folds = NULL,
                          g_bounds = NULL,
                          seed = 123) {
      self$data <- data
      self$W_nodes <- W_nodes
      self$A_node <- A_node
      self$Y_node <- Y_node
      self$W <- data[, W_nodes, drop = FALSE]
      self$A <- data[[A_node]]
      self$Y <- data[[Y_node]]
      self$n_folds <- n_folds

      # make folds
      self$folds <- create_folds(n = nrow(data),
                                 n_folds = n_folds,
                                 seed = seed)

      # compute g-bound
      if (is.null(g_bounds)) {
        self$g_bounds <- c(5/sqrt(nrow(data))/log(nrow(data)), 1)
      }
    },

    fit_Q = function(Q_method,
                     family,
                     cross_fit,
                     discrete_SL,
                     hal_args,
                     parallel) {

      if (family == "gaussian") {
        loss_fn <- loss_squared_error
        outcome_type <- "continuous"
      } else if (family == "binomial") {
        loss_fn <- loss_loglik_binomial
        outcome_type <- "binomial"
      }
      if (length(Q_method) == 1) {
        lrnr <- Q_method[[1]]
      } else {
        if (discrete_SL) {
          metalearner <- Lrnr_cv_selector$new(eval_function = loss_fn)
        } else {
          metalearner <- Lrnr_nnls$new()
        }
        lrnr_stack <- Stack$new(Q_method)
        lrnr <- Lrnr_sl$new(learners = lrnr_stack, metalearner = metalearner)
      }
      # TODO: use Lrnr_cv to generate cross fitted predictions
      Q_task <- sl3_Task$new(data = self$data,
                             covariates = c(self$W_nodes, self$A_node),
                             outcome = self$Y_node,
                             outcome_type = outcome_type)
      self$Q_fit <- lrnr$train(Q_task)

      # make counterfactual data and predictions
      data_A1 <- self$data; data_A1[[self$A_node]] <- 1
      data_A0 <- self$data; data_A0[[self$A_node]] <- 0
      Q1_task <- sl3_Task$new(data = data_A1,
                              covariates = c(self$W_nodes, self$A_node),
                              outcome = self$Y_node,
                              outcome_type = outcome_type)
      Q0_task <- sl3_Task$new(data = data_A0,
                              covariates = c(self$W_nodes, self$A_node),
                              outcome = self$Y_node,
                              outcome_type = outcome_type)
      self$Q$QAW <- self$Q_fit$predict(Q_task)
      self$Q$Q0W <- self$Q_fit$predict(Q0_task)
      self$Q$Q1W <- self$Q_fit$predict(Q1_task)
    },

    fit_g = function(g_method,
                     cross_fit,
                     discrete_SL,
                     hal_args,
                     parallel) {

      if (length(g_method) == 1) {
        lrnr <- g_method[[1]]
      } else {
        if (discrete_SL) {
          metalearner <- Lrnr_cv_selector$new(eval_function = loss_loglik_binomial)
        } else {
          metalearner <- Lrnr_nnls$new()
        }
        lrnr_stack <- Stack$new(g_method)
        lrnr <- Lrnr_sl$new(learners = lrnr_stack, metalearner = metalearner)
      }
      # TODO: use Lrnr_cv to generate cross fitted predictions
      g_task <- sl3_Task$new(data = self$data,
                             covariates = self$W_nodes,
                             outcome = self$A_node)
      self$g_fit <- lrnr$train(g_task)
      self$g1W <- self$g_fit$predict(g_task)

    },

    fit_cate = function(hal_args,
                        parallel) {
      # R-learner
      self$cate_fit$pseudo_outcome <- ifelse(abs(self$A-self$g1W) < 1e-10, 0,
                                             (self$Y-self$theta_pred)/(self$A-self$g1W))
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
      R <- self$Y-self$theta_pred-(self$A-self$g1W)*tau
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
                                        theta = self$theta_pred,
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

    estimate = function(Q_method,
                        g_method,
                        family,
                        n_lambda = 1,
                        cross_fit_Q = TRUE,
                        discrete_SL_Q = TRUE,
                        cross_fit_g = TRUE,
                        discrete_SL_g = TRUE,
                        cate_hal_args = list(max_degree = 3L,
                                             smoothness_orders = 1L,
                                             num_knots = 20L),
                        target_method = "tmle",
                        hal_args = NULL,
                        parallel = FALSE,
                        browse = FALSE) {
      if (browse) browser()
      self$fit_Q(Q_method = Q_method,
                 family = family,
                 cross_fit = cross_fit_Q,
                 discrete_SL = discrete_SL_Q,
                 hal_args = hal_args,
                 parallel = parallel)
      self$fit_g(g_method = g_method,
                 cross_fit = cross_fit_g,
                 discrete_SL = discrete_SL_g,
                 hal_args = hal_args,
                 parallel = parallel)
      self$theta_pred <- self$Q$Q1W*self$g1W+self$Q$Q0W*(1-self$g1W)
      self$fit_cate(hal_args = cate_hal_args,
                    parallel = parallel)
      self$target(n_lambda = n_lambda,
                  method = target_method)
      self$inference()
    }
  )
)
