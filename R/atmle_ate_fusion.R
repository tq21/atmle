#' Adaptive TMLE for RCT + RWD Data Fusion
#'
#' R6 class implementing an adaptive TMLE procedure for combining randomized
#' trial and real-world data.
#'
#' @details
#' The class estimates two ATE parameters:
#'
#' \itemize{
#'   \item Average over the pooled covariate distribution.
#'   \item Average over the trial-only covariate distribution (\eqn{S=1}).
#' }
#'
#' The implementation uses a direct nuisance strategy, then targets HAL working
#' model coefficients for \eqn{\tau_A(W)} and \eqn{\tau_S(W, A)}.
#'
#' @section Initialization:
#' Create the object with:
#'
#' `atmle_ate_fusion$new(data, S_node, W_nodes, A_node, Y_node, family, n_folds)`
#'
#' where:
#'
#' \itemize{
#'   \item `data`: input data frame.
#'   \item `S_node`: trial participation indicator column.
#'   \item `W_nodes`: baseline covariate columns.
#'   \item `A_node`: treatment column.
#'   \item `Y_node`: outcome column.
#'   \item `family`: outcome family (`"gaussian"` or `"binomial"`).
#'   \item `n_folds`: number of cross-fitting folds.
#' }
#'
#' @section Main Method (`$run()`):
#' The primary workflow is `fit$run(...)`.
#'
#' Key arguments:
#'
#' \itemize{
#'   \item `g_bar_method`, `theta_method`, `Pi_method`, `Q_bar_method`:
#'   sl3 nuisance-learning specifications.
#'   \item `A_cate_args`, `S_cate_args`: HAL basis tuning lists for the
#'   CATE/CARE working models.
#'   \item `target_method`: targeting routine for working-model coefficients.
#'   \item `target_gwt`, `max_iter`, `n_lambda`: iterative targeting controls.
#'   \item `g_bar_bound`, `theta_bound`, `Pi_bound`, `Q_bar_bound`:
#'   nuisance bounding controls.
#' }
#'
#' `A_cate_args` and `S_cate_args` should include HAL basis settings such as
#' `max_degree`, `smoothness_orders`, and `num_knots`. For the CARE working
#' model (`S_cate_args`), `force_A = TRUE` forces the main-effect basis for
#' treatment `A` into the HAL fit.
#'
#' @section Outputs:
#' After `run()`, the object contains:
#'
#' \itemize{
#'   \item `results`: final inference table for both estimands.
#'   \item `tau_A`, `tau_S`, `tau_A_star`, `tau_S_star`: intermediate and
#'   targeted working-model objects.
#' }
#'
#' @return An R6 class generator object. Use `$new(...)` to instantiate.
#'
#' @importFrom glmnet cv.glmnet
#' @importFrom hal9001 enumerate_basis make_design_matrix
#' @importFrom Matrix colSums
#' @importFrom origami make_folds folds2foldvec fold_from_foldvec
#' @importFrom purrr map map_dfr
#' @importFrom R6 R6Class
#' @importFrom sl3 sl3_Task
#' @importFrom stats coef gaussian glm plogis qlogis qnorm var
#' @export
atmle_ate_fusion <- R6::R6Class(
  classname = "A-TMLE for RCT + RWD",
  public = list(

    data = NULL,
    W_nodes = NULL,
    A_node = NULL,
    Y_node = NULL,
    W = NULL,
    A = NULL,
    Y = NULL,
    family = NULL,
    n_folds = NULL,
    Q_fit = NULL,
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
    g_bar_fit = NULL,
    theta_fit = NULL,
    Pi_fit = NULL,
    Q_bar_fit = NULL,
    g = list(S1 = NULL, S0 = NULL),
    Q = list(S1 = NULL, S1A1 = NULL, S1A0 = NULL, S0A1 = NULL, S0A0 = NULL),
    Pi = list(A = NULL, A1 = NULL, A0 = NULL),
    Pi_bar = NULL,
    g_bar = NULL,
    g_bar0 = NULL,
    Q_bar = list(A = NULL, A1 = NULL, A0 = NULL),
    theta = NULL,
    results = NULL,
    folds = NULL,
    foldsid = NULL,
    folds_S1 = NULL,
    folds_S0 = NULL,
    g_S1_fit = NULL,
    g_S0_fit = NULL,

    initialize = function(data,
                          S_node,
                          W_nodes,
                          A_node,
                          Y_node,
                          family,
                          n_folds = NULL) {

      family <- match.arg(family, c("gaussian", "binomial"))

      self$data <- data
      self$W_nodes <- W_nodes
      self$A_node <- A_node
      self$Y_node <- Y_node
      self$W <- data[, W_nodes, drop = FALSE]
      self$A <- data[[A_node]]
      self$Y <- data[[Y_node]]
      self$family <- family
      self$n_folds <- n_folds
      self$S_node <- S_node
      self$S <- data[[S_node]]
      self$Delta <- rep(1, nrow(data))
      self$weights <- rep(1, nrow(data))
    },

    run_init_est = function(g_bar_method,
                            theta_method,
                            Pi_method,
                            Q_bar_method,
                            g_bar_bound = NULL,
                            theta_bound = NULL,
                            Pi_bound = NULL,
                            Q_bar_bound = NULL) {
      # direct nuisance estimation for psi_tilde and psi_pound:
      # g_bar(W)=P(A=1|W), theta(W)=E(Y|W), Pi(A,W)=P(S=1|A,W), Q_bar(A,W)=E(Y|A,W)
      data_A1 <- self$data; data_A1[[self$A_node]] <- 1
      data_A0 <- self$data; data_A0[[self$A_node]] <- 0

      # g_bar(W)=P(A=1|W)
      self$g_bar_fit <- fit_regression(
        data = self$data,
        method = g_bar_method,
        folds = self$folds,
        covariate_nodes = self$W_nodes,
        outcome_node = self$A_node,
        bound = g_bar_bound
      )
      g_bar_task <- sl3_Task$new(
        data = self$data,
        covariates = self$W_nodes,
        folds = self$folds,
        outcome = self$A_node
      )
      self$g_bar <- as.numeric(self$g_bar_fit$predict(g_bar_task))
      self$g_bar0 <- 1 - self$g_bar

      # theta(W)=E(Y|W)
      self$theta_fit <- fit_regression(
        data = self$data,
        method = theta_method,
        folds = self$folds,
        covariate_nodes = self$W_nodes,
        outcome_node = self$Y_node,
        bound = theta_bound
      )
      theta_task <- sl3_Task$new(
        data = self$data,
        covariates = self$W_nodes,
        folds = self$folds,
        outcome = self$Y_node
      )
      self$theta <- as.numeric(self$theta_fit$predict(theta_task))

      # Pi(A,W)=P(S=1|A,W)
      self$Pi_fit <- fit_regression(
        data = self$data,
        method = Pi_method,
        folds = self$folds,
        covariate_nodes = c(self$W_nodes, self$A_node),
        outcome_node = self$S_node,
        bound = Pi_bound
      )
      Pi_A_task <- sl3_Task$new(
        data = self$data,
        covariates = c(self$W_nodes, self$A_node),
        folds = self$folds,
        outcome = self$S_node
      )
      Pi_A1_task <- sl3_Task$new(
        data = data_A1,
        covariates = c(self$W_nodes, self$A_node),
        folds = self$folds,
        outcome = self$S_node
      )
      Pi_A0_task <- sl3_Task$new(
        data = data_A0,
        covariates = c(self$W_nodes, self$A_node),
        folds = self$folds,
        outcome = self$S_node
      )
      self$Pi$A <- .bound(as.numeric(self$Pi_fit$predict(Pi_A_task)), c(1e-3, 1-1e-3))
      self$Pi$A1 <- .bound(as.numeric(self$Pi_fit$predict(Pi_A1_task)), c(1e-3, 1-1e-3))
      self$Pi$A0 <- .bound(as.numeric(self$Pi_fit$predict(Pi_A0_task)), c(1e-3, 1-1e-3))

      # Pi_bar(W)=P(S=1|W)
      self$Pi_bar_fit <- fit_regression(
        data = self$data,
        method = Pi_method,
        folds = self$folds,
        covariate_nodes = self$W_nodes,
        outcome_node = self$S_node
      )
      Pi_bar_task <- sl3_Task$new(
        data = self$data,
        covariates = self$W_nodes,
        folds = self$folds,
        outcome = self$S_node
      )
      self$Pi_bar <- as.numeric(self$Pi_bar_fit$predict(Pi_bar_task))

      # Q_bar(A,W)=E(Y|A,W)
      self$Q_bar_fit <- fit_regression(data = self$data,
                                       method = Q_bar_method,
                                       folds = self$folds,
                                       covariate_nodes = c(self$W_nodes, self$A_node),
                                       outcome_node = self$Y_node,
                                       bound = Q_bar_bound)
      Q_bar_A_task <- sl3_Task$new(data = self$data,
                                   covariates = c(self$W_nodes, self$A_node),
                                   folds = self$folds,
                                   outcome = self$Y_node)
      Q_bar_A1_task <- sl3_Task$new(data = data_A1,
                                    covariates = c(self$W_nodes, self$A_node),
                                    folds = self$folds,
                                    outcome = self$Y_node)
      Q_bar_A0_task <- sl3_Task$new(data = data_A0,
                                    covariates = c(self$W_nodes, self$A_node),
                                    folds = self$folds,
                                    outcome = self$Y_node)
      self$Q_bar$A <- as.numeric(self$Q_bar_fit$predict(Q_bar_A_task))
      self$Q_bar$A1 <- as.numeric(self$Q_bar_fit$predict(Q_bar_A1_task))
      self$Q_bar$A0 <- as.numeric(self$Q_bar_fit$predict(Q_bar_A0_task))

    },

    filter_hal_basis = function(blist,
                                phi,
                                force_keep_idx = integer(0)) {
      n <- nrow(phi)
      min_support_prop <- 1/sqrt(n)
      keep <- as.vector(Matrix::colSums(phi != 0) / n >= min_support_prop)
      if (length(force_keep_idx) > 0L) {
        keep[force_keep_idx] <- TRUE
      }

      # glmnet requires at least two columns; if filtering is too aggressive,
      # keep the original basis rather than failing downstream.
      if (sum(keep) < 2L) {
        keep <- rep(TRUE, length(blist))
      }

      return(list(
        blist = blist[keep],
        phi = phi[, keep, drop = FALSE],
        keep = keep
      ))
    },

    find_hal_main_effect_basis = function(blist,
                                          col_idx,
                                          smoothness_orders) {
      order_idx <- if (length(smoothness_orders) == 1L) {
        as.integer(smoothness_orders)
      } else {
        as.integer(smoothness_orders[col_idx])
      }

      basis_idx <- which(vapply(
        blist,
        function(b) {
          length(b$cols) == 1L &&
            identical(as.integer(b$cols), as.integer(col_idx)) &&
            length(b$cutoffs) == 1L &&
            isTRUE(all.equal(as.numeric(b$cutoffs), 0)) &&
            length(b$orders) == 1L &&
            identical(as.integer(b$orders), order_idx)
        },
        logical(1)
      ))

      if (length(basis_idx) != 1L) {
        stop("Could not uniquely identify the main-effect HAL basis for A.")
      }

      return(basis_idx)
    },

    fit_cate = function(W,
                        A,
                        Y,
                        g1W,
                        theta,
                        cate_args,
                        parallel,
                        force_A = FALSE,
                        default_hal_args = list(max_degree = 3L,
                                                smoothness_orders = 1L,
                                                num_knots = 20L)) {
      # R-learner
      cate_fit <- list()
      denom_bound <- 5/sqrt(nrow(self$data))/log(nrow(self$data))
      denom <- A-g1W
      denom[denom >= 0] <- .bound(denom[denom >= 0], c(denom_bound, 1))
      denom[denom < 0] <- .bound(denom[denom < 0], c(-1, -denom_bound))
      cate_fit$pseudo_outcome <- (Y-theta)/denom
      cate_fit$pseudo_weights <- (A-g1W)^2

      if (is.null(cate_args)) {
        cate_args <- list()
      }
      hal_args <- utils::modifyList(default_hal_args, cate_args)

      cate_fit$blist <- enumerate_basis(x = as.matrix(W),
                                        max_degree = hal_args$max_degree,
                                        smoothness_orders = hal_args$smoothness_orders,
                                        num_knots = hal_args$num_knots)
      cate_fit$phi <- make_design_matrix(X = as.matrix(W), blist = cate_fit$blist)
      force_keep_idx <- integer(0)
      if (force_A) {
        force_keep_idx <- self$find_hal_main_effect_basis(
          blist = cate_fit$blist,
          col_idx = ncol(as.matrix(W)),
          smoothness_orders = hal_args$smoothness_orders
        )
      }
      filtered_basis <- self$filter_hal_basis(blist = cate_fit$blist,
                                              phi = cate_fit$phi,
                                              force_keep_idx = force_keep_idx)
      cate_fit$blist <- filtered_basis$blist
      cate_fit$phi <- filtered_basis$phi
      penalty_factor <- rep(1, ncol(cate_fit$phi))
      if (force_A) {
        force_basis_idx <- self$find_hal_main_effect_basis(
          blist = cate_fit$blist,
          col_idx = ncol(as.matrix(W)),
          smoothness_orders = hal_args$smoothness_orders
        )
        penalty_factor[force_basis_idx] <- 0
      }

      cate_fit$fit <- cv.glmnet(x = cate_fit$phi,
                                y = cate_fit$pseudo_outcome,
                                weights = cate_fit$pseudo_weights,
                                family = "gaussian",
                                alpha = 1,
                                foldid = self$foldsid,
                                penalty.factor = penalty_factor,
                                maxit = 1e+06,
                                parallel = parallel)
      cate_fit$fit_type <- "glmnet"

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
            wt <- self$Pi_bar/mean(self$S)*wt
          }
        } else {
          wt <- rep(1, length(self$A))
          H0W <- -(1-self$A)/self$g_bar0*cate_fit$cate_W0
          if (avg_over_S1) {
            H0W <- self$Pi_bar/mean(self$S)*H0W
          }
        }

        # logistic submodel
        epsilon <- as.numeric(coef(glm(self$S ~ -1+offset(qlogis(self$Pi_star$A))+H0W,
                                       family = "quasibinomial", weights = wt)))
        epsilon[is.na(epsilon)] <- 0
        if (length(epsilon) == 0) {
          epsilon <- 0
        }

        # update
        if (self$target_gwt) {
          self$Pi_star$A0 <- plogis(qlogis(self$Pi_star$A0)+epsilon[1]*H0W)
          self$Pi_star$A[self$A == 0] <- self$Pi_star$A0[self$A == 0]
        } else {
          if (avg_over_S1) {
            self$Pi_star$A0 <- plogis(qlogis(self$Pi_star$A0)+epsilon[1]*(-self$Pi_bar/mean(self$S)/self$g_bar0*cate_fit$cate_W0))
          } else {
            self$Pi_star$A0 <- plogis(qlogis(self$Pi_star$A0)+epsilon[1]*(-1/self$g_bar0*cate_fit$cate_W0))
          }
          self$Pi_star$A[self$A == 0] <- self$Pi_star$A0[self$A == 0]
        }
      } else {
        # both treated and controls in external data
        if (self$target_gwt) {
          wt <- self$A/self$g_bar+(1-self$A)/self$g_bar0
          H1W <- cate_fit$cate_W1*self$A
          H0W <- cate_fit$cate_W0*(1-self$A)
          if (avg_over_S1) {
            wt <- self$Pi_bar/mean(self$S)*wt
          }
        } else {
          wt <- rep(1, length(self$A))
          H1W <- self$A/self$g_bar*cate_fit$cate_W1
          H0W <- (1-self$A)/self$g_bar0*cate_fit$cate_W0
          if (avg_over_S1) {
            H1W <- self$Pi_bar/mean(self$S)*H1W
            H0W <- self$Pi_bar/mean(self$S)*H0W
          }
        }

        # logistic submodel
        epsilon <- as.numeric(coef(glm(self$S ~ -1+offset(qlogis(self$Pi_star$A))+H0W+H1W,
                                       family = "quasibinomial", weights = wt)))
        epsilon[is.na(epsilon)] <- 0
        if (length(epsilon) < 2) {
          epsilon <- c(epsilon, rep(0, 2-length(epsilon)))
        }

        # updates
        if (self$target_gwt) {
          self$Pi_star$A <- plogis(qlogis(self$Pi_star$A)+epsilon[1]*H0W+epsilon[2]*H1W)
          self$Pi_star$A0 <- plogis(qlogis(self$Pi_star$A0)+epsilon[1]*cate_fit$cate_W0)
          self$Pi_star$A1 <- plogis(qlogis(self$Pi_star$A1)+epsilon[2]*cate_fit$cate_W1)
        } else {
          self$Pi_star$A <- plogis(qlogis(self$Pi_star$A)+epsilon[1]*H0W+epsilon[2]*H1W)
          if (avg_over_S1) {
            self$Pi_star$A0 <- plogis(qlogis(self$Pi_star$A0)+epsilon[1]*self$Pi_bar/mean(self$S)/self$g_bar0*cate_fit$cate_W0)
            self$Pi_star$A1 <- plogis(qlogis(self$Pi_star$A1)+epsilon[2]*self$Pi_bar/mean(self$S)/self$g_bar*cate_fit$cate_W1)
          } else {
            self$Pi_star$A0 <- plogis(qlogis(self$Pi_star$A0)+epsilon[1]*cate_fit$cate_W0/self$g_bar0)
            self$Pi_star$A1 <- plogis(qlogis(self$Pi_star$A1)+epsilon[2]*cate_fit$cate_W1/self$g_bar)
          }
        }
      }

      pi_bounds <- c(1e-3, 1-1e-3)
      self$Pi_star$A <- .bound(self$Pi_star$A, pi_bounds)
      self$Pi_star$A0 <- .bound(self$Pi_star$A0, pi_bounds)
      self$Pi_star$A1 <- .bound(self$Pi_star$A1, pi_bounds)

      # update relevant parts of tau_S
      denom_bound <- 5/sqrt(nrow(self$data))/log(nrow(self$data))
      denom <- self$S-self$Pi_star$A
      denom[denom >= 0] <- .bound(denom[denom >= 0], c(denom_bound, 1))
      denom[denom < 0] <- .bound(denom[denom < 0], c(-1, -denom_bound))
      pseudo_outcome <- (self$Y-self$Q_bar$A)/denom
      pseudo_weights <- (self$S-self$Pi_star$A)^2

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
        if (tau_A) {
          obj <- self$target_tau_A_relaxed(phi_W = phi,
                                           pseudo_outcome = cate_fit$pseudo_outcome,
                                           pseudo_weights = cate_fit$pseudo_weights,
                                           avg_over_S1 = avg_over_S1)
        } else {
          obj <- self$target_tau_S_relaxed(phi_WA = phi,
                                           phi_W1 = cate_fit$phi_W1,
                                           phi_W0 = cate_fit$phi_W0,
                                           pseudo_outcome = cate_fit$pseudo_outcome,
                                           pseudo_weights = cate_fit$pseudo_weights,
                                           avg_over_S1 = avg_over_S1)
        }
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

    # TMLE targeting of beta for tau_A
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
      eic_obj <- private$get_tau_A_eic(phi_W = phi_W,
                                       cate_W = tau_star,
                                       avg_over_S1 = avg_over_S1,
                                       IM_inv = IM_inv)

      return(list(beta = beta,
                  eic = eic_obj$eic,
                  cate_W = tau_star))

    },

    # TMLE targeting of beta for tau_S
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
      eic_obj <- private$get_tau_S_eic(phi_WA = phi_WA,
                                       phi_W1 = phi_W1,
                                       phi_W0 = phi_W0,
                                       beta = beta,
                                       avg_over_S1 = avg_over_S1,
                                       IM_inv = IM_inv)

      return(list(beta = beta,
                  eic = eic_obj$eic,
                  cate_W = eic_obj$cate_WA))

    },

    # Relaxed targeting of beta for tau_A
    target_tau_A_relaxed = function(phi_W,
                                    pseudo_outcome,
                                    pseudo_weights,
                                    avg_over_S1) {

      phi_W <- as.matrix(phi_W)
      relax_fit <- glm.fit(x = phi_W,
                           y = pseudo_outcome,
                           weights = pseudo_weights,
                           family = gaussian())
      beta <- as.numeric(relax_fit$coefficients)
      beta[is.na(beta)] <- 0
      tau_star <- as.numeric(phi_W %*% beta)

      # compute EIC
      eic_obj <- private$get_tau_A_eic(phi_W = phi_W,
                                       cate_W = tau_star,
                                       avg_over_S1 = avg_over_S1)

      return(list(beta = beta,
                  eic = eic_obj$eic,
                  cate_W = tau_star))
    },

    # Relaxed targeting of beta for tau_S
    target_tau_S_relaxed = function(phi_WA,
                                    phi_W1,
                                    phi_W0,
                                    pseudo_outcome,
                                    pseudo_weights,
                                    avg_over_S1) {

      phi_WA <- as.matrix(phi_WA)
      phi_W1 <- as.matrix(phi_W1)
      phi_W0 <- as.matrix(phi_W0)
      relax_fit <- glm.fit(x = phi_WA,
                           y = pseudo_outcome,
                           weights = pseudo_weights,
                           family = gaussian())
      beta <- as.numeric(relax_fit$coefficients)
      beta[is.na(beta)] <- 0

      # compute EIC
      IM <- t(phi_WA) %*% diag(self$Pi_star$A*(1-self$Pi_star$A)) %*% phi_WA / nrow(phi_WA)
      IM_inv <- mat_inverse(IM)
      eic_obj <- private$get_tau_S_eic(phi_WA = phi_WA,
                                       phi_W1 = phi_W1,
                                       phi_W0 = phi_W0,
                                       beta = beta,
                                       avg_over_S1 = avg_over_S1,
                                       IM_inv = IM_inv)

      return(list(beta = beta,
                  eic = eic_obj$eic,
                  cate_W = eic_obj$cate_WA))
    },

    # Target beta_A over a sequence of working models
    target_beta_A_seq = function(n_lambda,
                                 avg_over_S1) {
      lambda_seq <- self$tau_A$fit$lambda
      lambda_cv <- self$tau_A$fit$lambda.min
      lambda_seq <- lambda_seq[lambda_seq <= lambda_cv]
      lambda_seq <- lambda_seq[1:min(n_lambda, length(lambda_seq))]
      beta_cv <- as.numeric(coef(self$tau_A$fit, s = lambda_cv))
      non_zero_cv <- which(beta_cv != 0)

      res_list <- lapply(seq_along(lambda_seq), function(.j) {
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

    # Iterative targeting of Pi and beta_S
    target_Pi_beta_S = function(n_lambda,
                                avg_over_S1,
                                max_iter,
                                verbose) {

      # target in a sequence of working models (or cv selected WM if n_lambda = 1)
      if (identical(self$tau_S$fit_type, "glm")) {
        lambda_seq <- NA_real_
        beta_cv <- as.numeric(self$tau_S$fit$coefficients)
        beta_cv[is.na(beta_cv)] <- 0
        non_zero_cv <- seq_along(beta_cv)
      } else {
        cv_lambda <- self$tau_S$fit$lambda.min
        lambda_seq <- self$tau_S$fit$lambda
        lambda_seq <- lambda_seq[lambda_seq <= cv_lambda]
        lambda_seq <- lambda_seq[1:min(n_lambda, length(lambda_seq))]
        beta_cv <- as.numeric(coef(self$tau_S$fit, s = cv_lambda))
        non_zero_cv <- which(beta_cv != 0)
      }
      phi_W1_cv <- cbind(1, self$tau_S$phi_W1)[, non_zero_cv, drop = FALSE]
      phi_W0_cv <- cbind(1, self$tau_S$phi_W0)[, non_zero_cv, drop = FALSE]
      phi_WA_cv <- cbind(1, self$tau_S$phi)[, non_zero_cv, drop = FALSE]
      cate_W1_cv <- as.numeric(phi_W1_cv %*% beta_cv[non_zero_cv])
      cate_W0_cv <- as.numeric(phi_W0_cv %*% beta_cv[non_zero_cv])
      cate_WA_cv <- as.numeric(phi_WA_cv %*% beta_cv[non_zero_cv])

      res_list <- lapply(seq_along(lambda_seq), function(.j) {

        # extract info on current working model
        lambda <- lambda_seq[.j]
        if (identical(self$tau_S$fit_type, "glm")) {
          non_zero <- seq_along(beta_cv)
        } else {
          non_zero <- which(as.numeric(coef(self$tau_S$fit, s = lambda)) != 0)
        }
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
                  if (!is.finite(PnEIC) || !is.finite(sn)) {
                    break
                  }
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

      df_psi <- purrr::map_dfr(self$tau_A_star, function(.tau_A) {
        purrr::map_dfr(self$tau_S_star, function(.tau_S) {
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
      })

      df_psi_avg_over_S1 <- purrr::map_dfr(self$tau_A_star_avg_over_S1, function(.tau_A) {
        purrr::map_dfr(self$tau_S_star_avg_over_S1, function(.tau_S) {
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
      })

      self$results <- rbind(df_psi, df_psi_avg_over_S1)

      return(invisible(self$results))

    },

    run = function(g_bar_method,
                   theta_method,
                   Pi_method,
                   Q_bar_method,
                   A_cate_args = list(max_degree = 3L,
                                      smoothness_orders = 1L,
                                      num_knots = 20L),
                   S_cate_args = list(max_degree = 3L,
                                      smoothness_orders = 1L,
                                      num_knots = 20L,
                                      force_A = FALSE),
                   target_method = "tmle",
                   target_gwt = TRUE,
                   max_iter = 50,
                   n_lambda = 10,
                   g_bar_bound = NULL,
                   theta_bound = NULL,
                   Pi_bound = NULL,
                   Q_bar_bound = NULL,
                   parallel = FALSE,
                   verbose = TRUE,
                   browse = FALSE) {

      if (browse) browser()
      target_method <- match.arg(target_method, c("tmle", "relaxed"))

      self$controls_only <- all(self$A[self$S == 0] == 0)

      if (is.null(self$n_folds)) {
        n <- nrow(self$data)
        if (n <= 30) {
          self$n_folds <- n
        } else if (n <= 500) {
          self$n_folds <- 20
        } else if (n <= 1000) {
          self$n_folds <- 10
        } else if (n <= 10000) {
          self$n_folds <- 5
        } else {
          self$n_folds <- 3
        }
      }

      # cross fitting schemes --------------------------------------------------
      if (self$family == "binomial") {
        strata_ids <- paste0(self$S, self$Y)
      } else if (self$family == "gaussian") {
        strata_ids <- self$S
      }
      self$folds <- make_folds(n = nrow(self$data), V = self$n_folds,
                               strata_ids = strata_ids)
      foldid <- folds2foldvec(self$folds)
      self$foldsid <- foldid
      foldid_S1 <- foldid[self$S == 1]
      foldid_S0 <- foldid[self$S == 0]
      self$folds_S1 <- purrr::map(seq(self$n_folds), function(v) {
        fold_from_foldvec(v = v, folds = foldid_S1)
      })
      self$folds_S0 <- purrr::map(seq(self$n_folds), function(v) {
        fold_from_foldvec(v = v, folds = foldid_S0)
      })

      # initial estimation -----------------------------------------------------
      self$run_init_est(
        g_bar_method = g_bar_method,
        theta_method = theta_method,
        Pi_method = Pi_method,
        Q_bar_method = Q_bar_method,
        g_bar_bound = g_bar_bound,
        theta_bound = theta_bound,
        Pi_bound = Pi_bound,
        Q_bar_bound = Q_bar_bound
      )

      # obtain CATE working model ----------------------------------------------
      self$tau_A <- self$fit_cate(W = self$W,
                                  A = self$A,
                                  Y = self$Y,
                                  g1W = self$g_bar,
                                  theta = self$theta,
                                  cate_args = A_cate_args,
                                  parallel = parallel,
                                  force_A = FALSE)

      # obtain CARE working model ----------------------------------------------
      self$tau_S <- self$fit_cate(W = as.matrix(cbind(self$W, A = self$A)),
                                  A = self$S,
                                  Y = self$Y,
                                  g1W = self$Pi$A,
                                  theta = self$Q_bar$A,
                                  cate_args = S_cate_args,
                                  parallel = parallel,
                                  force_A = isTRUE(S_cate_args$force_A))
      self$tau_S$phi_W1 <- make_design_matrix(X = as.matrix(cbind(self$W, A = 1)),
                                              blist = self$tau_S$blist)
      self$tau_S$phi_W0 <- make_design_matrix(X = as.matrix(cbind(self$W, A = 0)),
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
      Pi_star_tmp <- purrr::map(self$tau_S_star, "Pi_star")
      self$Pi_star <- NULL

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
      self$Pi_star_avg_over_S1 <- purrr::map(self$tau_S_star_avg_over_S1, "Pi_star")
      self$Pi_star <- Pi_star_tmp

      # point estimate and inference -------------------------------------------
      self$inference()

      return(invisible(self))

    }
  ),
  private = list(
    get_tau_A_eic = function(phi_W,
                             cate_W,
                             avg_over_S1,
                             IM_inv = NULL) {
      phi_W <- as.matrix(phi_W)
      if (is.null(IM_inv)) {
        IM <- t(phi_W) %*% diag(self$g_bar * self$g_bar0) %*% phi_W / nrow(phi_W)
        IM_inv <- mat_inverse(IM)
      }

      if (avg_over_S1) {
        beta_comp <- as.vector(phi_W %*% IM_inv %*% colMeans((self$S / mean(self$S)) * phi_W) *
                                 (self$A - self$g_bar) *
                                 (self$Y - self$theta - (self$A - self$g_bar) * cate_W))
        W_comp <- self$S / mean(self$S) * (cate_W - mean(self$S / mean(self$S) * cate_W))
      } else {
        beta_comp <- as.vector(phi_W %*% IM_inv %*% colMeans(phi_W) *
                                 (self$A - self$g_bar) *
                                 (self$Y - self$theta - (self$A - self$g_bar) * cate_W))
        W_comp <- cate_W - mean(cate_W)
      }

      return(list(W_comp = W_comp,
                  beta_comp = beta_comp,
                  eic = W_comp + beta_comp))
    },

    get_tau_S_eic = function(phi_WA,
                             phi_W1,
                             phi_W0,
                             beta,
                             avg_over_S1,
                             IM_inv = NULL) {
      phi_WA <- as.matrix(phi_WA)
      phi_W1 <- as.matrix(phi_W1)
      phi_W0 <- as.matrix(phi_W0)
      if (is.null(IM_inv)) {
        IM <- t(phi_WA) %*% diag(self$Pi_star$A * (1 - self$Pi_star$A)) %*% phi_WA / nrow(phi_WA)
        IM_inv <- mat_inverse(IM)
      }

      tau_S <- list(
        phi_WA = phi_WA,
        phi_W1 = phi_W1,
        phi_W0 = phi_W0,
        cate_WA = as.numeric(phi_WA %*% beta),
        cate_W1 = as.numeric(phi_W1 %*% beta),
        cate_W0 = as.numeric(phi_W0 %*% beta)
      )
      Y_tmp <- self$Y
      Y_tmp[is.na(self$Y)] <- 0

      if (self$controls_only) {
        if (avg_over_S1) {
          psi_pound_est <- mean(self$S / mean(self$S) * (1 - self$Pi_star$A0) * tau_S$cate_W0)
          W_comp <- self$S / mean(self$S) * ((1 - self$Pi_star$A0) * tau_S$cate_W0 - psi_pound_est)
          Pi_comp <- -self$Pi_bar / mean(self$S) * (1 - self$A) / self$g_bar0 * tau_S$cate_W0 * (self$S - self$Pi_star$A)
        } else {
          psi_pound_est <- mean((1 - self$Pi_star$A0) * tau_S$cate_W0)
          W_comp <- (1 - self$Pi_star$A0) * tau_S$cate_W0 - psi_pound_est
          Pi_comp <- -(1 - self$A) / self$g_bar0 * tau_S$cate_W0 * (self$S - self$Pi_star$A)
        }
        D <- tau_S$phi_WA %*% IM_inv * (self$S - self$Pi_star$A) *
          (Y_tmp - self$Q_bar$A - (self$S - self$Pi_star$A) * tau_S$cate_WA) * self$weights
        if (avg_over_S1) {
          if (ncol(D) > 1) {
            beta_comp <- rowSums(D %*% diag(colMeans(self$S / mean(self$S) * (1 - self$Pi_star$A0) * tau_S$phi_W0)))
          } else {
            beta_comp <- rowSums(D * colMeans(self$S / mean(self$S) * (1 - self$Pi_star$A0) * tau_S$phi_W0))
          }
        } else {
          if (ncol(D) > 1) {
            beta_comp <- rowSums(D %*% diag(colMeans((1 - self$Pi_star$A0) * tau_S$phi_W0)))
          } else {
            beta_comp <- rowSums(D * colMeans((1 - self$Pi_star$A0) * tau_S$phi_W0))
          }
        }
      } else {
        if (avg_over_S1) {
          psi_pound_est <- mean(self$S / mean(self$S) *
                                  ((1 - self$Pi_star$A0) * tau_S$cate_W0 - (1 - self$Pi_star$A1) * tau_S$cate_W1))
          W_comp <- self$S / mean(self$S) *
            ((1 - self$Pi_star$A0) * tau_S$cate_W0 - (1 - self$Pi_star$A1) * tau_S$cate_W1 - psi_pound_est)
          Pi_comp <- self$Pi_bar / mean(self$S) *
            (self$A / self$g_bar * tau_S$cate_W1 - (1 - self$A) / self$g_bar0 * tau_S$cate_W0) *
            (self$S - self$Pi_star$A)
        } else {
          psi_pound_est <- mean((1 - self$Pi_star$A0) * tau_S$cate_W0 - (1 - self$Pi_star$A1) * tau_S$cate_W1)
          W_comp <- (1 - self$Pi_star$A0) * tau_S$cate_W0 - (1 - self$Pi_star$A1) * tau_S$cate_W1 - psi_pound_est
          Pi_comp <- (self$A / self$g_bar * tau_S$cate_W1 - (1 - self$A) / self$g_bar0 * tau_S$cate_W0) *
            (self$S - self$Pi_star$A)
        }
        D <- tau_S$phi_WA %*% IM_inv * (self$S - self$Pi_star$A) *
          (Y_tmp - self$Q_bar$A - (self$S - self$Pi_star$A) * tau_S$cate_WA) * self$weights
        if (avg_over_S1) {
          if (ncol(D) > 1) {
            beta_comp <- rowSums(D %*% diag(colMeans(self$S / mean(self$S) * (1 - self$Pi_star$A0) * tau_S$phi_W0))) -
              rowSums(D %*% diag(colMeans(self$S / mean(self$S) * (1 - self$Pi_star$A1) * tau_S$phi_W1)))
          } else {
            beta_comp <- rowSums(D * colMeans(self$S / mean(self$S) * (1 - self$Pi_star$A0) * tau_S$phi_W0)) -
              rowSums(D * colMeans(self$S / mean(self$S) * (1 - self$Pi_star$A1) * tau_S$phi_W1))
          }
        } else {
          if (ncol(D) > 1) {
            beta_comp <- rowSums(D %*% diag(colMeans((1 - self$Pi_star$A0) * tau_S$phi_W0))) -
              rowSums(D %*% diag(colMeans((1 - self$Pi_star$A1) * tau_S$phi_W1)))
          } else {
            beta_comp <- rowSums(D * colMeans((1 - self$Pi_star$A0) * tau_S$phi_W0)) -
              rowSums(D * colMeans((1 - self$Pi_star$A1) * tau_S$phi_W1))
          }
        }
      }

      return(list(W_comp = W_comp,
                  Pi_comp = Pi_comp,
                  beta_comp = beta_comp,
                  eic = W_comp + Pi_comp + beta_comp,
                  cate_WA = tau_S$cate_WA,
                  cate_W1 = tau_S$cate_W1,
                  cate_W0 = tau_S$cate_W0))
    }
  )
)
