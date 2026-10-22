#' TMLE with Pooled Outcome Regression
#'
#' R6 class implementing TMLE-PR using an initial outcome regression fitted on
#' all observations, including the binary source or subgroup indicator `S`.
#' Predictions evaluate this regression at `S = 1`. Both binary and continuous
#' outcomes are supported.
#'
#' @details
#' This ports the original TMLE-PR implementation from `atmleFusion`,
#' preserving its targeting, denominator truncation, inference, and result labels.
#'
#' Two targets are returned, in this order:
#' \itemize{
#'   \item `ATE (avg. over pooled)`: the S=1 treatment contrast averaged over
#'   all observed baseline covariates.
#'   \item `ATE (avg. over RCT)`: the S=1 treatment contrast averaged over
#'   baseline covariates among observations with S=1. For subgroup analysis,
#'   this is the subgroup-specific effect; the historical label is retained.
#' }
#'
#' @section Initialization:
#' Create an object using
#' `tmle_pr$new(data, S_node, W_nodes, A_node, Y_node, family,
#' n_folds = NULL, seed = 123)`.
#' `data` is a data frame with fully observed outcomes; `S_node` and `A_node`
#' identify numeric 0/1 columns, `W_nodes` identifies baseline covariates, and
#' `Y_node` identifies the outcome. `family` is `"binomial"` or `"gaussian"`,
#' or a corresponding family object. `n_folds` controls nuisance-learner
#' cross-validation; when omitted, the original sample-size-based rule is used.
#' `seed` sets the random seed when `run()` constructs folds.
#'
#' @section Fitting:
#' Call `fit$run(Q_method, g_method, Pi_bar_method,
#' g_bound = NULL, browse = FALSE)`.
#' Each method is an sl3 specification containing a `learners` list and, when
#' using multiple learners, a `metalearner`.
#' The outcome regression uses S, W, and A on the pooled sample; the treatment
#' regression uses W among S=1; the membership regression uses W on all rows.
#' Membership probabilities are needed for the pooled-covariate target.
#'
#' `g_bound` floors the products of the treatment probabilities and either the
#' membership probabilities (pooled target) or empirical S=1 proportion
#' (S=1 target). Its default is `5 / (sqrt(n) * log(n))`, where n is the total
#' sample size. This is the original implementation's product truncation rule.
#' `browse = TRUE` enters the R debugger before fitting.
#'
#' Known treatment probabilities can be supplied before `run()` by assigning
#' a length-n vector to `fit$g11W`; the treatment model is then skipped.
#' `fit$run_init_est()` fits the initial nuisance models, `fit$target_Q()` returns
#' targeted predictions for either population, and `fit$inference()` computes
#' results from the stored targeted predictions. The complete workflow is
#' provided by `run()`. After fitting, `fit$inference(alpha = 0.1)` computes
#' 90% Wald intervals; supply the same `g_bound` used for fitting if customized.
#'
#' @section Outputs:
#' `fit$results` contains `param`, `psi`, `lower`, `upper`, `se`, `alpha`,
#' `PnEIC`, and an `eic` list column. Standard errors use the sample variance
#' of the influence curve divided by the total sample size. `fit$eic` contains
#' the same curves as `pooled_W` and `rct_W`, each aligned with the input rows.
#' Initial predictions are in `fit$Q`; targeted predictions are in `fit$Q_star`
#' and `fit$Q_star_S1`.
#'
#' @md
#' @return An R6 class generator. Use `$new(...)` to construct an estimator.
#' @examples
#' set.seed(2026)
#' n <- 400
#' dat <- data.frame(W = rnorm(n), S = rbinom(n, 1, 0.4),
#'                   A = rbinom(n, 1, 0.5))
#' dat$Y <- with(dat, rbinom(n, 1, plogis(-1 + 0.4*W + 0.3*S + 0.5*A)))
#' method <- list(learners = list(sl3::Lrnr_glm$new()))
#' fit <- tmle_pr$new(dat, "S", "W", "A", "Y", "binomial", n_folds = 3)
#' fit$run(Q_method = method, g_method = method, Pi_bar_method = method)
#' fit$results[, c("param", "psi", "se", "lower", "upper")]
#'
#' @importFrom R6 R6Class
#' @importFrom origami make_folds folds2foldvec fold_from_foldvec
#' @importFrom purrr map
#' @importFrom sl3 sl3_Task
#' @importFrom stats coef glm offset plogis qlogis qnorm var weighted.mean
#' @export
tmle_pr <- R6::R6Class(
  classname = "TMLE for RCT + RWD",
  public = list(

    data = NULL,
    S_node = NULL,
    W_nodes = NULL,
    A_node = NULL,
    Y_node = NULL,
    family = NULL,
    n_folds = NULL,
    seed = NULL,
    S = NULL,
    W = NULL,
    A = NULL,
    Y = NULL,
    Delta = NULL,
    weights = NULL,
    folds = NULL,
    folds_S1 = NULL,
    Q = NULL,
    Q_star = NULL,
    Q_star_S1 = NULL,
    QSWA_fit = NULL,
    g_fit = NULL,
    g11W_fit = NULL,
    g11W = NULL,
    g01W = NULL,
    Pi_bar_fit = NULL,
    Pi_bar = NULL,
    eic = NULL,
    results = NULL,

    initialize = function(data,
                          S_node,
                          W_nodes,
                          A_node,
                          Y_node,
                          family,
                          n_folds = NULL,
                          seed = 123) {

      if (inherits(family, "family")) {
        family <- family$family
      }
      if (!is.character(family) ||
          length(family) != 1 ||
          !(family %in% c("gaussian", "binomial"))) {
        stop("`family` must be one of 'gaussian' or 'binomial'.")
      }

      self$data <- data
      self$S_node <- S_node
      self$W_nodes <- W_nodes
      self$A_node <- A_node
      self$Y_node <- Y_node
      self$family <- family
      self$n_folds <- n_folds
      self$seed <- seed
      self$S <- data[[S_node]]
      self$W <- data[, W_nodes, drop = FALSE]
      self$A <- data[[A_node]]
      self$Y <- data[[Y_node]]
      self$Delta <- rep(1, nrow(data))
      self$weights <- rep(1, nrow(data))
      self$Q <- NULL
      self$Q_star <- NULL
      self$Q_star_S1 <- NULL

    },

    target_Q = function(avg_over_S1, g_bound = NULL) {
      # scale to (0,1)
      if (self$family == "gaussian") {
        min_Y <- min(self$Y, self$Q$Q1WA, self$Q$Q1W1, self$Q$Q1W0)-0.001
        max_Y <- max(self$Y, self$Q$Q1WA, self$Q$Q1W1, self$Q$Q1W0)+0.001
        Y <- (self$Y-min_Y)/(max_Y-min_Y)
        Q1WA <- (self$Q$Q1WA-min_Y)/(max_Y-min_Y)
        Q1W1 <- (self$Q$Q1W1-min_Y)/(max_Y-min_Y)
        Q1W0 <- (self$Q$Q1W0-min_Y)/(max_Y-min_Y)
      } else if (self$family == "binomial") {
        Y <- self$Y; Q1WA <- self$Q$Q1WA; Q1W1 <- self$Q$Q1W1; Q1W0 <- self$Q$Q1W0
      } else {
        stop("Unsupported family. Use 'gaussian' or 'binomial'.")
      }

      # clever covariate
      if (is.null(g_bound)) {
        g_bound <- 5/sqrt(nrow(self$data))/log(nrow(self$data))
      }
      if (avg_over_S1) {
        denom_A1 <- .bound(mean(self$S) * self$g11W, c(g_bound, 1))
        denom_A0 <- .bound(mean(self$S) * self$g01W, c(g_bound, 1))
        HSAW <- self$S * (self$A/denom_A1-(1-self$A)/denom_A0)
        HS1W <- self$S/denom_A1
        HS0W <- -self$S/denom_A0
      } else {
        denom_A1 <- .bound(self$Pi_bar * self$g11W, c(g_bound, 1))
        denom_A0 <- .bound(self$Pi_bar * self$g01W, c(g_bound, 1))
        HSAW <- self$S * (self$A/denom_A1-(1-self$A)/denom_A0)
        HS1W <- 1/denom_A1
        HS0W <- -1/denom_A0
      }

      # logistic submodel
      epsilon <- as.numeric(
        coef(glm(Y ~ -1+offset(qlogis(Q1WA))+HSAW, family = "quasibinomial"))
      )
      epsilon[is.na(epsilon)] <- 0

      # TMLE updates
      Q_star <- list(Q1WA = plogis(qlogis(Q1WA)+epsilon*HSAW),
                     Q1W1 = plogis(qlogis(Q1W1)+epsilon*HS1W),
                     Q1W0 = plogis(qlogis(Q1W0)+epsilon*HS0W))

      # scale back
      if (self$family == "gaussian") {
        Q_star$Q1WA <- Q_star$Q1WA*(max_Y-min_Y)+min_Y
        Q_star$Q1W1 <- Q_star$Q1W1*(max_Y-min_Y)+min_Y
        Q_star$Q1W0 <- Q_star$Q1W0*(max_Y-min_Y)+min_Y
      }

      return(Q_star)

    },

    inference = function(alpha = 0.05,
                         g_bound = NULL) {
      # parameter that avg. over pooled W
      psi_pooled_W <- mean(self$Q_star$Q1W1-self$Q_star$Q1W0)
      eic_pooled_W <- get_np_eic_pooled_W(Q = self$Q_star,
                                          Pi = self$Pi_bar,
                                          g11W = self$g11W,
                                          S = self$S,
                                          A = self$A,
                                          Y = self$Y,
                                          psi = psi_pooled_W,
                                          g_bound = g_bound)
      pn_eic_pooled_W <- mean(eic_pooled_W, na.rm = TRUE)
      se_pooled_W <- sqrt(var(eic_pooled_W, na.rm = TRUE)/length(eic_pooled_W))
      lower_pooled_W <- psi_pooled_W+qnorm(alpha/2)*se_pooled_W
      upper_pooled_W <- psi_pooled_W+qnorm(1-alpha/2)*se_pooled_W

      # parameter that avg. over RCT W
      psi_rct_W_tmp <- mean(self$S/mean(self$S)*(self$Q_star_S1$Q1W1-self$Q_star_S1$Q1W0))
      psi_rct_W <- weighted.mean(self$Q_star_S1$Q1W1[self$S==1]-self$Q_star_S1$Q1W0[self$S==1],
                                 w = (self$S/mean(self$S))[self$S==1])
      eic_rct_W <- get_np_eic_rct_W(Q = self$Q_star_S1,
                                    pS = mean(self$S),
                                    g11W = self$g11W,
                                    S = self$S,
                                    A = self$A,
                                    Y = self$Y,
                                    psi = psi_rct_W,
                                    g_bound = g_bound)
      pn_eic_rct_W <- mean(eic_rct_W, na.rm = TRUE)
      se_rct_W <- sqrt(var(eic_rct_W, na.rm = TRUE)/length(eic_rct_W))
      lower_rct_W <- psi_rct_W+qnorm(alpha/2)*se_rct_W
      upper_rct_W <- psi_rct_W+qnorm(1-alpha/2)*se_rct_W

      self$results <- data.frame(param = c("ATE (avg. over pooled)", "ATE (avg. over RCT)"),
                                 psi = c(psi_pooled_W, psi_rct_W),
                                 lower = c(lower_pooled_W, lower_rct_W),
                                 upper = c(upper_pooled_W, upper_rct_W),
                                 se = c(se_pooled_W, se_rct_W),
                                 alpha = alpha,
                                 PnEIC = c(pn_eic_pooled_W, pn_eic_rct_W))
      self$results$eic <- I(list(eic_pooled_W, eic_rct_W))
      self$eic <- list(
        pooled_W = eic_pooled_W,
        rct_W = eic_rct_W
      )

      return(invisible(self$results))

    },

    run_init_est = function(Q_method,
                            g_method,
                            Pi_bar_method) {

      # Q(S,W,A)=E(Y|S,W,A)
      if (is.null(self$Q)) {
        self$QSWA_fit <- fit_regression(data = self$data,
                                        method = Q_method,
                                        folds = self$folds,
                                        covariate_nodes = c(self$S_node,
                                                            self$W_nodes,
                                                            self$A_node),
                                        outcome_node = self$Y_node)
        data_S1 <- self$data; data_S1[[self$S_node]] <- 1
        data_S1A1 <- data_S1; data_S1A1[[self$A_node]] <- 1
        data_S1A0 <- data_S1; data_S1A0[[self$A_node]] <- 0
        task_Q1WA <- sl3_Task$new(data = data_S1,
                                  covariates = c(self$S_node, self$W_nodes, self$A_node),
                                  outcome = self$Y_node,
                                  folds = self$folds)
        task_Q1W1 <- sl3_Task$new(data = data_S1A1,
                                  covariates = c(self$S_node, self$W_nodes, self$A_node),
                                  outcome = self$Y_node,
                                  folds = self$folds)
        task_Q1W0 <- sl3_Task$new(data = data_S1A0,
                                  covariates = c(self$S_node, self$W_nodes, self$A_node),
                                  outcome = self$Y_node,
                                  folds = self$folds)
        self$Q$Q1WA <- self$QSWA_fit$predict(task_Q1WA)
        self$Q$Q1W1 <- self$QSWA_fit$predict(task_Q1W1)
        self$Q$Q1W0 <- self$QSWA_fit$predict(task_Q1W0)
      }

      # g(1|S=1,W)=P(A=1|S=1,W)
      if (is.null(self$g11W)) {
        self$g_fit <- fit_regression(data = self$data,
                                     method = g_method,
                                     folds = self$folds_S1,
                                     covariate_nodes = self$W_nodes,
                                     outcome_node = self$A_node,
                                     subset = which(self$S == 1))
        task_g11W <- sl3_Task$new(data = self$data,
                                  covariates = self$W_nodes,
                                  outcome = self$A_node,
                                  folds = self$folds)
        self$g11W <- self$g_fit$predict(task_g11W)
      }
      self$g01W <- 1-self$g11W

      # Pi_bar(W)=P(S=1|W)
      if (is.null(self$Pi_bar)) {
        self$Pi_bar_fit <- fit_regression(data = self$data,
                                          method = Pi_bar_method,
                                          folds = self$folds,
                                          covariate_nodes = self$W_nodes,
                                          outcome_node = self$S_node)
        self$Pi_bar <- self$Pi_bar_fit$predict()
      }

      return(invisible(self))

    },

    run = function(Q_method,
                   g_method,
                   Pi_bar_method,
                   g_bound = NULL,
                   browse = FALSE) {

      if (browse) browser()

      if (is.null(self$n_folds)) {
        n_eff <- nrow(self$data)
        if (n_eff <= 30) {
          self$n_folds <- n_eff
        } else if (n_eff <= 500) {
          self$n_folds <- 20
        } else if (n_eff <= 1000) {
          self$n_folds <- 10
        } else if (n_eff <= 10000) {
          self$n_folds <- 5
        } else {
          self$n_folds <- 2
        }
      }

      # cross fitting schemes
      set.seed(self$seed)
      self$folds <- make_folds(n = nrow(self$data), V = self$n_folds,
                               strata_ids = self$S)
      foldid <- folds2foldvec(self$folds)
      foldid_S1 <- foldid[self$S == 1]
      self$folds_S1 <- map(seq(self$n_folds), function(v) {
        fold_from_foldvec(v = v, folds = foldid_S1)
      })

      # initial estimation -----------------------------------------------------
      self$run_init_est(Q_method = Q_method,
                        g_method = g_method,
                        Pi_bar_method = Pi_bar_method)

      # TMLE targeting of the outcome regression Q -----------------------------
      self$Q_star <- self$target_Q(avg_over_S1 = FALSE, g_bound = g_bound)

      # TMLE targeting of the outcome regression Q (avg. over S=1 only) --------
      self$Q_star_S1 <- self$target_Q(avg_over_S1 = TRUE, g_bound = g_bound)

      # point estimate and inference -------------------------------------------
      self$inference(g_bound = g_bound)

      return(invisible(self))

    }
  )
)
