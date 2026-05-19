#' TMLE Pooled Regression for RCT + RWD Data Fusion
#'
#' R6 class implementing a pooled-regression TMLE for average treatment effects
#' using randomized trial and real-world data.
#'
#' @details
#' The estimator fits a pooled outcome regression \eqn{E(Y \mid S, W, A)}, a
#' trial treatment mechanism \eqn{P(A = 1 \mid S = 1, W)}, and a trial
#' participation mechanism \eqn{P(S = 1 \mid W)}. It then targets the outcome
#' regression for two ATE parameters:
#'
#' \itemize{
#'   \item Average over the pooled covariate distribution.
#'   \item Average over the trial-only covariate distribution (\eqn{S=1}).
#' }
#'
#' @section Initialization:
#' Create the object with:
#'
#' `tmle_pr$new(data, S_node, W_nodes, A_node, Y_node, family, n_folds, seed)`
#'
#' @section Main Method (`$run()`):
#' The primary workflow is `fit$run(Q_method, g_method, Pi_bar_method, ...)`,
#' where each method is an `sl3` learner specification accepted by
#' `fit_regression()`.
#'
#' @return An R6 class generator object. Use `$new(...)` to instantiate.
#'
#' @importFrom origami make_folds folds2foldvec fold_from_foldvec
#' @importFrom purrr map
#' @importFrom R6 R6Class
#' @importFrom sl3 sl3_Task
#' @importFrom stats coef glm plogis qlogis qnorm var
#' @export
tmle_pr <- R6::R6Class(
  classname = "TMLE pooled regression for RCT + RWD",
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
      family <- match.arg(family, c("gaussian", "binomial"))

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
      if (self$family == "gaussian") {
        min_Y <- min(self$Y, self$Q$Q1WA, self$Q$Q1W1, self$Q$Q1W0) - 0.001
        max_Y <- max(self$Y, self$Q$Q1WA, self$Q$Q1W1, self$Q$Q1W0) + 0.001
        Y <- (self$Y - min_Y) / (max_Y - min_Y)
        Q1WA <- (self$Q$Q1WA - min_Y) / (max_Y - min_Y)
        Q1W1 <- (self$Q$Q1W1 - min_Y) / (max_Y - min_Y)
        Q1W0 <- (self$Q$Q1W0 - min_Y) / (max_Y - min_Y)
      } else {
        Y <- self$Y
        Q1WA <- .bound(self$Q$Q1WA, c(1e-6, 1 - 1e-6))
        Q1W1 <- .bound(self$Q$Q1W1, c(1e-6, 1 - 1e-6))
        Q1W0 <- .bound(self$Q$Q1W0, c(1e-6, 1 - 1e-6))
      }

      if (is.null(g_bound)) {
        g_bound <- 5 / sqrt(nrow(self$data)) / log(nrow(self$data))
      }
      if (avg_over_S1) {
        denom_A1 <- .bound(mean(self$S) * self$g11W, c(g_bound, 1))
        denom_A0 <- .bound(mean(self$S) * self$g01W, c(g_bound, 1))
        HSAW <- self$S * (self$A / denom_A1 - (1 - self$A) / denom_A0)
        HS1W <- self$S / denom_A1
        HS0W <- -self$S / denom_A0
      } else {
        denom_A1 <- .bound(self$Pi_bar * self$g11W, c(g_bound, 1))
        denom_A0 <- .bound(self$Pi_bar * self$g01W, c(g_bound, 1))
        HSAW <- self$S * (self$A / denom_A1 - (1 - self$A) / denom_A0)
        HS1W <- 1 / denom_A1
        HS0W <- -1 / denom_A0
      }

      epsilon <- as.numeric(coef(
        glm(Y ~ -1 + offset(qlogis(Q1WA)) + HSAW,
            family = stats::quasibinomial())
      ))
      epsilon[is.na(epsilon)] <- 0

      Q_star <- list(
        Q1WA = plogis(qlogis(Q1WA) + epsilon * HSAW),
        Q1W1 = plogis(qlogis(Q1W1) + epsilon * HS1W),
        Q1W0 = plogis(qlogis(Q1W0) + epsilon * HS0W)
      )

      if (self$family == "gaussian") {
        Q_star$Q1WA <- Q_star$Q1WA * (max_Y - min_Y) + min_Y
        Q_star$Q1W1 <- Q_star$Q1W1 * (max_Y - min_Y) + min_Y
        Q_star$Q1W0 <- Q_star$Q1W0 * (max_Y - min_Y) + min_Y
      }

      return(Q_star)
    },

    inference = function(alpha = 0.05,
                         g_bound = NULL) {
      psi_pooled_W <- mean(self$Q_star$Q1W1 - self$Q_star$Q1W0)
      eic_pooled_W <- tmle_pr_eic_pooled_W(
        Q = self$Q_star,
        Pi = self$Pi_bar,
        g11W = self$g11W,
        S = self$S,
        A = self$A,
        Y = self$Y,
        psi = psi_pooled_W,
        g_bound = g_bound
      )
      pn_eic_pooled_W <- mean(eic_pooled_W, na.rm = TRUE)
      se_pooled_W <- sqrt(var(eic_pooled_W, na.rm = TRUE) / length(eic_pooled_W))
      lower_pooled_W <- psi_pooled_W + qnorm(alpha / 2) * se_pooled_W
      upper_pooled_W <- psi_pooled_W + qnorm(1 - alpha / 2) * se_pooled_W

      psi_rct_W <- mean(
        self$Q_star_S1$Q1W1[self$S == 1] - self$Q_star_S1$Q1W0[self$S == 1]
      )
      eic_rct_W <- tmle_pr_eic_rct_W(
        Q = self$Q_star_S1,
        pS = mean(self$S),
        g11W = self$g11W,
        S = self$S,
        A = self$A,
        Y = self$Y,
        psi = psi_rct_W,
        g_bound = g_bound
      )
      pn_eic_rct_W <- mean(eic_rct_W, na.rm = TRUE)
      se_rct_W <- sqrt(var(eic_rct_W, na.rm = TRUE) / length(eic_rct_W))
      lower_rct_W <- psi_rct_W + qnorm(alpha / 2) * se_rct_W
      upper_rct_W <- psi_rct_W + qnorm(1 - alpha / 2) * se_rct_W

      self$results <- data.frame(
        param = c("ATE (avg. over pooled)", "ATE (avg. over RCT)"),
        psi = c(psi_pooled_W, psi_rct_W),
        lower = c(lower_pooled_W, lower_rct_W),
        upper = c(upper_pooled_W, upper_rct_W),
        se = c(se_pooled_W, se_rct_W),
        alpha = alpha,
        PnEIC = c(pn_eic_pooled_W, pn_eic_rct_W)
      )
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
      data_S1 <- self$data
      data_S1[[self$S_node]] <- 1
      data_S1A1 <- data_S1
      data_S1A1[[self$A_node]] <- 1
      data_S1A0 <- data_S1
      data_S1A0[[self$A_node]] <- 0

      if (is.null(self$Q)) {
        self$QSWA_fit <- fit_regression(
          data = self$data,
          method = Q_method,
          folds = self$folds,
          covariate_nodes = c(self$S_node, self$W_nodes, self$A_node),
          outcome_node = self$Y_node
        )
        task_Q1WA <- sl3_Task$new(
          data = data_S1,
          covariates = c(self$S_node, self$W_nodes, self$A_node),
          outcome = self$Y_node,
          folds = self$folds
        )
        task_Q1W1 <- sl3_Task$new(
          data = data_S1A1,
          covariates = c(self$S_node, self$W_nodes, self$A_node),
          outcome = self$Y_node,
          folds = self$folds
        )
        task_Q1W0 <- sl3_Task$new(
          data = data_S1A0,
          covariates = c(self$S_node, self$W_nodes, self$A_node),
          outcome = self$Y_node,
          folds = self$folds
        )
        self$Q$Q1WA <- as.numeric(self$QSWA_fit$predict(task_Q1WA))
        self$Q$Q1W1 <- as.numeric(self$QSWA_fit$predict(task_Q1W1))
        self$Q$Q1W0 <- as.numeric(self$QSWA_fit$predict(task_Q1W0))
      }

      if (is.null(self$g11W)) {
        self$g_fit <- fit_regression(
          data = self$data,
          method = g_method,
          folds = self$folds_S1,
          covariate_nodes = self$W_nodes,
          outcome_node = self$A_node,
          subset = which(self$S == 1)
        )
        task_g11W <- sl3_Task$new(
          data = self$data,
          covariates = self$W_nodes,
          outcome = self$A_node,
          folds = self$folds
        )
        self$g11W <- as.numeric(self$g_fit$predict(task_g11W))
      }
      self$g01W <- 1 - self$g11W

      if (is.null(self$Pi_bar)) {
        self$Pi_bar_fit <- fit_regression(
          data = self$data,
          method = Pi_bar_method,
          folds = self$folds,
          covariate_nodes = self$W_nodes,
          outcome_node = self$S_node
        )
        task_Pi_bar <- sl3_Task$new(
          data = self$data,
          covariates = self$W_nodes,
          outcome = self$S_node,
          folds = self$folds
        )
        self$Pi_bar <- as.numeric(self$Pi_bar_fit$predict(task_Pi_bar))
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

      set.seed(self$seed)
      self$folds <- make_folds(
        n = nrow(self$data),
        V = self$n_folds,
        strata_ids = self$S
      )
      foldid <- folds2foldvec(self$folds)
      foldid_S1 <- foldid[self$S == 1]
      self$folds_S1 <- map(seq(self$n_folds), function(v) {
        fold_from_foldvec(v = v, folds = foldid_S1)
      })

      self$run_init_est(
        Q_method = Q_method,
        g_method = g_method,
        Pi_bar_method = Pi_bar_method
      )

      self$Q_star <- self$target_Q(avg_over_S1 = FALSE, g_bound = g_bound)
      self$Q_star_S1 <- self$target_Q(avg_over_S1 = TRUE, g_bound = g_bound)
      self$inference(g_bound = g_bound)

      return(invisible(self))
    }
  )
)

tmle_pr_eic_pooled_W <- function(Q,
                                 Pi,
                                 g11W,
                                 S,
                                 A,
                                 Y,
                                 psi,
                                 g_bound = NULL) {
  if (is.null(g_bound)) {
    g_bound <- 5 / sqrt(length(Y)) / log(length(Y))
  }
  g01W <- 1 - g11W
  denom_A1 <- .bound(Pi * g11W, c(g_bound, 1))
  denom_A0 <- .bound(Pi * g01W, c(g_bound, 1))
  W_comp <- Q$Q1W1 - Q$Q1W0 - psi
  Q_comp <- S * (A / denom_A1 - (1 - A) / denom_A0) * (Y - Q$Q1WA)
  W_comp + Q_comp
}

tmle_pr_eic_rct_W <- function(Q,
                              pS,
                              g11W,
                              S,
                              A,
                              Y,
                              psi,
                              g_bound = NULL) {
  if (is.null(g_bound)) {
    g_bound <- 5 / sqrt(length(Y)) / log(length(Y))
  }
  g01W <- 1 - g11W
  denom_A1 <- .bound(pS * g11W, c(g_bound, 1))
  denom_A0 <- .bound(pS * g01W, c(g_bound, 1))
  W_comp <- (S / pS) * (Q$Q1W1 - Q$Q1W0 - psi)
  Q_comp <- S * (A / denom_A1 - (1 - A) / denom_A0) * (Y - Q$Q1WA)
  W_comp + Q_comp
}
