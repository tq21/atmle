#' @title TMLE R6 Class
#'
#' @import R6
#' @importFrom origami make_folds
tmle_R6 <- R6Class(
  "TMLE",
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
    foldid = NULL,
    g1W_bounds = NULL,
    g0W_bounds = NULL,
    family = NULL,
    Q_fit = NULL,
    Q = list(QAW = NA, Q0W = NA, Q1W = NA),
    Q_star = list(QAW = NA, Q0W = NA, Q1W = NA),
    g_fit = NULL,
    g1W = NULL,
    g1W_bd = NULL,
    g0W = NULL,
    g0W_bd = NULL,
    psi = NULL,
    eic = NULL,
    results = NULL,

    initialize = function(data,
                          W_nodes,
                          A_node,
                          Y_node,
                          family,
                          n_folds = NULL) {

      self$data <- data
      self$W_nodes <- W_nodes
      self$A_node <- A_node
      self$Y_node <- Y_node
      self$family <- family
      self$W <- data[, W_nodes, drop = FALSE]
      self$A <- data[[A_node]]
      self$Y <- data[[Y_node]]
      self$n_folds <- n_folds

    },

    create_folds = function(n_eff,
                            n_folds,
                            strata_ids,
                            seed) {

      set.seed(seed)

      # rule of thumb from tmle R package
      if (is.null(n_folds)) {
        if (n_eff <= 30){
          n_folds <- n_eff
        } else if (n_eff <= 500) {
          n_folds <- 20
        } else if (n_eff <= 1000) {
          n_folds <- 10
        } else if (n_eff <= 10000){
          n_folds <- 5
        } else {
          n_folds <- 2
        }
      }

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

    target = function(g1W_bounds,
                      g0W_bounds,
                      browse = FALSE) {

      if (browse) browser()

      # compute g-bound
      if (is.null(g1W_bounds)) {
        g1W_bounds <- c(5/sqrt(nrow(data))/log(nrow(data)), 1)
      }
      if (is.null(g0W_bounds)) {
        g0W_bounds <- c(5/sqrt(nrow(data))/log(nrow(data)), 1)
      }
      self$g0W <- 1-self$g1W
      self$g1W_bd <- .bound(self$g1W, g1W_bounds)
      self$g0W_bd <- .bound(self$g0W, g0W_bounds)

      # compute clever covariates
      HQ1W <- self$A/self$g1W_bd
      HQ0W <- (1-self$A)/self$g0W_bd
      HQAW <- HQ1W-HQ0W

      # bound Q
      if (self$family == "gaussian") {
        Y_bounds <- range(c(self$Y, self$Q$Q1W, self$Q$Q0W), na.rm = TRUE)+c(-0.001, 0.001)
        Y_bd <- (self$Y-Y_bounds[1])/(Y_bounds[2]-Y_bounds[1])
        QAW <- (self$Q$QAW-Y_bounds[1])/(Y_bounds[2]-Y_bounds[1])
        Q1W <- (self$Q$Q1W-Y_bounds[1])/(Y_bounds[2]-Y_bounds[1])
        Q0W <- (self$Q$Q0W-Y_bounds[1])/(Y_bounds[2]-Y_bounds[1])
      } else if (self$family == "binomial") {
        Y_bd <- self$Y
        QAW <- self$Q$QAW
        Q1W <- self$Q$Q1W
        Q0W <- self$Q$Q0W
      }

      # TMLE targeting
      fit_tmle <- glm(Y_bd ~ -1+offset(qlogis(QAW))+HQAW,
                      family = "quasibinomial")
      eps <- as.numeric(coef(fit_tmle))
      self$Q_star$QAW <- plogis(qlogis(QAW)+eps[1]*HQAW)
      self$Q_star$Q1W <- plogis(qlogis(Q1W)+eps[1]*HQ1W)
      self$Q_star$Q0W <- plogis(qlogis(Q0W)+eps[1]*HQ0W)

      if (self$family == "gaussian") {
        # transform back to original scale
        self$Q_star$QAW <- Y_bounds[1]+(Y_bounds[2]-Y_bounds[1])*self$Q_star$QAW
        self$Q_star$Q1W <- Y_bounds[1]+(Y_bounds[2]-Y_bounds[1])*self$Q_star$Q1W
        self$Q_star$Q0W <- Y_bounds[1]+(Y_bounds[2]-Y_bounds[1])*self$Q_star$Q0W
      }

    },

    run = function(Q_method,
                   g_method,
                   family,
                   g1W_bounds = NULL,
                   g0W_bounds = NULL,
                   seed = 123,
                   browse = FALSE) {

      if (browse) browser()

      # make folds -------------------------------------------------------------
      if (self$family == "binomial") {
        strata_ids <- paste0(self$A, "-", self$Y)
      } else {
        strata_ids <- paste0(self$A)
      }
      self$create_folds(n_eff = nrow(self$data),
                        n_folds = self$n_folds,
                        strata_ids = strata_ids,
                        seed = seed)

      # estimate outcome regression Q ------------------------------------------
      self$Q_fit <- self$fit_regression(method = Q_method,
                                        folds = self$folds,
                                        covariate_nodes = c(self$W_nodes,
                                                            self$A_node),
                                        outcome_node = self$Y_node)
      self$Q$QAW <- self$Q_fit$predict()

      # make counterfactual data and predictions
      data_A1 <- self$data; data_A1[[self$A_node]] <- 1
      data_A0 <- self$data; data_A0[[self$A_node]] <- 0
      Q1_task <- sl3_Task$new(data = data_A1,
                              covariates = c(self$W_nodes, self$A_node),
                              outcome = self$Y_node)
      Q0_task <- sl3_Task$new(data = data_A0,
                              covariates = c(self$W_nodes, self$A_node),
                              outcome = self$Y_node)
      self$Q$Q0W <- self$Q_fit$predict(Q0_task)
      self$Q$Q1W <- self$Q_fit$predict(Q1_task)

      # estimate treatment mechanism g -----------------------------------------
      self$g_fit <- self$fit_regression(method = g_method,
                                        folds = self$folds,
                                        covariate_nodes = self$W_nodes,
                                        outcome_node = self$A_node)
      self$g1W <- self$g_fit$predict()

      # TMLE targeting of the outcome regression Q -----------------------------
      self$target(g1W_bounds = g1W_bounds,
                  g0W_bounds = g0W_bounds)

      # point estimate and inference -------------------------------------------
      self$psi <- mean(self$Q_star$Q1W-self$Q_star$Q0W)
      self$inference()
    },

    inference = function(alpha = 0.05,
                         browse = FALSE) {

      if (browse) browser()
      self$psi <- mean(self$Q_star$Q1W-self$Q_star$Q0W)
      self$eic <- eic_ate(Y = self$Y,
                          A = self$A,
                          Q = self$Q_star,
                          g1W = self$g1W_bd,
                          g0W = self$g0W_bd,
                          psi = self$psi)
      se <- sqrt(var(self$eic, na.rm = TRUE)/nrow(self$data))
      lower <- self$psi+qnorm(alpha/2)*se
      upper <- self$psi+qnorm(1-alpha/2)*se

      self$results <- data.frame(psi = self$psi,
                                 lower = lower,
                                 upper = upper,
                                 se = se,
                                 PnEIC = mean(self$eic, na.rm = TRUE))

      return(invisible(self$results))

    }
  )
)
