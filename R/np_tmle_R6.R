library(R6)
np_tmle_R6 <- R6Class(
  classname = "TMLE for RCT + RWD",
  inherit = tmle_R6,
  public = list(

    data = NULL,
    S_node = NULL,
    W_nodes = NULL,
    A_node = NULL,
    Y_node = NULL,
    family = NULL,
    n_folds = NULL,
    S = NULL,
    W = NULL,
    A = NULL,
    Y = NULL,
    Delta = NULL,
    weights = NULL,
    folds = NULL,
    folds_S1 = NULL,
    Q_star_S1 = NULL,
    QSWA_fit = NULL,
    g11W_fit = NULL,
    g11W = NULL,
    g10W = NULL,
    Pi_bar_fit = NULL,
    Pi_bar = NULL,
    results = NULL,

    initialize = function(data,
                          S_node,
                          W_nodes,
                          A_node,
                          Y_node,
                          family,
                          n_folds = NULL,
                          seed = 123) {

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

    target_Q = function(avg_over_S1) {
      # scale to (0,1)
      if (self$family == "gaussian") {
        min_Y <- min(self$Y, self$Q$Q1WA, self$Q$Q1W1, self$Q$Q1W0)-0.001
        max_Y <- max(self$Y, self$Q$Q1WA, self$Q$Q1W1, self$Q$Q1W0)+0.001
        Y <- (self$Y-min_Y)/(max_Y-min_Y)
        Q1WA <- (self$Q$Q1WA-min_Y)/(max_Y-min_Y)
        Q1W1 <- (self$Q$Q1W1-min_Y)/(max_Y-min_Y)
        Q1W0 <- (self$Q$Q1W0-min_Y)/(max_Y-min_Y)
      } else {
        Y <- self$Y; Q1WA <- self$Q$Q1WA; Q1W1 <- self$Q$Q1W1; Q1W0 <- self$Q$Q1W0
      }

      # clever covariate
      if (avg_over_S1) {
        HSAW <- (self$S/mean(self$S))*(self$A/self$g11W-(1-self$A)/self$g10W)
        HS1W <- (self$S/mean(self$S))*(1/self$g11W)
        HS0W <- (self$S/mean(self$S))*(-1/self$g10W)
      } else {
        HSAW <- (self$S/self$Pi_bar)*(self$A/self$g11W-(1-self$A)/self$g10W)
        HS1W <- (self$S/self$Pi_bar)*(1/self$g11W)
        HS0W <- (self$S/self$Pi_bar)*(-1/self$g10W)
      }

      # logistic submodel
      epsilon <- coef(glm(Y ~ -1+offset(qlogis(Q1WA))+HSAW,
                          family = "quasibinomial"))
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

    inference = function(alpha = 0.05) {
      # parameter that avg. over pooled W
      psi_pooled_W <- mean(self$Q_star$Q1W1-self$Q_star$Q1W0)
      eic_pooled_W <- get_np_eic_pooled_W(Q = self$Q_star,
                                          Pi = self$Pi_bar,
                                          g11W = self$g11W,
                                          S = self$S,
                                          A = self$A,
                                          Y = self$Y,
                                          psi = psi_pooled_W)
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
                                    psi = psi_rct_W)
      se_rct_W <- sqrt(var(eic_rct_W, na.rm = TRUE)/length(eic_rct_W))
      lower_rct_W <- psi_rct_W+qnorm(alpha/2)*se_rct_W
      upper_rct_W <- psi_rct_W+qnorm(1-alpha/2)*se_rct_W

      self$results <- data.frame(param = c("ATE (avg. over pooled)", "ATE (avg. over RCT)"),
                                 psi = c(psi_pooled_W, psi_rct_W),
                                 lower = c(lower_pooled_W, lower_rct_W),
                                 upper = c(upper_pooled_W, upper_rct_W),
                                 se = c(se_pooled_W, se_rct_W),
                                 alpha = alpha)

      return(invisible(self$results))

    },

    run_init_est = function(Q_method,
                            g_method,
                            Pi_bar_method,
                            g_bound = NULL) {

      # Q(S,W,A)=E(Y|S,W,A)
      if (is.null(self$Q)) {
        self$QSWA_fit <- self$fit_regression(method = Q_method,
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
        self$g_fit <- self$fit_regression(method = g_method,
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
      if (is.null(g_bound)) {
        g_bound <- 5/sqrt(nrow(self$data))/log(nrow(self$data))
      }
      self$g11W <- .bound(self$g11W, c(g_bound, 1))
      self$g01W <- .bound(self$g01W, c(g_bound, 1))

      # Pi_bar(W)=P(S=1|W)
      if (is.null(self$Pi_bar)) {
        self$Pi_bar_fit <- self$fit_regression(method = Pi_bar_method,
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

      # cross fitting schemes
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
                        Pi_bar_method = Pi_bar_method,
                        g_bound = g_bound)

      # TMLE targeting of the outcome regression Q -----------------------------
      self$Q_star <- self$target_Q(avg_over_S1 = FALSE)

      # TMLE targeting of the outcome regression Q (avg. over S=1 only) --------
      self$Q_star_S1 <- self$target_Q(avg_over_S1 = TRUE)

      # point estimate and inference -------------------------------------------
      self$inference()

      return(self)

    }
  )
)
